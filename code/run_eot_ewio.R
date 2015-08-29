library(raster)
library(latticeExtra)
library(remote)
# library(RghcnV3)

#### global settings
data_path_local <- "/media/ede/tims_ex/eot_kili_analysis/"
# data_path_server <- "/media/memory01/data/casestudies/kilimanjaro_eot/"
setwd(data_path_local)
# setwd(data_path_server)
dir.create("tmp", showWarnings = FALSE) # for temporary raster files
#rasterOptions(maxmemory = 1e+06)
rasterOptions(tmpdir = "tmp")

n_eots <- 3 # number of eots to calculate
chunk_size <- 5 # chunk size in years
frq <- 12 # frequency of data stacks
ext <- "ewio"

exp_var <- 80

rsp_ptrn <- paste("*", ext, exp_var, ".rda", sep = "")
prd_ptrn <- paste("*", "dns", exp_var, ".rda", sep = "")

### load data
rsp_dat <- list.files("five_yr_chunks/data", pattern = glob2rx(rsp_ptrn),
                      full.names = TRUE)
prd_dat <- list.files("five_yr_chunks/data", pattern = glob2rx(prd_ptrn),
                      full.names = TRUE)
load(rsp_dat) # kc
load(prd_dat) # sst

for (lag in seq.int(0, 12, 1)) {

  lngth <- frq - lag

  ## lagalize stacks
  dat_lst <- lagalize(x = stck_sst_dsn_dns, y = stck_gpcp_dsn_dns_ewio,
                      lag = lag, freq = frq)

  ## cut off overhead at the end to ensure similar length irrespective of lag
  dat_lst[[1]] <- dat_lst[[1]][[1:(nlayers(dat_lst[[1]]) - lngth)]]
  dat_lst[[2]] <- dat_lst[[2]][[1:(nlayers(dat_lst[[2]]) - lngth)]]
  names(dat_lst) <- c("sst", "gpcp")

  ## chunk start and end indeces
  st <- seq(1, nlayers(dat_lst$sst) - chunk_size * 12 + 1, 12)
  nd <- seq(chunk_size * 12, nlayers(dat_lst$gpcp), 12)

  ### now loop over chunks
  for (i in seq(st)) {

    prd <- dat_lst$sst[[st[i]:nd[i]]]
    rsp <- dat_lst$gpcp[[st[i]:nd[i]]]

    interval <- paste(substr(names(prd)[1], 5, 8),
                      substr(names(prd)[length(names(prd))], 5, 8),
                      sep = "-")

    lag_name <- paste("lag", sprintf("%02.f", lag), sep = "")

    out_name <- paste("sst_gpcp_dns", exp_var, "_", ext, "_",
                      interval, "_", lag_name, ".rda", sep = "")
    out_path <- paste("five_yr_chunks/results/modes", out_name, sep = "/")

    cat("\n",
        "*****************************************************************",
        "\n", "processing", ext, "chunk", interval, lag_name, "...", "\n",
        "*****************************************************************",
        "\n\n")

    modes <- eot(x = prd,
                 y = rsp,
                 n = n_eots,
                 standardised = TRUE,
                 reduce.both = TRUE)

    save(modes, file = out_path)

  }
  cat("finished lag", lag, "\n")
}

cat("finished all calculations...")
