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
ext <- "kili"
exp_var <- 80

rsp_ptrn <- paste("*", "chirps", "*", ext, exp_var, ".rda", sep = "")
prd_ptrn <- paste("*", "dns", exp_var, ".rda", sep = "")

### load data
rsp_dat <- list.files("five_yr_chunks/data", pattern = glob2rx(rsp_ptrn),
                      full.names = TRUE)
prd_dat <- list.files("five_yr_chunks/data", pattern = glob2rx(prd_ptrn),
                      full.names = TRUE)
load(rsp_dat) # kc
load(prd_dat) # sst

for (lag in seq.int(0, 12, 1)) {

  lag_lst <- lagalize(x = stck_sst_dsn_dns,
                      y = stck_prcp_kili_dsn_dns,
                      lag = lag,
                      freq = frq)

  prd <- lag_lst[[1]]
  rsp <- lag_lst[[2]]

  lag_name <- paste("lag", sprintf("%02.f", lag), sep = "")

  out_name <- paste("sst_prcp_dsn_dns", exp_var, "_", ext, "_",
                    "1982-2010", "_", lag_name, ".rda", sep = "")
  out_path <- paste("five_yr_chunks/results/modes", out_name, sep = "/")

  cat("\n",
      "*****************************************************************",
      "\n", "processing", ext, lag_name, "...", "\n",
      "*****************************************************************",
      "\n\n")

  modes <- eot(x = prd,
               y = rsp,
               n = n_eots,
               standardised = TRUE,
               reduce.both = TRUE)

  save(modes, file = out_path)

}

cat("finished all calculations...")
