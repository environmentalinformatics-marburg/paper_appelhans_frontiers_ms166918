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

### delete first year of chirps data for lagging
stck_prcp_kili_dsn_dns <- subset(stck_prcp_kili_dsn_dns,
                                 13:nlayers(stck_prcp_kili_dsn_dns))

seasons <- list(JF = c("01", "02"),
                MAM = c("03", "04", "05"),
                JJAS = c("06", "07", "08", "09"),
                OND = c("10", "11", "12"))

ssn_lag_lu <- list(JF = list(Jan = c(1, 12:1),
                             Feb = c(2:1, 12:2)),
                   MAM = list(Mar = c(3:1, 12:3),
                              Apr = c(4:1, 12:4),
                              May = c(5:1, 12:5)),
                   JJAS = list(Jun = c(6:1, 12:6),
                               Jul = c(7:1, 12:7),
                               Aug = c(8:1, 12:8),
                               Sep = c(9:1, 12:9)),
                   OND = list(Oct = c(10:1, 12:10),
                              Nov = c(11:1, 12:11),
                              Dec = c(12:1, 12)))

sst_nms <- sapply(seq(names(stck_sst_dsn_dns)), function(i) {
  strsplit(names(stck_sst_dsn_dns)[i], "_")[[1]][3]
})

chirps_nms <- sapply(seq(names(stck_prcp_kili_dsn_dns)), function(i) {
  strsplit(names(stck_prcp_kili_dsn_dns)[i], "\\.")[[1]][5]
})

chirps_ssns <- lapply(seq(seasons), function(i) {
  ind <- which(chirps_nms %in% seasons[[i]])
  subset(stck_prcp_kili_dsn_dns, ind)
})

sst_subs <- list(13:(nlayers(stck_sst_dsn_dns) - 0),
                 12:(nlayers(stck_sst_dsn_dns) - 1),
                 11:(nlayers(stck_sst_dsn_dns) - 2),
                 10:(nlayers(stck_sst_dsn_dns) - 3),
                 9:(nlayers(stck_sst_dsn_dns) - 4),
                 8:(nlayers(stck_sst_dsn_dns) - 5),
                 7:(nlayers(stck_sst_dsn_dns) - 6),
                 6:(nlayers(stck_sst_dsn_dns) - 7),
                 5:(nlayers(stck_sst_dsn_dns) - 8),
                 4:(nlayers(stck_sst_dsn_dns) - 9),
                 3:(nlayers(stck_sst_dsn_dns) - 10),
                 2:(nlayers(stck_sst_dsn_dns) - 11),
                 1:(nlayers(stck_sst_dsn_dns) - 12))

for (season in seq(seasons)) {

  rsp <- chirps_ssns[[season]]

  for (lag in seq.int(0, 12, 1)) {

    lag_name <- paste("lag", sprintf("%02.f", lag), sep = "")

    out_name <- paste("sst_prcp_dsn_dns", exp_var, "_", ext, "_",
                      names(seasons)[season], "_",
                      lag_name, ".rda", sep = "")
    out_path <- paste("five_yr_chunks/results/seasonal/modes",
                      out_name, sep = "/")

    cat("\n",
        "*****************************************************************",
        "\n", "processing", ext, names(seasons)[season], lag_name, "...", "\n",
        "*****************************************************************",
        "\n\n")

    pred <- subset(stck_sst_dsn_dns, sst_subs[[lag + 1]])

    ind <- sapply(seq(ssn_lag_lu[[season]]), function(i) {
      ssn_lag_lu[[season]][[i]][lag + 1]
    })

    pred_nms <- sapply(seq(names(pred)), function(i) {
      strsplit(names(pred)[i], "_")[[1]][3]
    })

    indx <- which(pred_nms %in% sprintf("%02.f", ind))
    prd <- subset(pred, indx)

    modes <- eot(x = prd,
                 y = rsp,
                 n = n_eots,
                 standardised = TRUE,
                 reduce.both = TRUE)

    save(modes, file = out_path)

  }
  cat("finished season", names(seasons)[season], "\n")
}

cat("finished all calculations...")
