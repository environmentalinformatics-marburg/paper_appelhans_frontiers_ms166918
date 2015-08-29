library(raster)
library(latticeExtra)
library(remote)
library(Rsenal)
# library(RghcnV3)

#### global settings
data_path_local <- "/media/ede/tims_ex/eot_kili_analysis/"
# data_path_server <- "/media/memory01/data/casestudies/kilimanjaro_eot/"
setwd(data_path_local)
# setwd(data_path_server)

#### create input data stacks #############################################
### read file names
exp_var <- 0.8

files_sst <- list.files("source_data/sst_1982_2010",
                        pattern = glob2rx("*.rst"),
                        full.names = TRUE)

# files_prcp <- list.files("source_data/gpcp_1982_2010",
#                          pattern = glob2rx("*.rst"),
#                          full.names = TRUE)

files_prcp <- list.files("source_data/chirps_1982_2010",
                         pattern = glob2rx("*.tif"),
                         full.names = TRUE)

stck_sst <- stack(files_sst)
nms_sst <- names(stck_sst)
stck_prcp <- stack(files_prcp)
nms_prcp <- names(stck_prcp)

### get kili extent from hemp dem
dem <- raster("../kiliDEM/in/dem_hemp_ll_wgs84.tif")
ext <- extent(dem)
stck_prcp_kili <- crop(stck_prcp, ext, snap = "out")

stck_sst_dsn <- deseason(stck_sst, cycle.window = 12)
stck_prcp_dsn <- deseason(stck_prcp, cycle.window = 12)
stck_prcp_kili_dsn <- deseason(stck_prcp_kili, cycle.window = 12)

stck_sst_dsn_dns <- denoise(stck_sst_dsn, expl.var = exp_var)
stck_prcp_dsn_dns <- denoise(stck_prcp_dsn, expl.var = exp_var)
stck_prcp_kili_dsn_dns <- denoise(stck_prcp_kili_dsn, exp_var)

load("source_data/lmsk1Deg.RData")
lm <- lmsk1Deg

stck_sst_dsn_dns <- stck_sst_dsn_dns * lm
names(stck_sst_dsn_dns) <- nms_sst

ext_sst <- extent(c(-180, 180, -60, 60))
ext_prcp_kc <- extent(c(27.5, 45, -12.5, 5))
ext_prcp_ewio <- extent(c(27.5, 62.5, -17.5, 17.5))

stck_sst_dsn_dns <- crop(stck_sst_dsn_dns, ext_sst)
stck_prcp_dsn_dns_ewio <- crop(stck_prcp_dsn_dns, ext_prcp_ewio)
stck_prcp_dsn_dns_kc <- crop(stck_prcp_dsn_dns, ext_prcp_kc)

save(stck_sst_dsn_dns,
     file = paste("five_yr_chunks/data/stck_sst_dsn_dns",
                  exp_var * 100, ".rda", sep = ""))
# save(stck_prcp_dsn_dns_ewio,
#      file = paste("five_yr_chunks/data/stck_gpcp_dsn_dns_ewio",
#                   exp_var * 100, ".rda", sep = ""))
# save(stck_prcp_dsn_dns_kc,
#      file = paste("five_yr_chunks/data/stck_gpcp_dsn_dns_kc",
#                   exp_var * 100, ".rda", sep = ""))

save(stck_prcp_dsn_dns_ewio,
     file = paste("five_yr_chunks/data/stck_chirps_dsn_dns_ewio",
                  exp_var * 100, ".rda", sep = ""))
save(stck_prcp_dsn_dns_kc,
     file = paste("five_yr_chunks/data/stck_chirps_dsn_dns_kc",
                  exp_var * 100, ".rda", sep = ""))

save(stck_prcp_kili_dsn_dns,
     file = paste("five_yr_chunks/data/stck_chirps_dsn_dns_kili",
                  exp_var * 100, ".rda", sep = ""))
###########################################################################
