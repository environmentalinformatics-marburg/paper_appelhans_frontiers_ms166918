library(raster)
library(latticeExtra)
library(remote)
library(Rsenal)

#### global settings
data_path_local <- "/media/ede/tims_ex/eot_kili_analysis/"
# data_path_server <- "/media/memory01/data/casestudies/kilimanjaro_eot/"
setwd(data_path_local)
# setwd(data_path_server)
dir.create("tmp", showWarnings = FALSE) # for temporary raster files
#rasterOptions(maxmemory = 1e+06)
rasterOptions(tmpdir = "tmp")

ssns <- list("JF", "MAM", "JJAS", "OND")
lst_pos <- list(seq(1, 52, 4),
                seq(2, 52, 4),
                seq(3, 52, 4),
                seq(4, 52, 4))

out <- list()

for (i in seq(ssns)) {

  fls <- list.files("five_yr_chunks/results/seasonal/modes",
                    pattern = glob2rx(paste0("*", ssns[[i]], "*")),
                    full.names = TRUE)

  for (j in seq(fls)) {
    load(fls[j])
    tmp <- spplot(modes@modes$mode_02@rsq_predictor,
                  col.regions = envinmrPalette(1000),
                  at = seq(0, 1, length.out = 1000))

    out[[lst_pos[[i]][[j]]]] <- tmp
    #return(out)
  }
  #return(out)
}

png("five_yr_chunks/graphs/seasonal_modes_02.png",
    width = 30, height = 30, units = "cm", res = 300)
latticeCombineGrid(out, layout = c(4, 13))
dev.off()
