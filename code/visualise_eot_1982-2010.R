library(raster)
library(remote)
library(reshape2)
library(latticeExtra)
library(Rsenal)

path_local <- "/media/ede/tims_ex/eot_kili_analysis"
setwd(path_local)

fls <- list.files("five_yr_chunks/results/modes",
                  pattern = glob2rx("*kili_1982-2010*.rda"),
                  full.names = TRUE)


### mode 1
stck1 <- do.call(stack, lapply(fls, function(i) {
  load(i)
  rst <- modes@modes$mode_01@rsq_predictor
  return(rst)
  rm(modes)
}))

plt1 <- lapply(seq(nlayers(stck1)), function(i) {
  spplot(stck1[[i]], col.regions = envinmrPalette(1000),
         at = seq(0, 1, length.out = 1000))
})

# latticeCombineGrid(plt[2:13])


### mode2
stck2 <- do.call(stack, lapply(fls, function(i) {
  load(i)
  rst <- modes@modes$mode_02@rsq_predictor
  return(rst)
  rm(modes)
}))

plt2 <- lapply(seq(nlayers(stck2)), function(i) {
  spplot(stck2[[i]], col.regions = envinmrPalette(1000),
         at = seq(0, 1, length.out = 1000))
})

# latticeCombineGrid(plt[2:13])


### mode3
stck3 <- do.call(stack, lapply(fls, function(i) {
  load(i)
  rst <- modes@modes$mode_03@rsq_predictor
  return(rst)
  rm(modes)
}))

plt3 <- lapply(seq(nlayers(stck3)), function(i) {
  spplot(stck3[[i]], col.regions = envinmrPalette(1000),
         at = seq(0, 1, length.out = 1000))
})

# latticeCombineGrid(plt[2:13])


plt_lst <- lapply(seq(plt1), function(i) {
  c(plt1[[i]], plt2[[i]], plt3[[i]])
})

latticeCombineGrid(plt_lst, layout = c(3, 13))
