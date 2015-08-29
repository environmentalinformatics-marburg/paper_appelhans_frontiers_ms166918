library(raster)
library(reshape2)
library(latticeExtra)
library(Rsenal)
library(remote)
library(gridExtra)

path_local <- "/media/ede/tims_ex/eot_kili_analysis"
setwd(path_local)

ind <- sprintf("%02.f", 0:12)

fls_list <- lapply(ind, function(j) {
  ptrn <- paste("*ewio*", "lag", j, ".rda", sep = "")

  fls <- list.files("five_yr_chunks/results/modes",
                    pattern = glob2rx(ptrn),
                    full.names = TRUE)

})

# fls_nms <- list.files("five_yr_chunks/results/modes",
#                       pattern = glob2rx("*ewio*.rda"),
#                       full.names = FALSE)

## settings
md <- "mode_01"
ext <- extent(c(27.5, 62.5, -17.5, 17.5))

calcMeanPixRsq <- function(flist, mode = "mode_01") {

  rsq_all <- lapply(flist, function(i) {

    tmp <- do.call("stack", lapply(i, function(j) {

      load(j)
      rsq <- modes@modes[[mode]]@rsq_sums_predictor /
        ncell(modes@modes[[mode]]@r_response)
      return(rsq)
      rm(modes)
    }))

    calc(tmp, fun = mean, na.rm = TRUE)
  })

  return(rsq_all)

}

rsq_eot1 <- calcMeanPixRsq(fls_list, md)

plotMeanPixRsq <- function(rst.list, mode = "mode_01", limit = 0.025,
                           ext = NULL) {

  p_lst <- lapply(seq(rst.list), function(i) {

    r <- rst.list[[i]]
    r[r <= limit] <- NA

    p <- spplot(r, col.regions = envinmrPalette(1000),
                at = seq(0, 0.1, 0.001),
                panel = function(x, y, z, ...) {
                  panel.levelplot(x, y, z, ...)
                  panel.text(x = 80, y = 55, labels = paste(mode, " ", "lag ",
                                                            ind[i], sep = ""))
                })

    m <- layer(Rsenal:::panel.map("world", border = "grey70"))

    if (!is.null(ext)) {

      l <- layer(panel.extent(extent(ext), border = "darkred", lwd = 2))
      p + m + l
    } else p + m

  })

  return(p_lst)

}

p_lst <- plotMeanPixRsq(rsq_eot1, md, ext = ext)

p_112 <- latticeCombineGrid(p_lst[2:13], layout = c(3, 4))
p_0 <- p_lst[[1]]

png("five_yr_chunks/graphs/mean_exp_var_ewio_lags.png", width = 30,
    height = 15, units = "cm", res = 300)
grid.newpage()
print(p_112)
dev.off()
###########################################################################




#mean_rsq_eot1 <- calc(rsqsums_eot1, fun = mean, na.rm = TRUE)

s <- mean_rsq_eot1
s[s < 0.02] <- NA

exv_p <- spplot(s, col.regions = envinmrPalette(200),
                main = "EWIO mean exp. var. per pixel")

ext_gpcp_ewio <- extent(c(27.5, 62.5, -17.5, 17.5))

exv_p + layer(panel.extent(ext = ext_gpcp_ewio))
