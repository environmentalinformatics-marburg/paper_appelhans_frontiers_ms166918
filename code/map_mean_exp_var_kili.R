library(raster)
library(reshape2)
library(latticeExtra)
library(Rsenal)
library(remote)
library(gridExtra)
library(grid)

path_local <- "/media/ede/tims_ex/eot_kili_analysis"
setwd(path_local)

ind <- sprintf("%02.f", 0:12)

fls_list <- lapply(ind, function(j) {
  ptrn <- paste("*kili*", "lag", j, ".rda", sep = "")

  fls <- list.files("five_yr_chunks/results/modes",
                    pattern = glob2rx(ptrn),
                    full.names = TRUE)
  fls <- fls[-2]

})

### mean pixel importance
md <- "mode_01"
ext <- NULL
calcMeanPixRsq <- function(flist, mode = "mode_01") {

  rsq_all <- lapply(flist, function(i) {

    tmp <- do.call("stack", lapply(i, function(j) {
      print(j)
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

clr <- envinmrPalette
clr <- colorRampPalette(c("black", rev(brewer.pal(9, "YlGnBu")), "white"))

plotMeanPixRsq <- function(rst.list, mode = "mode_01", limit = 0.025,
                           ext = NULL) {

  p_lst <- lapply(seq(rst.list), function(i) {

    r <- rst.list[[i]]
    r[r <= limit] <- NA

    p <- spplot(r, col.regions = clr(1000),
                at = seq(0, 0.12, 0.001),
                colorkey = FALSE,
                panel = function(x, y, z, ...) {
                  panel.levelplot(x, y, z, ...)
                  panel.text(x = -165, y = 45, col = "white", font = 2,
                             labels = paste0(letters[i+1], ")"))
                })

    m <- layer(panel.map("world", border = "white", col = "grey60", lwd = 0.5))

    if (!is.null(ext)) {

      l <- layer(panel.extent(extent(ext), border = "darkred", lwd = 2))
      p + m + l
    } else p + m

  })

  return(p_lst)

}

p_lst <- plotMeanPixRsq(rsq_eot1[2:13], md, ext = ext, limit = 0)

p_112 <- latticeCombineGrid(p_lst, layout = c(3, 4))
p_0 <- spplot(rsq_eot1[[1]], col.regions = clr(1000),
              at = seq(0, 0.12, 0.001),
              main = "Mean coefficient of determination",
              colorkey = list(space = "top", width = 1, height = 0.75),
              panel = function(x, y, z, ...) {
                panel.levelplot(x, y, z, ...)
                panel.text(x = -160, y = 45, col = "white", font = 2,
                           labels = paste0(letters[1], ")"))
              }) +
  layer(panel.map("world", border = "white", col = "grey60", lwd = 0.5))


# png("five_yr_chunks/graphs/mean_pix_exp_var.png",
#      width = 250, height = 220, units = "mm", res = 600)
# grid.newpage()
#
# vp1 <- viewport(x = 0, y = 1, width = 1, height = 3/7,
#                 just = c("left", "top"), name = "topvp")
# pushViewport(vp1)
#
# print(p_0, newpage = FALSE)
#
# upViewport(0)
#
# vp2 <- viewport(x = 0, y = 0, width = 1, height = 4/7,
#                 just = c("left", "bottom"), name = "bottomvp")
# pushViewport(vp2)
#
# print(update(p_112, colorkey = list(draw = FALSE)), newpage = FALSE)
# dev.off()


tiff("five_yr_chunks/graphs/figure03.tif",
    width = 250, height = 220, units = "mm", res = 600,
    compression = "lzw")
grid.newpage()

vp1 <- viewport(x = 0, y = 1, width = 1, height = 3/7,
                just = c("left", "top"), name = "topvp")
pushViewport(vp1)

print(p_0, newpage = FALSE)

upViewport(0)

vp2 <- viewport(x = 0, y = 0, width = 1, height = 4/7,
                just = c("left", "bottom"), name = "bottomvp")
pushViewport(vp2)

print(update(p_112, colorkey = list(draw = FALSE)), newpage = FALSE)
dev.off()


### plot locations for each lag
locs <- lapply(seq(fls_list), function(i) {

  tmp <- do.call("rbind", lapply(seq(fls_list[[i]]), function(j) {

    load(fls_list[[i]][j])
    loc <- modes@modes[[md]]@coords_bp
    return(loc)

  }))

})

bg <- modes@modes$mode_01@r_predictor

plts <- lapply(seq(locs), function(i) {

  pts <- data.frame(locs[[i]])
  pts$chunk <- 1984:2007
  coordinates(pts) <- ~ x + y

  clrs <- colorRampPalette(brewer.pal(3, "YlOrBr"))

  spplot(bg, col.regions = "grey40", colorkey = FALSE) +
    layer(panel.map("world", col = "grey80", border = "white")) +
    as.layer(spplot(pts, col.regions = clrs(24),
                    cex = 2, edge.col = "black"))

})

latticeCombineGrid(plts[2:13], layout = c(3, 4))



### plot locations for each chunk
ind <- 1982:2005

fls_list <- lapply(ind, function(j) {
  ptrn <- paste("*kili*", j, "-*", ".rda", sep = "")

  fls <- list.files("five_yr_chunks/results/modes",
                    pattern = glob2rx(ptrn),
                    full.names = TRUE)
  fls <- fls[1:13]

})


locs <- lapply(seq(fls_list), function(i) {

  tmp <- do.call("rbind", lapply(seq(fls_list[[i]]), function(j) {

    load(fls_list[[i]][j])
    loc <- modes@modes[[md]]@coords_bp
    return(loc)

  }))

})

bg <- modes@modes$mode_01@r_predictor

plts <- lapply(seq(locs), function(i) {

  pts <- data.frame(locs[[i]])
  pts$lag <- 0:12
  coordinates(pts) <- ~ x + y

  clrs <- colorRampPalette(brewer.pal(3, "YlOrBr"))

  spplot(bg, col.regions = "grey40", colorkey = FALSE) +
    layer(panel.map("world", col = "grey80", border = "white")) +
    as.layer(spplot(pts, col.regions = clrs(13),
                    cex = 2, edge.col = "black"))

})

latticeCombineGrid(plts, layout = c(4, 6))


###########################################################################
# fls_nms <- list.files("five_yr_chunks/results/modes",
#                       pattern = glob2rx("*ewio*.rda"),
#                       full.names = FALSE)

## settings
md <- "mode_01"
ext <- extent(c(27.5, 45, -12.5, 5))

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

    p <- spplot(r, col.regions = clr(1000),
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

png("five_yr_chunks/graphs/mean_exp_var_kc_lag0.png", width = 30,
    height = 15, units = "cm", res = 300)
grid.newpage()
print(p_0)
dev.off()
###########################################################################


#mean_rsq_eot1 <- calc(rsqsums_eot1, fun = mean, na.rm = TRUE)

s <- mean_rsq_eot1
s[s < 0.02] <- NA

exv_p <- spplot(s, col.regions = envinmrPalette(200),
                main = "EWIO mean exp. var. per pixel")

ext_gpcp_ewio <- extent(c(27.5, 62.5, -17.5, 17.5))

exv_p + layer(panel.extent(ext = ext_gpcp_ewio))
