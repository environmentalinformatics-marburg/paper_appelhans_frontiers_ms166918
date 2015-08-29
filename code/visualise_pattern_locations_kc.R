library(raster)
library(latticeExtra)
library(Rsenal)
library(mapdata)
library(maptools)
library(MASS)
library(gridExtra)

path <- "/media/ede/tims_ex/eot_kili_example_data"
setwd(path)

### visualise respective pattern locations
domain <- "kc"
eot <- "EOT_01"

patterns <- list(p1 = list(lags = c("lag00", "lag01", "lag02",
                                    "lag03", "lag04", "lag05"),
                           chunks = c("1993-1997", "1994-1998",
                                      "1995-1999", "1996-2000",
                                      "1997-2001"),
                           pal = "Reds"),
                 p2 = list(lags = c("lag08", "lag09", "lag10"),
                           chunks = c("1993-1997", "1994-1998",
                                      "1995-1999", "1996-2000",
                                      "1997-2001"),
                           pal = "Purples"),
                 p3 = list(lags = c("lag00", "lag01", "lag02"),
                           chunks = c("2002-2006", "2003-2007", "2004-2008"),
                           pal = "Oranges"),
                 p4 = list(lags = c("lag10"),
                           chunks = c("2002-2006", "2003-2007",
                                      "2004-2008", "2005-2009"),
                           pal = "Greens"))



dat_patterns <- lapply(patterns, function(p) {
  fls <- list.files("five_yr_chunks/results/modes",
                    pattern = glob2rx("*kc*.rda"), full.names = TRUE,
                    recursive = TRUE)

  fls_lags <- do.call("c", lapply(p$lags, function(i) {
    grep(i, fls, value = TRUE)
  }))

  fls_lags_chunks <- do.call("c", lapply(p$chunks, function(i) {
    grep(i, fls_lags, value = TRUE)
  }))

  tst <- do.call("rbind", lapply(seq(fls_lags_chunks), function(i) {
    current <- fls_lags_chunks[i]

    load(current)

    df <- data.frame(x = modes@modes$mode_01@coords_bp[, 1],
                     y = modes@modes$mode_01@coords_bp[, 2],
                     eot = "EOT_01")
    rownames(df) <- NULL
    return(df)
  }))

  return(tst)

})



locs.lst <- dat_patterns[[domain]]

wm <- map("world", plot = FALSE, fill = TRUE)

dens.rst <- raster(nrows = 360, ncols = 120,
                   xmn = -60, xmx = 60, ymn = -180, ymx = 180)
dens.rst[] <- rnorm(360 * 120, 2, 1)

#dens.gr <- matrix(dens$z, ncol = length(dens$x), nrow = length(dens$y))

clrs <- colorRampPalette(c("white", rev(brewer.pal(9, "Spectral"))))
levs <- seq(min(dens.rst[]), max(dens.rst[]), length.out = 999)

p.locs <- spplot(t(dens.rst), mm = wm, at = levs, #main = png.name,
                 colorkey = FALSE, #list(space = "top", height = 0.75, width = 1),
                 col.regions = "grey30", panel = function(..., mm) {
                   panel.levelplot(...)
                   panel.polygon(mm$x, mm$y, lwd = 0.5,
                                 border = "grey20", col = "grey80")
                   #                    panel.rect(ext.gpcp.ls[[domain.nr]]@xmin,
                   #                               ext.gpcp.ls[[domain.nr]]@ymin,
                   #                               ext.gpcp.ls[[domain.nr]]@xmax,
                   #                               ext.gpcp.ls[[domain.nr]]@ymax,
                   #                               border = "black", lwd = 1.5, lty = 1)
                 })

png.name <- paste("overview_patterns_eot1_", domain, ".png", sep = "")

png(paste("graphs", png.name, sep = "/"),
    width = 20, height = 10, units = "cm", res = 300)
grid.newpage()
print(p.locs)
downViewport(trellis.vpname(name = "figure"))

for (p in seq(patterns)) {
  eot1 <- dat_patterns[[p]]

  clr <- colorRampPalette(brewer.pal(9, patterns[[p]]$pal))

  for (k in 1:nrow(eot1)) {

    xpos <- eot1$x[k] / 360 + 0.5
    ypos <- eot1$y[k] / 120 + 0.5

    vp <- viewport(x = xpos, y = ypos,
                   height = 0.05, width = 0.05,
                   just = c("centre", "centre"))

    pushViewport(vp)

    grid.circle(gp = gpar(fill = clr(nrow(eot1))[k]))
    #grid.text(label = eot1$chunk[k], gp = gpar(cex = 0.3))

    upViewport()

  }

}

dev.off()
