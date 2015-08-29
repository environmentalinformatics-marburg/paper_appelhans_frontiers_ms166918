library("rworldmap")
library("rgdal")
library("rgeos")
library("gridExtra")
library("RColorBrewer")
library("oce")
library("mapproj")
library("raster")
library("OpenStreetMap")
library("Rsenal")
library("latticeExtra")


path <- "/media/ede/tims_ex/eot_kili_analysis"
setwd(path)

### sst data
files_sst <- list.files("source_data/sst_1982_2010",
                        pattern = glob2rx("*.rst"),
                        full.names = TRUE)
load("source_data/lmsk1Deg.RData")
lm <- lmsk1Deg

sst_stck <- stack(files_sst)
sst_stck <- sst_stck * lm
ext <- extent(c(-180, 180, -60, 60))
sst_mean <- crop(calc(sst_stck, mean, na.rm = TRUE), ext)

tza <- getData("GADM", country = "TZA", level = 0)

sst_p <- spplot(sst_mean, col.regions = envinmrPalette(1000),
                at = seq(-2.5, 32.5, 0.1), main = "Temperature [°C]",
                scales = list(alternating = 3,
                              x = list(at = seq(-180, 180, 60)),
                              y = list(at = seq(-60, 60, 30))),
                colorkey = list(space = "top", height = 0.5, width = 1)) +
  layer(panel.map("world", border = "white", col = "grey70", lwd = 0.5)) +
  as.layer(spplot(tza, zcol = "ID_0",
                  col.regions = "grey50", col = "transparent"))




### chirps data
files_chirps <- list.files("source_data/chirps_1982_2010",
                           pattern = glob2rx("*.tif"),
                           full.names = TRUE)
chirps_stck <- stack(files_chirps)
dem <- raster("../kiliDEM/in/dem_hemp_ll_wgs84.tif")
ext_dem <- extent(dem)
chirps_stck_kili <- crop(chirps_stck, ext_dem, snap = "out")
#extent(dem) <- extent(chirps_stck_kili)

indx_st <- seq(1, nlayers(chirps_stck), 12)
indx_nd <- seq(12, nlayers(chirps_stck), 12)

chirps_mean <- calc(stack(lapply(seq(indx_st), function(i) {
  tmp <- subset(chirps_stck_kili, indx_st[i]:indx_nd[i])
  tmp_sum <- calc(tmp, sum, na.rm = TRUE)
  return(tmp_sum)
})), mean, na.rm = TRUE)

dem_flipped <- flip(dem, "y")
x_dem <- coordinates(dem_flipped)[, 1]
y_dem <- coordinates(dem_flipped)[, 2]
z_dem <- dem_flipped[]

clrs <- colorRampPalette(brewer.pal(9, "YlGnBu")[2:9])
clrs_conts <- colorRampPalette(rev(brewer.pal(9, "Greys")[4:9]))

prcp <- spplot(chirps_mean, col.regions = clrs(1000),
               main = "Precipitation [mm/a]",
               colorkey = list(space = "top", height = 0.75, width = 1),
               scales = list(draw = TRUE, alternating = 3),
               at = seq(300, 2600, 10), panel = function(x, y, z, ...) {
                 panel.levelplot(x, y, z, ...)
                 Rsenal:::panel.filledcontour(x_dem, y_dem, z_dem,
                                              at.contour = seq(500, 6000, 500),
                                              fill.contours = FALSE,
                                              col.contours = clrs_conts(15))
               })

# cp <- levelplot(z ~ x + y, col.regions = "transparent",
#                 panel = Rsenal:::panel.filledcontour,
#                 at.contour = seq(500, 6000, 500),
#                 fill.contours = FALSE, col.contours = "grey60")

sst_p_kili <- sst_p +
  as.layer(xyplot(mean(coordinates(dem_flipped)[, 2], na.rm = TRUE) ~
                    mean(coordinates(dem_flipped)[, 1], na.rm = TRUE),
                  pch = 24, cex = 1, fill = "white", col = "black"))

sst_p_kili

png("figure01.png", width = 27, height = 15, units = "cm",
     res = 600)
plot.new()
print(sst_p_kili) #+ layer(panel.extent(extent(dem)))
dev.off()

png("figure02.png", width = 20, height = 20, units = "cm",
     res = 600)
plot.new()
print(prcp) #+ layer(panel.extent(extent(dem)))
dev.off()

# tiff("figure01.tif", width = 27, height = 15, units = "cm",
#      res = 300, compression = "lzw")
# plot.new()
# print(sst_p_kili) #+ layer(panel.extent(extent(dem)))
# dev.off()
#
# tiff("figure02.tif", width = 20, height = 20, units = "cm",
#      res = 300, compression = "lzw")
# plot.new()
# print(prcp) #+ layer(panel.extent(extent(dem)))
# dev.off()




###########################################################################





# ### projection & extent of fulldisk
# p_sst <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +over"
# #p_sst <- "+proj=eck4 +over"
# ext_glob <- extent(c(-180, 180, -90, 90))
# ext_sst <- extent(-35, 125, -80, 80)
#
#
#
# ### get bing data global
# glob <- openmap(upperLeft = c(80, -180), lowerRight = c(-80, 180),
#                 minNumTiles = 10, type = "mapquest-aerial")
# #glob_rob <- openproj(glob, projection = p_sst)
# glob_rst <- raster(glob)
#
# # data(wrld_simpl)
# # world_ext = projectExtent(wrld_simpl, crs = p_sst)
# # glob_rst_crop = crop(x = glob_rst, y = world_ext, snap = 'out')
# #
# # glob_rst_crop_longlat = projectRaster(glob_rst_crop, crs = '+proj=longlat')
# # glob_rst_crop_rob <- projectRaster(glob_rst_crop_longlat,
# #                                    crs = p_sst, over = TRUE)
# # glob_rst_crop_rob_tr <- trim(glob_rst_crop_rob)
#
#
# lout <- rgb2spLayout(glob_rst)
# pt <- data.frame(ID = "point",
#                  x = mean(coordinates(glob_rst)[, 1]),
#                  y = mean(coordinates(glob_rst)[, 2]))
# coordinates(pt) <- ~x+y
# projection(pt) <- CRS(projection(glob_rst))
# pt@bbox[1, 1] <- extent(glob_rst)@xmin
# pt@bbox[1, 2] <- extent(glob_rst)@xmax
# pt@bbox[2, 1] <- extent(glob_rst)@ymin
# pt@bbox[2, 2] <- extent(glob_rst)@ymax
# omap_plt <- spplot(pt, sp.layout = lout)
#
# ### sst data
# files_sst <- list.files("/media/ede/tims_ex/eot_kili_analysis/source_data/sst_1982_2010",
#                         pattern = glob2rx("*.rst"),
#                         full.names = TRUE)
# load("/media/ede/tims_ex/eot_kili_analysis/source_data/lmsk1Deg.RData")
# lm <- lmsk1Deg
#
# sst_stck <- stack(files_sst)
# sst_stck <- sst_stck * lm
# sst_stck <- crop(sst_stck, extent(c(-180, 180, -60, 60)))
#
# sst_mean <- calc(sst_stck, mean, na.rm = TRUE)
#
# sst_mean_rob <- projectRaster(sst_mean, crs = projection(glob_rst))
# # sst_mean_rob_tr <- trim(sst_mean_rob)
# # sst_mean_rob_tr <- crop(sst_mean_rob_tr, extent(pt))
#
# omap_plt + as.layer(spplot(sst_mean_rob,
#                            col.regions = envinmrPalette(1000),
#                            at = seq(-5, 35, 0.1), alpha = 0.6))
#
#
#
#
# ### get data TZA
# tza <- raster::getData('GADM', country = 'TZA', level = 0)
# tza_ul <- c(bbox(tza)[2, 2], bbox(tza)[1, 1])
# tza_lr <- c(bbox(tza)[2, 1], bbox(tza)[1, 2])
# tza_rst <- openmap(upperLeft = tza_ul, lowerRight = tza_lr,
#                    minNumTiles = 15, type = "bing")
# tza_rst <- raster(tza_rst)
# tza_rst_utm <- projectRaster(tza_rst, crs = "+init=epsg:32737")
# tza_utm <- spTransform(tza, projection(tza_rst_utm))
#
# ### mask raster to TZA borders & convert to sp.layout
# tza_rst_utm_msk <- mask(tza_rst_utm, tza_utm)
# lout_tza_utm <- rgb2spLayout(tza_rst_utm_msk)
# lout_tza_utm[[2]] <- gsub("#010101", "transparent", lout_tza_utm[[2]])
#
# ### plot TZA with bing backgraound
# tza_p <- spplot(tza_utm, zcol = "PID", fill = "transparent",
#                 colorkey = FALSE, col = "black", lwd = 1.5,
#                 sp.layout = lout_tza_utm,
#                 par.settings = list(axis.line = list(col = 0)))
#
#
#
#
#
#
# ### fulldisk plot
# africa_rst_geos <- projectRaster(africa_rst, crs = p)
# tza_geos <- spTransform(tza, p)
# tza_geos@bbox[1, 1] <- extent(africa_rst_geos)@xmin
# tza_geos@bbox[1, 2] <- extent(africa_rst_geos)@xmax
# tza_geos@bbox[2, 1] <- extent(africa_rst_geos)@ymin
# tza_geos@bbox[2, 2] <- extent(africa_rst_geos)@ymax
#
# lout_africa_geos <- rgb2spLayout(africa_rst_geos)
# lout_africa_geos[[2]] <- gsub("#010101", "transparent",
#                               lout_africa_geos[[2]])
#
# fdsk_p <- spplot(tza_geos, zcol = "PID", fill = "white", alpha.regions = 0.5,
#                  colorkey = FALSE, col = "black", alpha = 1, lwd = 1.2,
#                  sp.layout = lout_africa_geos,
#                  par.settings = list(axis.line = list(col = 0))) +
#   layer(panel.text(x = coordinates(tza_geos)[, 1],
#                    y = coordinates(tza_geos)[, 2],
#                    labels = "TZA", cex = 1.2))
#
#
#
#
#
# ### kili aerial
# kili <- kiliAerial(projection = "+init=epsg:32737",
#                    minNumTiles = 20, rasterize = TRUE,
#                    template = extent(c(36.99033, 37.74099,
#                                        -3.465092, -2.83096)))
#
# lout_kili <- rgb2spLayout(kili)
#
#
# ### KiLi plots
# # lyr <- ogrListLayers("/home/ede/tappelhans/uni/temp/data/plots_shp/PlotPoles_ARC1960_mod_20140807_final.shp")
# # plots <- readOGR("/home/ede/tappelhans/uni/temp/data/plots_shp/PlotPoles_ARC1960_mod_20140807_final.shp",
# #                  layer = lyr)
# # proj4string(plots) <- "+init=epsg:21037"
# # plots <- spTransform(plots, CRS(projection(kili)))
# # plots <- subset(plots, PoleType == "AMP")
# # plots <- plots[plots$PlotID %in% c("cof3", "foc0", "fod2", "gra1",
# #                                    "hom4", "mai0", "sav0"), ]
# # plots@data$X <- coordinates(plots)[, 1]
# # plots@data$Y <- coordinates(plots)[, 2]
#
#
# ### locations of met sites
# # met_sites <- data.frame(name = c("KIA", "MOSHI",
# #                                  as.character(substr(plots@data$PlotID, 1, 3))),
# #                         x = c(284634, 314296, coordinates(plots)[, 1]),
# #                         y = c(9620955, 9628204, coordinates(plots)[, 2]))
# # met_sites$name <- factor(toupper(met_sites$name),
# #                          levels = toupper(as.character(met_sites$name)))
# # coordinates(met_sites) <- ~x+y
# # met_sites@proj4string <- CRS(projection(kili))
# # met_sites@bbox[1, 1] <- extent(kili)@xmin
# # met_sites@bbox[1, 2] <- extent(kili)@xmax
# # met_sites@bbox[2, 1] <- extent(kili)@ymin
# # met_sites@bbox[2, 2] <- extent(kili)@ymax
#
#
# ### merge point data
#
#
# c(xFromCol(kili, col = round(ncol(kili) / 2, 0)),
#   xFromCol(kili, col = round(ncol(kili) / 2, 0))),
# c(yFromRow(kili, row = round(nrow(kili) / 2, 0)),
#   yFromRow(kili, row = round(nrow(kili) / 2, 0)))
#
# met_sites <- data.frame(x = xFromCol(kili, col = round(ncol(kili) / 2, 0)),
#                         y = yFromRow(kili, row = round(nrow(kili) / 2, 0)),
#                         name = "kili")
# coordinates(met_sites) <- ~ x + y
#
# ### kili plot
# # cols <- c("white", "white", rep("grey70", 7)) #colorRampPalette(brewer.pal(7, "Paired"))(7))
#
# cols <- "transparent"
#
# num_xmin <- xmin(kili) - 1000
# num_xmax <- xmax(kili) + 1000
# num_xlim <- c(num_xmin, num_xmax)
#
# num_ymin <- ymin(kili) - 1000
# num_ymax <- ymax(kili) + 1000
# num_ylim <- c(num_ymin, num_ymax)
#
# kili_p <- spplot(met_sites, zcol = "name", sp.layout = lout_kili,
#                  col.regions = cols,
#                  pch = c(22, 22, rep(21, 7)),
#                  cex = c(1, 1, rep(0.7, 7)),
#                  edge.col = "transparent",
#                  xlim = num_xlim, ylim = num_ylim,
#                  par.settings = list(axis.line = list(col = 0)),
#                  auto.key = FALSE)
# #                  key.space = list(x = 0.075, y = 0.8, corner = c(0, 0),
# #                                   background = "white", border = "black",
# #                                   alpha.background = 0.5,
# #                                   cex = 0.75, columns = 3,
# #                                   between = 0.5, between.columns = 0.5))
#
# kili_p
#
#
#
# ### update tza_p with panel.extent kili
# tza_p2 <- tza_p +
#   #   layer(panel.extent(extent(kili),
#   #                      col = "white", alpha = 0.5)) +
#   layer(panel.extent(extent(kili)))
#
# # map_data <- list(globe_geos = africa_rst_geos,
# #                  tza_geos = tza_geos,
# #                  tza_utm = tza_utm,
# #                  lout_tza_utm = lout_tza_utm,
# #                  lout_africa_geos = lout_africa_geos,
# #                  kili = kili,
# #                  lout_kili = lout_kili,
# #                  met_sites = met_sites)
# #
# # save(map_data, file = "software/r_stuff/bitsnpieces/map_data.rda")
#
#
# ### if map_data.rda exists, can start from here:
# # load("software/r_stuff/bitsnpieces/map_data.rda")
# #
# # tza_utm <- map_data$tza_utm
# # globe_geos <- map_data$globe_geos
# # tza_geos <- map_data$tza_geos
# # lout_tza_utm <- map_data$lout_tza_utm
# # lout_africa_geos <- map_data$lout_africa_geos
# # kili <- map_data$kili
# # lout_kili <- map_data$lout_kili
#
#
# png("software/r_stuff/bitsnpieces/kili_overview_bing_globe.png",
#     width = 7, height = 7, units = "in", res = 600)
# grid.newpage()
# print(fdsk_p, newpage = FALSE)
#
# downViewport(trellis.vpname("figure"))
# # grid.rect(x = x_tz, y = y_tz, width = wdth_tz, height = hgt_tz,
# #           gp = gpar(col = "black", fill = "transparent"),
# #           just = c("left", "top"))
#
# vp1 <- viewport(x = 1, y = 0.275, width = 0.5, height = 0.5,
#                 just = c("right", "centre"), name = "insetTZA")
#
# pushViewport(vp1)
# grid.rect(gp = gpar(col = "black", fill = "white", alpha = 0.75))
# print(tza_p2, newpage = FALSE)
#
# # downViewport(trellis.vpname("figure"))
# # grid.rect(gp = gpar(fill = "transparent"))
#
# upViewport(0)
#
# moveToGrob(x = 1, y = 1, name = panel.extent(extent(kili))$name, vp = "insetTZA")
# lineToGrob(x = 1, y = 0, name = "tst", vp = "insetTZA")
# grid.line.to(1, 0, name = "tst")
#
# upViewport(1)
#
# vp2 <- viewport(x = 1, y = 0.975, width = 0.5, height = 0.425,
#                 just = c("right", "top"), name = "insetKili")
#
# pushViewport(vp2)
# grid.rect(gp = gpar(col = "black", fill = "white", alpha = 0.75))
#
# # upViewport(1)
# #
# # vp3 <- viewport(x = 0.75, y = 0.7625, width = 0.575, height = 0.5,
# #                 just = c("centre", "centre"), name = "inset3")
# #
# # pushViewport(vp3)
#
# mar_theme <- lattice.options()
# mar_theme$layout.widths$left.padding$x <- 2
# mar_theme$layout.widths$left.padding$units <- "points"
# mar_theme$layout.widths$right.padding$x <- 2
# mar_theme$layout.widths$right.padding$units <- "points"
# mar_theme$layout.heights$top.padding$x <- 2
# mar_theme$layout.heights$top.padding$units <- "points"
# mar_theme$layout.heights$bottom.padding$x <- 2
# mar_theme$layout.heights$bottom.padding$units <- "points"
#
#
# print(update(kili_p, lattice.options = mar_theme), newpage = FALSE)
# # downViewport(trellis.vpname(name = "figure"))
# # offsetGridText(x = coordinates(met_sites), xlim = num_xlim, ylim = num_ylim,
# #                labels = met_sites$name, stext = TRUE, offset = .0175,
# #                gp = gpar(fontsize = 6), r = 0.05)
#
# upViewport(1)
# #pushViewport(vp1)
#
# urx <- 1 - (extent(kili)@xmax / extent(tza_utm)@xmax)
# ury <- extent(kili)@ymax / extent(tza_utm)@ymax
#
# grid.move.to(x = urx, y = ury, default.units = "native", name = NULL,
#              draw = TRUE)
# grid.line.to(x = 1, y = 0, draw = TRUE)#, vp = "insetKili")
#
#
#
# dev.off()
#
