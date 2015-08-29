library(raster)
library(remote)
library(reshape2)
library(latticeExtra)
library(Rsenal)
library(dplyr)
library(grid)
library(lubridate)

path_local <- "/media/ede/tims_ex/eot_kili_analysis"
setwd(path_local)

load("five_yr_chunks/results/exp_var_overview_dns80_kili.rda")

out$mat <- out$mat[24:1, ]

### plot exp var matrix
mat_anom <- out$mat * 100 #- mean(out$mat)
df <- melt(mat_anom)
names(df) <- c("chunk", "lag", "value")
df$chunk <- factor(df$chunk, levels = rev(levels(df$chunk)))
mat <- acast(df, chunk ~ lag, value.var = "value")
mx <- 60

clrs <- colorRampPalette(c("black", rev(brewer.pal(9, "YlGnBu")), "white"))

mat_plt <- levelplot(t(mat), col.regions = clrs(1000),
                     ylab = "Chunk", xlab = "Lag",
                     #main = "             Coefficient of determination",
                     colorkey = FALSE, # list(space = "top", draw = TRUE,
                                     # height = 0.75, width = 1),
                     scales = list(x = list(rot = 45)),
                     at = seq(0, mx, 10))

x_at <- c(3, 1, 1, 3, 4, 6, 7, 12, 10, 11, 11)
y_at <- c(13, 12, 11, 11, 11, 11, 11, 13, 12, 12, 11)
labs <- c("a)", "b)", "b)", "c)", "c)", "d)", "d)", "e)", "e)", "e)", "e)")

mat_plt <- mat_plt + layer(panel.text(x = x_at, y = y_at, label = labs))

mat_upr <- mat
q95 <- quantile(mat_upr, 0.95)
mat_upr[mat_upr < q95] <- NA
mat_plt_upr <- levelplot(t(mat_upr), #col.regions = clrs(1000),
                         ylab = "Chunk", xlab = "Lag",
                         region = TRUE, border = "black",
                         border.lwd = 1,
                         col.regions = "transparent",
                         colorkey = list(space = "top", draw = TRUE,
                                         height = 0.75, width = 1),
                         scales = list(x = list(rot = 45)),
                         at = seq(0, mx, length.out = 1000))

mat_lwr <- mat
q05 <- quantile(mat_lwr, 0.05)
mat_lwr[mat_lwr > q05] <- NA
mat_plt_lwr <- levelplot(t(mat_lwr), #col.regions = clrs(1000),
                         ylab = "Chunk", xlab = "Lag",
                         region = TRUE, border = "white",
                         border.lwd = 1, border.lty = 3,
                         col.regions = "transparent",
                         colorkey = list(space = "top", draw = TRUE,
                                         height = 0.75, width = 1),
                         scales = list(x = list(rot = 45)),
                         at = seq(0, mx, length.out = 1000))

exv_plt <- mat_plt + as.layer(mat_plt_upr) + as.layer(mat_plt_lwr)

#### correlation matrix
mat_r <- out$mat_r
clrs_r <- colorRampPalette(brewer.pal(9, "RdBu"))

mat_r_plt <- levelplot(t(mat_r), col.regions = clrs_r(1000),
                       ylab = "Chunk", xlab = "Lag",
                       #main = "             Coefficient of determination",
                       colorkey = FALSE, #list(space = "top", draw = TRUE,
                                       #height = 0.75, width = 1),
                       scales = list(x = list(rot = 45)),
                       at = seq(-1, 1, 0.25))

mat_r_plt <- mat_r_plt + layer(panel.text(x = x_at, y = y_at, label = labs))

r_plt <- mat_r_plt + as.layer(mat_plt_upr) + as.layer(mat_plt_lwr)

plt_fin <- print(latticeCombineGrid(list(exv_plt, r_plt), layout = c(2, 1)))
plt_fin <- update(plt_fin, scales = list(alternating = 1))

# png("five_yr_chunks/graphs/exp_var_r_overview_kili.png",
#     width = 20, height = 20, units = "cm", res = 600)
# grid.newpage()
# vp_ypos <- 1.1
#
# print(plt_fin)
# downViewport(trellis.vpname(name = "panel", column = 1, row = 1))
# #grid.rect(gp = gpar(fill = "black"))
# vp_key_left <- viewport(x = 0.5, y = vp_ypos, width = 0.75, height = 0.1,
#                         clip = "off", name = "colorkey_left")
# pushViewport(vp_key_left)
# #grid.rect(gp = gpar(fill = "black"))
#
# draw.colorkey(key = list(col = clrs(1000), width = 1,
#                          at = seq(0, mx, 10),
#                          space = "top"), draw = TRUE)
# grid.text("Explained variance [%]",
#           x = 0.5, y = 1, just = c("centre", "bottom"))
#
# upViewport(0)
#
# downViewport(trellis.vpname(name = "panel", column = 2, row = 1))
# #grid.rect(gp = gpar(fill = "black"))
# vp_key_right <- viewport(x = 0.5, y = vp_ypos, width = 0.75, height = 0.1,
#                         clip = "off", name = "colorkey_right")
# pushViewport(vp_key_right)
# #grid.rect(gp = gpar(fill = "black"))
#
# draw.colorkey(key = list(col = clrs_r(1000), width = 1,
#                          at = seq(-1, 1, 0.25),
#                          space = "top"), draw = TRUE)
# grid.text("Correlation coefficient",
#           x = 0.5, y = 1, just = c("centre", "bottom"))
# dev.off()

tiff("five_yr_chunks/graphs/exp_var_r_overview_kili.tif",
    width = 20, height = 20, units = "cm", res = 600,
    compression = "lzw")
grid.newpage()
vp_ypos <- 1.1

print(plt_fin)
downViewport(trellis.vpname(name = "panel", column = 1, row = 1))
#grid.rect(gp = gpar(fill = "black"))
vp_key_left <- viewport(x = 0.5, y = vp_ypos, width = 0.75, height = 0.1,
                        clip = "off", name = "colorkey_left")
pushViewport(vp_key_left)
#grid.rect(gp = gpar(fill = "black"))

draw.colorkey(key = list(col = clrs(1000), width = 1,
                         at = seq(0, mx, 10),
                         space = "top"), draw = TRUE)
grid.text("Explained variance [%]",
          x = 0.5, y = 1, just = c("centre", "bottom"))

upViewport(0)

downViewport(trellis.vpname(name = "panel", column = 2, row = 1))
#grid.rect(gp = gpar(fill = "black"))
vp_key_right <- viewport(x = 0.5, y = vp_ypos, width = 0.75, height = 0.1,
                         clip = "off", name = "colorkey_right")
pushViewport(vp_key_right)
#grid.rect(gp = gpar(fill = "black"))

draw.colorkey(key = list(col = clrs_r(1000), width = 1,
                         at = seq(-1, 1, 0.25),
                         space = "top"), draw = TRUE)
grid.text("Correlation coefficient",
          x = 0.5, y = 1, just = c("centre", "bottom"))
dev.off()
