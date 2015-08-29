library(raster)
library(remote)
library(reshape2)
library(latticeExtra)
library(Rsenal)
library(dplyr)
library(grid)

path_local <- "/media/ede/tims_ex/eot_kili_analysis"
setwd(path_local)

load("five_yr_chunks/results/exp_var_overview_dns80_kili.rda")

out$mat <- out$mat[24:1, ]

# rst <- raster(out$mat)
# rst33 <- focal(rst, w = standardFocalWindow(3, 1), fun = max, pad = TRUE, na.rm = TRUE)
#
# out$mat_focal <- as.matrix(rst33)
# rownames(out$mat_focal) <- rownames(out$mat)
# colnames(out$mat_focal) <- colnames(out$mat)

df <- melt(out$mat)
names(df) <- c("chunk", "lag", "value")
df <- df[order(df$chunk, df$lag), ]

df_ord <- df[order(df$value, decreasing = TRUE), ]
df_ord$chlag <- paste0(df_ord$chunk, df_ord$lag)

qupr <- quantile(df_ord$value, c(0.95, 0.9, 0.85))
qlwr <- quantile(df_ord$value, 0.05)

# dfupr <- df_ord[df_ord$value > qupr, ]


### mode 1
plts_upr <- list()
plts_upr_all <- list()


for (j in seq(qupr)) {
  print(j)
  if (j == 1) df_rest <- df_ord

  dfupr <- df_rest[df_rest$value > qupr[j], ]
  df_rest <- df_rest[!df_rest$chlag %in% dfupr$chlag, ]

  dfupr <- dfupr[order(dfupr$chunk, dfupr$lag), ]

  ptrns <- paste(dfupr$chunk, dfupr$lag, sep = "_")

  ### overview over all upper patterns
  mds_upr <- lapply(seq(nrow(dfupr)), function(i) {
    ptrn <- paste("*", paste("kili", dfupr[i, ]$chunk,
                             dfupr[i, ]$lag, sep = "_"),
                  ".rda", sep = "")
    print(ptrn)
    fl <- list.files("five_yr_chunks/results/modes",
                     pattern = ptrn, full.names = TRUE)
    load(fl)
    md <- modes@modes$mode_01@rsq_predictor

    plt <- spplot(md, col.regions = envinmrPalette(1000),
                  at = seq(0, 1, 0.01))


    return(list(mode = md,
                plt = plt))
  })

  extrct <- function(x) {
    function(y) lapply(seq(y), function(h) y[[h]][[x]])
  }

  extrct_plt <- extrct("plt")
  plts_upr_all[[j]] <- latticeCombineGrid(extrct_plt(mds_upr)) +
    layer(grid.text(ptrns[panel.number()], .5, .1))

  extrct_md <- extrct("mode")
  mds <- extrct_md(mds_upr)
  mds_upr_max <- calc(stack(mds), fun = median, na.rm = TRUE)
  plts_upr[[j]] <- spplot(mds_upr_max, col.regions = envinmrPalette(1000),
                          at = seq(0, 0.5, 0.01),
                          main = expression("(median) R"^2))

}

plts_upr_all[[1]]

latticeCombineGrid(plts_upr, layout = c(1, 3))

### mode 2
plts_upr2 <- list()

for (j in seq(qupr)) {
  print(j)
  if (j == 1) df_rest <- df_ord

  dfupr <- df_rest[df_rest$value > qupr[j], ]
  df_rest <- df_rest[!df_rest$chlag %in% dfupr$chlag, ]

  ### overview over all upper patterns
  mds_upr <- lapply(seq(nrow(dfupr)), function(i) {
    ptrn <- paste("*", paste("kili", dfupr[i, ]$chunk,
                             dfupr[i, ]$lag, sep = "_"),
                  ".rda", sep = "")
    print(ptrn)
    fl <- list.files("five_yr_chunks/results/modes",
                     pattern = ptrn, full.names = TRUE)
    load(fl)
    modes@modes$mode_02@rsq_predictor
  })

  mds_upr_max <- calc(stack(mds_upr), fun = median, na.rm = TRUE)
  plts_upr2[[j]] <- spplot(mds_upr_max, col.regions = envinmrPalette(1000),
                          at = seq(0, 0.25, 0.01))

}

latticeCombineGrid(list(plts_upr[[1]], plts_upr2[[1]]), layout = c(1, 2))

### overview over all lower patterns
dflwr <- df_ord[df_ord$value < qlwr, ]
mds_lwr <- lapply(seq(nrow(dflwr)), function(i) {
  ptrn <- paste("*", paste("kili", dflwr[i, ]$chunk,
                           dflwr[i, ]$lag, sep = "_"),
                ".rda", sep = "")
  print(ptrn)
  fl <- list.files("five_yr_chunks/results/modes",
                   pattern = ptrn, full.names = TRUE)
  load(fl)
  modes@modes$mode_01@rsq_predictor
})

mds_lwr_max <- calc(stack(mds_lwr), fun = max, na.rm = TRUE)
spplot(mds_lwr_max, col.regions = envinmrPalette(1000), at = seq(0, 1, 0.01))

