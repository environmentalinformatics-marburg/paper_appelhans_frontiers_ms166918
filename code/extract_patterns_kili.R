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

### plot exp var matrix
mat_anom <- out$mat - mean(out$mat)
df <- melt(mat_anom)
names(df) <- c("chunk", "lag", "value")
df$chunk <- factor(df$chunk, levels = rev(levels(df$chunk)))
mat <- acast(df, chunk ~ lag, value.var = "value")
mx <- 0.35

clrs <- colorRampPalette(brewer.pal(9, "PRGn"))

mat_plt <- levelplot(t(mat), col.regions = clrs(1000),
                 ylab = "Chunk", xlab = "Lag",
                 #main = paste("Anomaly in explained variance (mode 1 + 2)"),
                 colorkey = list(space = "top", draw = TRUE,
                                 height = 0.75, width = 1),
                 scales = list(x = list(rot = 45)),
                 at = seq(-mx, mx, 0.1))

mat_upr <- mat
q95 <- quantile(mat_upr, 0.95)
mat_upr[mat_upr < q95] <- NA
mat_plt_upr <- levelplot(t(mat_upr), #col.regions = clrs(1000),
                         ylab = "Chunk", xlab = "Lag",
                         region = TRUE, border = "black",
                         border.lwd = 1,
                         col.regions = "transparent",
                         main = paste("Coefficient of determination (mode 1 + 2)"),
                         colorkey = list(space = "top", draw = TRUE,
                                         height = 0.75, width = 1),
                         scales = list(x = list(rot = 45)),
                         at = seq(0, 0.35, length.out = 1000))

mat_lwr <- mat
q05 <- quantile(mat_lwr, 0.05)
mat_lwr[mat_lwr > q05] <- NA
mat_plt_lwr <- levelplot(t(mat_lwr), #col.regions = clrs(1000),
                         ylab = "Chunk", xlab = "Lag",
                         region = TRUE, border = "black",
                         border.lwd = 1, border.lty = 3,
                         col.regions = "transparent",
                         main = paste("Coefficient of determination (mode 1 + 2)"),
                         colorkey = list(space = "top", draw = TRUE,
                                         height = 0.75, width = 1),
                         scales = list(x = list(rot = 45)),
                         at = seq(0, 0.35, length.out = 1000))

tiff("five_yr_chunks/graphs/exp_var_overview_kili.tif",
     width = 15, height = 20, units = "cm", res = 300, compression = "lzw")
mat_plt + as.layer(mat_plt_upr) + as.layer(mat_plt_lwr)
dev.off()


# rst <- raster(out$mat)
# rst33 <- focal(rst, w = standardFocalWindow(3, 1), fun = max, pad = TRUE, na.rm = TRUE)
#
# out$mat_focal <- as.matrix(rst33)
# rownames(out$mat_focal) <- rownames(out$mat)
# colnames(out$mat_focal) <- colnames(out$mat)

df <- melt(out$mat)
names(df) <- c("chunk", "lag", "value")
df <- df[order(df$chunk, df$lag), ]

### overview over all patterns
# mds <- lapply(seq(nrow(df)), function(i) {
#   ptrn <- paste("*", paste("kili", df[i, ]$chunk,
#                            df[i, ]$lag, sep = "_"),
#                 ".rda", sep = "")
#   print(ptrn)
#   fl <- list.files("five_yr_chunks/results/modes",
#                    pattern = ptrn, full.names = TRUE)
#   load(fl)
#   modes@modes$mode_01@rsq_predictor
# })
#
#
#
# plts <- lapply(seq(mds), function(i) {
#   print(i)
#   spplot(mds[[i]], col.regions = envinmrPalette(1000),
#          at = seq(0, 1, 0.01))
# })
#
# png("five_yr_chunks/graphs/all_modes_eot1.png",
#     width = 1189 * 1.5, height = 841 * 1.5, units = "mm", res = 300)
# grid.newpage()
# latticeCombineGrid(plts, c(13, 24))
#
# for (i in 1:13) {
#   print(i)
#   trellis.focus("panel", column = i, row = 1, clip.off = TRUE)
#   grid.text(colnames(out$mat)[i], x = unit(0.5, "npc"),
#             y = unit(1.1, "npc"))
# }
#
# for (i in 1:24) {
#   print(i)
#   trellis.focus("panel", column = 1, row = i, clip.off = TRUE)
#   grid.text(rownames(out$mat)[i], x = unit(-0.1, "npc"),
#             y = unit(0.5, "npc"), rot = 90)
# }
#
# trellis.unfocus()
#
# dev.off()

### > 0.40
df50 <- data.frame(chunk = df$chunk[df$value >= 0.50])
df50$lag <- df$lag[df$value >= 0.50]
df50 <- df50[order(df50$chunk, decreasing = TRUE), ]
df50$id <- 1:nrow(df50)

# ind_chunk_40 <- c("2002-2006", "2003-2007", "2004-2008", "2005-2009")
# ind_lag_40 <- c("lag04", "lag05", "lag06", "lag07")
# df50 <- df50[df50$chunk %in% ind_chunk_40, ]
# df50 <- df50[!df50$lag %in% ind_lag_40, ]


# ### > 0.25
# df25 <- data.frame(chunk = df$chunk[df$value >= 0.25])
# df25$lag <- df$lag[df$value >= 0.25]
#
# ind_chunk_25 <- c("1993-1997", "1994-1998", "1995-1999", "1996-2000")
# ind_lag_25 <- c("lag03", "lag04", "lag05", "lag06", "lag07", "lag08")
# df25 <- df25[df25$chunk %in% ind_chunk_25, ]
# df25 <- df25[!df25$lag %in% ind_lag_25, ]

# patterns <- list(df50[c(1, 2), ],
#                  df50[c(3:5, 9, 11), ],
#                  df50[c(8, 10), ],
#                  df50[c(13, 14, 17:19), ],
#                  df50[c(15, 16, 20), ])

# patterns <- list(df50[c(19, 20, 17, 9, 10, 2), ],
#                  df50[c(18, 11, 12, 4), ],
#                  df50[c(13, 14, 6, 7, 1), ],
#                  df50[c(22, 23, 24, 25, 26), ])#,
#                  df50[c(8, 15), ])

patterns <- list(df50[c(19, 20, 17), ],
                 df50[c(9, 10), ],
                 df50[c(18, 11, 12, 4), ])#,
                 df50[c(7, 1), ],
                 df50[c(15, 8), ])
###########################################################################



##### MODE 1 ##############################################################
### r response to check direction of relationship
rsp <- lapply(seq(patterns), function(i) {
  do.call("stack", lapply(seq(nrow(patterns[[i]])), function(j) {
    ptrn <- paste("*", paste("kili", patterns[[i]]$chunk[j],
                             patterns[[i]]$lag[j], sep = "_"),
                  ".rda", sep = "")
    print(ptrn)
    fl <- list.files("five_yr_chunks/results/modes",
                     pattern = ptrn, full.names = TRUE)
    load(fl)
    print(modes@modes$mode_01@r_response)
    modes@modes$mode_01@r_response
  }))
})

clrs_rsp <- colorRampPalette(brewer.pal(5, "PuOr"))

plts_rsp <- lapply(seq(rsp), function(i) {
  lst <- lapply(seq(nlayers(rsp[[i]])), function(j) {
    spplot(rsp[[i]][[j]], col.regions = clrs_rsp(1000),
           at = seq(-0.8, 0.8, 0.01))
  })

  return(latticeCombineGrid(lst))
})

plts_rsp[[1]]


plts_rsp_mn <- lapply(seq(rsp), function(i) {
  spplot(mean(rsp[[i]]), col.regions = clrs_rsp(1000),
         at = seq(-1, 1, 0.01),
         colorkey = FALSE) #+ #list(space = "left", height = 1/6*2, x = 0)) +
    #layer(panel.map("world", border = "white", col = "transparent", lwd = 0.5)) +
    #layer(panel.extent(extent(ext), border = "black", lwd = 2, lty = 1)) +
#     as.layer(spplot(mean(mds[[i]]), col.regions = "transparent",
#                     at = seq(tsh[i], 1, 0.2), contour = TRUE))
})

plts_rsp_mn[[1]]

################



### load modes
mds <- lapply(seq(patterns), function(i) {
  do.call("stack", lapply(seq(nrow(patterns[[i]])), function(j) {
    ptrn <- paste("*", paste("kili", patterns[[i]]$chunk[j],
                             patterns[[i]]$lag[j], sep = "_"),
                  ".rda", sep = "")
    print(ptrn)
    fl <- list.files("five_yr_chunks/results/modes",
                     pattern = ptrn, full.names = TRUE)
    load(fl)
    modes@modes$mode_01@rsq_predictor
  }))
})



plts <- lapply(seq(mds), function(i) {
  lst <- lapply(seq(nlayers(mds[[i]])), function(j) {
    spplot(mds[[i]][[j]], col.regions = envinmrPalette(1000),
           at = seq(0, 1, 0.01))
  })

  return(latticeCombineGrid(lst))
})

plts[[1]]




tsh <- c(0.7, 0.7, 0.7) # 0.6)

load("five_yr_chunks/data/stck_sst_dsn_dns80.rda")
load("five_yr_chunks/data/stck_chirps_dsn_dns_kili80.rda")

ext <- extent(stck_prcp_kili_dsn_dns)
clrs <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))
clrs <- envinmrPalette

plts_mn <- lapply(seq(mds), function(i) {
  spplot(mean(mds[[i]]), col.regions = clrs(1000),
         at = seq(0, 1, 0.01),
         colorkey = FALSE) + #list(space = "left", height = 1/6*2, x = 0)) +
    layer(panel.map("world", border = "white", col = "grey60", lwd = 0.5)) +
    layer(panel.extent(extent(ext), border = "black", lwd = 2, lty = 1)) +
    as.layer(spplot(mean(mds[[i]]), col.regions = "transparent",
                     at = seq(tsh[i], 1, 0.2), contour = TRUE))
})

plts_mn[[1]]

latticeCombineGrid(plts_mn)

ts_vals <- do.call("rbind", lapply(seq(patterns), function(i) {
  st <- substr(patterns[[i]]$chunk, 1, 4)
  nd <- as.character(as.numeric(st) + 4)
  return(data.frame(st, nd, stringsAsFactors = FALSE))
}))

sq <- seq(min(as.numeric(ts_vals$st)), max(as.numeric(ts_vals$nd)), 1)


###
chnks_sst <- lapply(seq(patterns), function(i) {
#   st <- substr(patterns[[i]]$chunk, 1, 4)
#   nd <- as.character(as.numeric(st) + 4)
#   sq <- seq(min(as.numeric(st)), max(as.numeric(nd)), 1)
#   print(sq)
  ptrn <- paste("*", sq, "*", ".rst", sep = "")
#   fls <- do.call("c", lapply(ptrn, function(j) {
#     list.files("source_data/sst_1982_2010",
#                pattern = glob2rx(j),
#                full.names = TRUE)
#   }))
  out <- stck_sst_dsn_dns[[which(substr(names(stck_sst_dsn_dns),
                                        5, 8) %in% sq)]]
  print(names(out))
  return(out)
})

chnks_gpcp <- lapply(seq(patterns), function(i) {
#   st <- substr(patterns[[i]]$chunk, 1, 4)
#   nd <- as.character(as.numeric(st) + 4)
#   sq <- seq(min(as.numeric(st)), max(as.numeric(nd)), 1)
#   print(sq)
  ptrn <- paste("*", sq, "*", ".rst", sep = "")
  #   fls <- do.call("c", lapply(ptrn, function(j) {
  #     list.files("source_data/sst_1982_2010",
  #                pattern = glob2rx(j),
  #                full.names = TRUE)
  #   }))
  out <- stck_prcp_kili_dsn_dns[[which(substr(names(stck_prcp_kili_dsn_dns),
                                        13, 16) %in% sq)]]
  print(names(out))
  return(out)
})



mds_cln_mask <- lapply(seq(mds), function(i) {
  out <- mean(mds[[i]])
  out[out < tsh[i]] <- NA
  out[out >= tsh[i]] <- 1
  return(out)
})


ts_plt_lst <- lapply(seq(patterns), function(i) {
  cells <- which(mds_cln_mask[[i]][] == 1)

  if (length(cells) == 0) return(NULL) else {

    sst_ts <- apply(X = extract(chnks_sst[[i]], cells), MARGIN = 2,
                    FUN = mean, na.rm = TRUE)
    gpcp_ts <- apply(chnks_gpcp[[i]][], 2, mean, na.rm = TRUE)
    x_by <- max(as.numeric(substr(patterns[[i]]$lag, 4, 5)))
    print(x_by)
    #   if (i == 6) {
    #     st <- as.numeric(which(sst_ts == min(sst_ts)))
    #   } else {
    #     st <- as.numeric(which(abs(sst_ts) == max(abs(sst_ts))))
    #   }
    sst_ts_srt <- sort(abs(sst_ts), decreasing = TRUE)
    st1 <- as.numeric(which(abs(sst_ts) == sst_ts_srt[1]))
    st2 <- as.numeric(which(abs(sst_ts) == sst_ts_srt[2]))
    nd1 <- st1 + x_by
    nd2 <- st2 + x_by
    labs_st <- substr(names(sst_ts)[1], 5, 8)
    labs_nd <- substr(names(sst_ts)[length(sst_ts)], 5, 8)
    xyplot(sst_ts ~ seq(sst_ts), type = "l", lwd = 2,
           col = "black", ylim = c(-3.5, 3.5),
           panel = function(...) {
             if (x_by > 0) {
               panel.polygon(x = c(st1, st1, nd1, nd1),
                             y = c(-3, 3, 3, -3),
                             col = "grey90", border = "grey30")
#                panel.polygon(x = c(st2, st2, nd2, nd2),
#                              y = c(-3, 3, 3, -3),
#                              col = "grey90", border = "grey30")
             } else {
               panel.abline(v = c(st1, st2), lwd = 1.5, col = "grey90")
             }
             #            panel.abline(v = seq(st_ovr, nd_ovr, 1), lty = 1,
             #                         lwd = 0.5, col = "grey")
             panel.xyplot(...)
             panel.text(x = rep(-1, 3),
                        y = c(-1, 0, 1),
                        labels = c("-1", "0", "1"))
             panel.text(x = 0, y = 1.8, adj = c(0, 0.5),
                        labels = paste(labs_st, labs_nd, sep = " - "))
           }) +
      as.layer(xyplot(gpcp_ts / 100 ~ seq(gpcp_ts),
                      type = "l", lty = 3, lwd = 2)) +
      layer_(panel.abline(h = 0, lty = 2, col = "grey")) +
      layer_(panel.abline(h = c(1, 2, 3), lty = 1, col = "grey")) +
      layer_(panel.abline(h = c(-1, -2, -3), lty = 1, col = "grey"))
  }
})


final <- latticeCombineGrid(list(plts_mn[[1]], plts_rsp_mn[[1]], ts_plt_lst[[1]],
                                 plts_mn[[2]], plts_rsp_mn[[2]], ts_plt_lst[[2]],
                                 plts_mn[[3]], plts_rsp_mn[[3]], ts_plt_lst[[3]]),#,
                                 #plts_mn[[4]], ts_plt_lst[[4]]),
                                 #plts_mn[[5]], ts_plt_lst[[5]]),
                            layout = c(3, 3))

tiff("five_yr_chunks/graphs/patterns_mode1_w_ts.tif",
     width = 174, height = 234, units = "mm", res = 600,
     compression = "lzw")
grid.newpage()
print(final)

downViewport(trellis.vpname("figure"))
vp1 <- viewport(x = 0.25, y = 1.02, height = 0.05, width = 0.3,
                just = c("centre", "bottom"),
                name = trellis.grobname("left", type = "colorkey"))

pushViewport(vp1)
# grid.rect()
draw.colorkey(key =  list(col = envinmrPalette(1000), width = 1,
                          at = seq(0, 1, 0.01),
                          space = "top"), draw = TRUE)

t <- current.vpTree(all = FALSE)
downViewport(t$children[[2]]$children[[1]]$parent$name)
grid.lines(x = 0.65, y = c(0,1), draw = TRUE)
upViewport(n = 2)

pushViewport(viewport(x = 0.5, y = 2, height = 1, width = 1,
                      just = c("centre", "top"), name = "main"))
grid.text(label = "Mean coefficient of determination")
upViewport(2)
vp2 <- viewport(x = 0.75, y = 1.02, height = 0.05, width = 0.3,
                just = c("centre", "bottom"), name = "right")

pushViewport(vp2)
# grid.rect()
draw.key(key = list(col = c("black", "cornflowerblue"),
                    lines = list(lty = c(1, 2), lwd = 2),
                    text = list(c("SST anomaly [K]",
                                  "GPCP anomaly [mm/d]"),
                                col = "black")), draw = TRUE)
upViewport(0)
dev.off()




### time series correlation
st <- 1993
nd <- 2001

ts_lst <- lapply(seq(patterns), function(i) {
  cells <- which(mds_cln_mask[[i]][] == 1)
  x_by <- max(as.numeric(substr(patterns[[i]]$lag, 4, 5)))
  sst_st <- stck_sst_dsn_dns[[which(substr(names(stck_prcp_kili_dsn_dns),
                                           13, 16) %in% sq) - x_by]]

  if (length(cells) == 0) return(NULL) else {

    sst_ts <- apply(X = extract(sst_st, cells), MARGIN = 2,
                    FUN = mean, na.rm = TRUE)
    gpcp_ts <- apply(chnks_gpcp[[i]][], 2, mean, na.rm = TRUE)

  }

  return(list(sst = sst_ts,
              prcp = gpcp_ts))
})

identical(ts_lst[[1]]$prcp, ts_lst[[2]]$prcp, ts_lst[[3]]$prcp)

ts1 <- ts_lst[[1]]$sst
ts2 <- ts_lst[[2]]$sst
ts3 <- ts_lst[[3]]$sst
ts4 <- ts_lst[[4]]$sst
rsp <- ts_lst[[1]]$prcp

library(forecast)

ts_reg <- lapply(seq(ts_lst[1:4]), function(i) {
  print(i)
  return(list(fit = auto.arima(rsp, xreg=ts_lst[[i]]$sst),#
              lm1 = lm(rsp ~ ts_lst[[i]]$sst)))
})

summary(ts_reg[[1]]$lm1)
summary(ts_reg[[2]]$lm1)
summary(ts_reg[[3]]$lm1)
summary(ts_reg[[4]]$lm1)
###########################################################################



##### MODE 2 ##############################################################
### load modes
# mds <- lapply(seq(patterns), function(i) {
#   do.call("stack", lapply(seq(nrow(patterns[[i]])), function(j) {
#     ptrn <- paste("*", paste("kili", patterns[[i]]$chunk[j],
#                              patterns[[i]]$lag[j], sep = "_"),
#                   ".rda", sep = "")
#     print(ptrn)
#     fl <- list.files("five_yr_chunks/results/modes",
#                      pattern = ptrn, full.names = TRUE)
#     load(fl)
#     modes@modes$mode_02@rsq_predictor
#   }))
# })
#
#
#
# plts <- lapply(seq(mds), function(i) {
#   lst <- lapply(seq(nlayers(mds[[i]])), function(j) {
#     spplot(mds[[i]][[j]], col.regions = envinmrPalette(1000),
#            at = seq(0, 1, 0.01))
#   })
#
#   return(latticeCombineGrid(lst))
# })
#
# plts[[1]]
#
# ext <- extent(c(27.5, 62.5, -17.5, 17.5))
#
# plts_mn <- lapply(seq(mds), function(i) {
#   spplot(mean(mds[[i]]), col.regions = envinmrPalette(1000),
#          at = seq(0, 1, 0.01),
#          colorkey = list(space = "left", height = 1/6*2, x = 0)) +
#     layer(panel.map("world", border = "grey40")) +
#     layer(panel.extent(extent(ext), border = "black", lwd = 2, lty = 1)) +
#     as.layer(spplot(mean(mds[[i]]), col.regions = "transparent",
#                     at = seq(0.3, 1, 0.7), contour = TRUE))
# })
#
# plts_mn[[1]]
#
# latticeCombineGrid(plts_mn)
#
#
#
# ### load data
# load("five_yr_chunks/data/stck_sst_dsn_dns80.rda")
# load("five_yr_chunks/data/stck_gpcp_dsn_dns_kili80.rda")
#
# chnks_sst <- lapply(seq(patterns), function(i) {
#   st <- substr(patterns[[i]]$chunk, 1, 4)
#   nd <- as.character(as.numeric(st) + 4)
#   sq <- seq(min(as.numeric(st)), max(as.numeric(nd)), 1)
#   print(sq)
#   ptrn <- paste("*", sq, "*", ".rst", sep = "")
#   #   fls <- do.call("c", lapply(ptrn, function(j) {
#   #     list.files("source_data/sst_1982_2010",
#   #                pattern = glob2rx(j),
#   #                full.names = TRUE)
#   #   }))
#   out <- stck_sst_dsn_dns[[which(substr(names(stck_sst_dsn_dns),
#                                         5, 8) %in% sq)]]
#   print(names(out))
#   return(out)
# })
#
# chnks_gpcp <- lapply(seq(patterns), function(i) {
#   st <- substr(patterns[[i]]$chunk, 1, 4)
#   nd <- as.character(as.numeric(st) + 4)
#   sq <- seq(min(as.numeric(st)), max(as.numeric(nd)), 1)
#   print(sq)
#   ptrn <- paste("*", sq, "*", ".rst", sep = "")
#   #   fls <- do.call("c", lapply(ptrn, function(j) {
#   #     list.files("source_data/sst_1982_2010",
#   #                pattern = glob2rx(j),
#   #                full.names = TRUE)
#   #   }))
#   out <- stck_gpcp_dsn_dns_kili[[which(substr(names(stck_gpcp_dsn_dns_kili),
#                                               6, 9) %in% sq)]]
#   print(names(out))
#   return(out)
# })
#
# mds_cln_mask <- lapply(seq(mds), function(i) {
#   out <- mean(mds[[i]])
#   out[out < 0.25] <- NA
#   out[out >= 0.25] <- 1
#   return(out)
# })
#
#
# ts_lst <- lapply(seq(patterns), function(i) {
#   cells <- which(mds_cln_mask[[i]][] == 1)
#
#   sst_ts <- apply(X = extract(chnks_sst[[i]], cells), MARGIN = 2,
#                   FUN = mean, na.rm = TRUE)
#   gpcp_ts <- apply(chnks_gpcp[[i]][], 2, mean, na.rm = TRUE)
#   x_by <- max(as.numeric(substr(patterns[[i]]$lag, 4, 5)))
#   print(x_by)
# #   if (i == 6) {
# #     st <- as.numeric(which(sst_ts == min(sst_ts)))
# #   } else {
# #     st <- as.numeric(which(abs(sst_ts) == max(abs(sst_ts))))
# #   }
#   sst_ts_srt <- sort(abs(sst_ts), decreasing = TRUE)
#   st1 <- as.numeric(which(abs(sst_ts) == sst_ts_srt[1]))
#   st2 <- as.numeric(which(abs(sst_ts) == sst_ts_srt[2]))
#   nd1 <- st1 + x_by
#   nd2 <- st2 + x_by
#   labs_st <- substr(names(sst_ts)[1], 5, 8)
#   labs_nd <- substr(names(sst_ts)[length(sst_ts)], 5, 8)
#   xyplot(sst_ts ~ seq(sst_ts), type = "l", lwd = 2,
#          col = "black", ylim = c(-2.2, 2.2),
#          panel = function(...) {
#            if (x_by > 0) {
#              panel.polygon(x = c(st1, st1, nd1, nd1),
#                            y = c(-2.2, 2.2, 2.2, -2.2),
#                            col = "grey90", border = "transparent")
#              panel.polygon(x = c(st2, st2, nd2, nd2),
#                            y = c(-2.2, 2.2, 2.2, -2.2),
#                            col = "grey90", border = "transparent")
#            } else {
#              panel.abline(v = c(st1, st2), lwd = 3, col = "grey90")
#            }
#            #            panel.abline(v = seq(st_ovr, nd_ovr, 1), lty = 1,
#            #                         lwd = 0.5, col = "grey")
#            panel.xyplot(...)
#            panel.text(x = rep(-1, 3),
#                       y = c(-1, 0, 1),
#                       labels = c("-1", "0", "1"))
#            panel.text(x = 0, y = 1.8, adj = c(0, 0.5),
#                       labels = paste(labs_st, labs_nd, sep = " - "))
#          }) +
#     as.layer(xyplot(gpcp_ts ~ seq(gpcp_ts),
#                     type = "l", lty = 3, lwd = 2)) +
#     layer_(panel.abline(h = 0, lty = 2, col = "grey")) +
#     layer_(panel.abline(h = 1, lty = 1, col = "grey")) +
#     layer_(panel.abline(h = -1, lty = 1, col = "grey"))
# })
#
#
#
# latticeCombineGrid(list(plts_mn[[1]], ts_lst[[1]],
#                         plts_mn[[2]], ts_lst[[2]],
#                         plts_mn[[3]], ts_lst[[3]],
#                         plts_mn[[4]], ts_lst[[4]],
#                         plts_mn[[5]], ts_lst[[5]],
#                         plts_mn[[6]], ts_lst[[6]]),
#                    layout = c(2, 6))
