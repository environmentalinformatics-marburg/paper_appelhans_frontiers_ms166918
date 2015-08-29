library(raster)
library(remote)
library(reshape2)
library(latticeExtra)
library(Rsenal)
library(dplyr)
library(grid)

path_local <- "/media/ede/tims_ex/eot_kili_analysis"
setwd(path_local)

load("five_yr_chunks/results/exp_var_overview_dns80_ewio.rda")

df <- melt(out$mat)
names(df) <- c("chunk", "lag", "value")

### > 0.15
df15 <- data.frame(chunk = df$chunk[df$value >= 0.15])
df15$lag <- df$lag[df$value >= 0.15]

ind_chunk_15 <- c("2002-2006", "2003-2007", "2004-2008", "2005-2009")
ind_lag_15 <- c("lag04", "lag05", "lag06", "lag07")
df15 <- df15[df15$chunk %in% ind_chunk_15, ]
df15 <- df15[!df15$lag %in% ind_lag_15, ]


### > 0.25
df25 <- data.frame(chunk = df$chunk[df$value >= 0.25])
df25$lag <- df$lag[df$value >= 0.25]

ind_chunk_25 <- c("1993-1997", "1994-1998", "1995-1999", "1996-2000")
ind_lag_25 <- c("lag03", "lag04", "lag05", "lag06", "lag07", "lag08")
df25 <- df25[df25$chunk %in% ind_chunk_25, ]
df25 <- df25[!df25$lag %in% ind_lag_25, ]

patterns <- list(df25[1:3, ],
                 df25[4:5, ],
                 df25[6:7, ],
                 df15[1:3, ],
                 df15[4:9, ],
                 df15[17:20, ])
###########################################################################



##### MODE 1 ##############################################################
### load modes
mds <- lapply(seq(patterns), function(i) {
  do.call("stack", lapply(seq(nrow(patterns[[i]])), function(j) {
    ptrn <- paste("*", paste("ewio", patterns[[i]]$chunk[j],
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

ext <- extent(c(27.5, 62.5, -17.5, 17.5))

plts_mn <- lapply(seq(mds), function(i) {
  spplot(mean(mds[[i]]), col.regions = envinmrPalette(1000),
         at = seq(0, 1, 0.01),
         colorkey = FALSE) + #list(space = "left", height = 1/6*2, x = 0)) +
    layer(panel.map("world", border = "grey40", lwd = 0.5)) +
    layer(panel.extent(extent(ext), border = "black", lwd = 2, lty = 1)) +
    as.layer(spplot(mean(mds[[i]]), col.regions = "transparent",
                     at = seq(0.65, 1, 0.35), contour = TRUE))
})

plts_mn[[1]]

latticeCombineGrid(plts_mn)



### load data
load("five_yr_chunks/data/stck_sst_dsn_dns80.rda")
load("five_yr_chunks/data/stck_gpcp_dsn_dns_ewio80.rda")

chnks_sst <- lapply(seq(patterns), function(i) {
  st <- substr(patterns[[i]]$chunk, 1, 4)
  nd <- as.character(as.numeric(st) + 4)
  sq <- seq(min(as.numeric(st)), max(as.numeric(nd)), 1)
  print(sq)
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
  st <- substr(patterns[[i]]$chunk, 1, 4)
  nd <- as.character(as.numeric(st) + 4)
  sq <- seq(min(as.numeric(st)), max(as.numeric(nd)), 1)
  print(sq)
  ptrn <- paste("*", sq, "*", ".rst", sep = "")
  #   fls <- do.call("c", lapply(ptrn, function(j) {
  #     list.files("source_data/sst_1982_2010",
  #                pattern = glob2rx(j),
  #                full.names = TRUE)
  #   }))
  out <- stck_gpcp_dsn_dns_ewio[[which(substr(names(stck_gpcp_dsn_dns_ewio),
                                        6, 9) %in% sq)]]
  print(names(out))
  return(out)
})

mds_cln_mask <- lapply(seq(mds), function(i) {
  out <- mean(mds[[i]])
  out[out < 0.65] <- NA
  out[out >= 0.65] <- 1
  return(out)
})


ts_lst <- lapply(seq(patterns), function(i) {
  cells <- which(mds_cln_mask[[i]][] == 1)

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
         col = "black", ylim = c(-2.2, 2.2),
         panel = function(...) {
           if (x_by > 0) {
             panel.polygon(x = c(st1, st1, nd1, nd1),
                           y = c(-2.2, 2.2, 2.2, -2.2),
                           col = "grey90", border = "transparent")
             panel.polygon(x = c(st2, st2, nd2, nd2),
                           y = c(-2.2, 2.2, 2.2, -2.2),
                           col = "grey90", border = "transparent")
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
    as.layer(xyplot(gpcp_ts ~ seq(gpcp_ts),
                    type = "l", lty = 3, lwd = 2)) +
    layer_(panel.abline(h = 0, lty = 2, col = "grey")) +
    layer_(panel.abline(h = 1, lty = 1, col = "grey")) +
    layer_(panel.abline(h = -1, lty = 1, col = "grey"))
})


final <- latticeCombineGrid(list(plts_mn[[1]], ts_lst[[1]],
                                 plts_mn[[2]], ts_lst[[2]],
                                 plts_mn[[3]], ts_lst[[3]],
                                 plts_mn[[4]], ts_lst[[4]],
                                 plts_mn[[5]], ts_lst[[5]],
                                 plts_mn[[6]], ts_lst[[6]]),
                            layout = c(2, 6))

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


###########################################################################



##### MODE 2 ##############################################################
### load modes
# mds <- lapply(seq(patterns), function(i) {
#   do.call("stack", lapply(seq(nrow(patterns[[i]])), function(j) {
#     ptrn <- paste("*", paste("ewio", patterns[[i]]$chunk[j],
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
# load("five_yr_chunks/data/stck_gpcp_dsn_dns_ewio80.rda")
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
#   out <- stck_gpcp_dsn_dns_ewio[[which(substr(names(stck_gpcp_dsn_dns_ewio),
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
