library(raster)
library(remote)
library(reshape2)
library(latticeExtra)
library(Rsenal)

path_local <- "/media/ede/tims_ex/eot_kili_analysis"
setwd(path_local)

fls <- list.files("five_yr_chunks/results/modes",
                  pattern = glob2rx("*kili*.rda"),
                  full.names = TRUE)
fls <- fls[-c(14:26)]
fls_nms <- list.files("five_yr_chunks/results/modes",
                      pattern = glob2rx("*kili*.rda"),
                      full.names = FALSE)
fls_nms <- fls_nms[-c(14:26)]

exp_var_eot <- do.call("c", lapply(fls, function(i) {
  load(i)
  exv <- modes@modes$mode_01@cum_exp_var
  return(exv)
  rm(modes)
}))

r_mn <- do.call("c", lapply(fls, function(i) {
  load(i)
  print(i)
  rmean <- cellStats(modes@modes$mode_01@r_response, mean, na.rm = TRUE)
  return(rmean)
  rm(modes)
}))

chunks <- substr(fls_nms, 25, 33)
lags <- paste("lag",
              sprintf("%02.f", rep(0:12, length(unique(chunks)))),
              sep = "")

df <- data.frame(exp_var = exp_var_eot,
                 r_mean = r_mn,
                 lag = lags,
                 chunk = chunks)
df$chunk <- factor(df$chunk, levels = rev(levels(df$chunk)))

mat <- acast(df, chunk ~ lag, value.var = "exp_var")
mat_r <- acast(df, chunk ~ lag, value.var = "r_mean")

# exp_var_eot_mat <- matrix(exp_var_eot, ncol = 13, nrow = 24,
#                            byrow = TRUE,
#                            dimnames = list(rev(unique(chunks)),
#                                            unique(lags)))

mx <- round(max(mat), 2) + 0.05

plt <- levelplot(t(mat), col.regions = envinmrPalette(1000),
                 ylab = "Chunk", xlab = "Lag",
                 main = "Exp. var KILI",
                 scales = list(x = list(rot = 45)),
                 at = seq(0, mx, length.out = 1000)) +
  as.layer(contourplot(t(mat), at = seq(0, mx, 0.25)))

clr_r <- colorRampPalette(brewer.pal(9, "RdBu"))
plt_r <- levelplot(t(mat_r), col.regions = clr_r(1000),
                   ylab = "Chunk", xlab = "Lag",
                   main = "Mean correlation between bp and response",
                   scales = list(x = list(rot = 45)),
                   at = seq(-1, 1, length.out = 1000))

out <- list(mat = mat, mat_r = mat_r, plt = plt, plt_r = plt_r)

save(out, file = "five_yr_chunks/results/exp_var_overview_dns80_kili.rda")

load("five_yr_chunks/results/exp_var_overview_dns80_kili.rda")

png("five_yr_chunks/graphs/exp_var_overview_kili.png",
    height = 20, width = 15, units = "cm", res = 300)
out$plt
dev.off()
