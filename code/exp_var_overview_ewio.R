library(raster)
library(Reot)
library(reshape2)
library(latticeExtra)
library(Rsenal)

path_local <- "/media/ede/tims_ex/eot_kili_analysis"
setwd(path_local)

fls <- list.files("five_yr_chunks/results/modes",
                  pattern = glob2rx("*ewio*.rda"),
                  full.names = TRUE)
fls_nms <- list.files("five_yr_chunks/results/modes",
                      pattern = glob2rx("*ewio*.rda"),
                      full.names = FALSE)

exp_var_eot1 <- do.call("c", lapply(fls, function(i) {
  load(i)
  exv <- modes@modes$mode_02@cum_exp_var
  return(exv)
  rm(modes)
}))

chunks <- substr(fls_nms, 21, 29)
lags <- paste("lag",
              sprintf("%02.f", rep(0:12, length(unique(chunks)))),
              sep = "")

df <- data.frame(exp_var = exp_var_eot1,
                 lag = lags,
                 chunk = chunks)
df$chunk <- factor(df$chunk, levels = rev(levels(df$chunk)))

mat <- acast(df, chunk ~ lag, value.var = "exp_var")

# exp_var_eot1_mat <- matrix(exp_var_eot1, ncol = 13, nrow = 24,
#                            byrow = TRUE,
#                            dimnames = list(rev(unique(chunks)),
#                                            unique(lags)))

mx <- round(max(mat), 2)

plt <- levelplot(t(mat), col.regions = envinmrPalette(1000),
                 ylab = "Chunk", xlab = "Lag",
                 main = "Exp. var EWIO",
                 scales = list(x = list(rot = 45)),
                 at = seq(0, mx, length.out = 1000)) +
  as.layer(contourplot(t(mat), at = seq(0.15, mx, 0.05)))

png("five_yr_chunks/graphs/exp_var_overview_ewio_md2.png",
    height = 20, width = 15, units = "cm", res = 300)
plt
dev.off()

out <- list(mat = mat, plt = plt)

save(out, file = "five_yr_chunks/results/exp_var_overview_dns80_ewio_md2.rda")
