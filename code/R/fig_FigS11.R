# install necessary packages
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# specify path to output files 
path_output_all_or_none <- '../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/model/'

# specify parameter ranges
sigma_het <- seq(from = 0, to = 3, by = 1)
EIR_equil <- c(1, 10, 100)
n_rep <- 200

# construct parameter grid 
param_sweep <- expand.grid(sigma_het = sigma_het,
                           EIR_equil = EIR_equil)

param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each  = n_rep), ]

# get the file indices for homogeneous biting scenarios 
file_indices <- lapply(EIR_equil, function(eir){which(param_sweep$EIR_equil == eir & param_sweep$sigma_het == 0)-1})

# load the necessary files 
EIR_over_time <- list()
for(ii in 1:length(file_indices))
{
  print(ii)
  mat <- matrix(NA, nrow = length(file_indices[[ii]]), ncol = 4015)
  for(ff in 1:length(file_indices[[ii]]))
  {
    print(ff)
    file <- paste(path_output_all_or_none, 'model_', file_indices[[ii]][ff], '.txt.bz2', sep = '')
    if(file.exists(file))
    {
      model_pred <- read.table(file)
      mat[ff,] <- model_pred[,29] * 365
    }
  }
  EIR_over_time[[ii]] <- mat
}

# generate plots 
time <- model_pred[,1]/365

jpeg(filename = '../../output/figs/fig_S11.jpg', width = 5, height = 8, units = 'in', res = 500)
layout(mat = matrix(1:3, nrow = 3))
par(mar = c(3.3,3.9,1.6,0.8))
plot(NA, NA, xlim = c(min(time), max(time)+0.5), ylim = c(0.8,1.5), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(v = 2025, lty = 2, col = '#222222')
EIR_CI <- apply(EIR_over_time[[1]], 2, quantile, probs = c(0,0.25,0.50,0.75,1), na.rm = T)
polygon(c(time, rev(time)), c(EIR_CI[1,], rev(EIR_CI[5,])),
        col = col2alpha('blue', alpha = 0.25), border = NA)
polygon(c(time, rev(time)), c(EIR_CI[2,], rev(EIR_CI[4,])),
        col = col2alpha('blue', alpha = 0.5), border = NA)
lines(time, EIR_CI[3,], col = 'blue', lwd = 1.5)
box()
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 2, line = 2.6, 'EIR (ibppy)')
mtext(side = 3, line = 0, bquote('Equilibrium EIR' == 1))
mtext(side = 3, line = 0, 'A', adj = 0, font = 2)
legend('topright', lty = c(2,1,NA,NA), col = c('#222222', 'blue', col2alpha('blue', alpha = 0.5), col2alpha('blue', alpha = 0.25)), 
       legend = c('Start of Trial Enrollment', 'Median', 'IQR', 'Range'), pch = c(NA,NA,15,15), pt.cex = 2, lwd = 1.5, bty = 'n')

plot(NA, NA, xlim = c(min(time), max(time)+0.5), ylim = c(9,11), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(v = 2025, lty = 2, col = '#222222')
EIR_CI <- apply(EIR_over_time[[2]], 2, quantile, probs = c(0,0.25,0.50,0.75,1), na.rm = T)
polygon(c(time, rev(time)), c(EIR_CI[1,], rev(EIR_CI[5,])),
        col = col2alpha('blue', alpha = 0.25), border = NA)
polygon(c(time, rev(time)), c(EIR_CI[2,], rev(EIR_CI[4,])),
        col = col2alpha('blue', alpha = 0.5), border = NA)
lines(time, EIR_CI[3,], col = 'blue', lwd = 1.5)
box()
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 2, line = 2.6, 'EIR (ibppy)')
mtext(side = 3, line = 0, bquote('Equilibrium EIR' == 10))
mtext(side = 3, line = 0, 'B', adj = 0, font = 2)

plot(NA, NA, xlim = c(min(time), max(time)+0.5), ylim = c(95,105), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(v = 2025, lty = 2, col = '#222222')
EIR_CI <- apply(EIR_over_time[[3]], 2, quantile, probs = c(0,0.25,0.50,0.75,1), na.rm = T)
polygon(c(time, rev(time)), c(EIR_CI[1,], rev(EIR_CI[5,])),
        col = col2alpha('blue', alpha = 0.25), border = NA)
polygon(c(time, rev(time)), c(EIR_CI[2,], rev(EIR_CI[4,])),
        col = col2alpha('blue', alpha = 0.5), border = NA)
lines(time, EIR_CI[3,], col = 'blue', lwd = 1.5)
box()
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Time')
mtext(side = 2, line = 2.6, 'EIR (ibppy)')
mtext(side = 3, line = 0, bquote('Equilibrium EIR' == 100))
mtext(side = 3, line = 0, 'C', adj = 0, font = 2)

dev.off()