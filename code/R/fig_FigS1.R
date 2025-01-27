if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# specify path to files 
path_in = '../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/survival/'

# specify parameter ranges
sigma_het <- seq(from = 0, to = 3, by = 1)
EIR_equil <- c(1, 10, 100)
n_rep <- 200

# construct parameter grid 
param_sweep <- expand.grid(sigma_het = sigma_het,
                           EIR_equil = EIR_equil)

param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each  = n_rep), ]

# loop through and load survival curves 
surv_placebo_PCR <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_treatment_PCR <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_placebo_LM <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_treatment_LM <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_placebo_D <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_treatment_D <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)

for(ii in 1:nrow(param_sweep))
{
  print(ii)
  file_surv <- paste(path_in, 'survival_', ii - 1, '.csv', sep = '')
  if(file.exists(file_surv))
  {
    surv <- read.csv(file_surv)
    surv_placebo_PCR[ii,] <- surv$surv_placebo_PCR
    surv_treatment_PCR[ii,] <- surv$surv_treatment_PCR
    
    surv_placebo_LM[ii,] <- surv$surv_placebo_LM
    surv_treatment_LM[ii,] <- surv$surv_treatment_LM
    
    surv_placebo_D[ii,] <- surv$surv_placebo_D
    surv_treatment_D[ii,] <- surv$surv_treatment_D
  }
}

# generate plot
palette <- brewer.pal(n = 9, name = 'YlOrRd')
palette <- palette[c(5,7,9)]
times <- 0:(ncol(surv_placebo_PCR)-1)

jpeg(filename = '../../output/figs/fig_S1.jpg', width = 10, height = 7.5, units = 'in', res = 500)
par(mar = c(3.3, 5.3,2.1,0.8))
layout(mat = matrix(1:12, nrow = 3, ncol = 4, byrow = T))
for(ee in 1:length(EIR_equil))
{
  for(ss in 1:length(sigma_het))
  {
    indices <- which(param_sweep$sigma_het == sigma_het[ss] & param_sweep$EIR_equil == EIR_equil[ee])
    plot(NA, NA, xlim = c(0, ncol(surv_placebo_PCR) - 1), ylim = c(0,1), axes = F,
         xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
    
    surv_placebo_PCR_IQR <- apply(surv_placebo_PCR[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_placebo_PCR_IQR[1,], rev(surv_placebo_PCR_IQR[3,])),
            col = col2alpha(palette[1]), border = NA)
    lines(surv_placebo_PCR_IQR[2,], col = palette[1], lty = 2)
    surv_treatment_PCR_IQR <- apply(surv_treatment_PCR[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_treatment_PCR_IQR[1,], rev(surv_treatment_PCR_IQR[3,])),
            col = col2alpha(palette[1]), border = NA)
    lines(surv_treatment_PCR_IQR[2,], col = palette[1], lty = 1)
    
    surv_placebo_LM_IQR <- apply(surv_placebo_LM[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_placebo_LM_IQR[1,], rev(surv_placebo_LM_IQR[3,])),
            col = col2alpha(palette[2]), border = NA)
    lines(surv_placebo_LM_IQR[2,], col = palette[2], lty = 2)
    surv_treatment_LM_IQR <- apply(surv_treatment_LM[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_treatment_LM_IQR[1,], rev(surv_treatment_LM_IQR[3,])),
            col = col2alpha(palette[2]), border = NA)
    lines(surv_treatment_LM_IQR[2,], col = palette[2], lty = 1)
    
    surv_placebo_D_IQR <- apply(surv_placebo_D[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_placebo_D_IQR[1,], rev(surv_placebo_D_IQR[3,])),
            col = col2alpha(palette[3]), border = NA)
    lines(surv_placebo_D_IQR[2,], col = palette[3], lty = 2)
    surv_treatment_D_IQR <- apply(surv_treatment_D[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_treatment_D_IQR[1,], rev(surv_treatment_D_IQR[3,])),
            col = col2alpha(palette[3]), border = NA)
    lines(surv_treatment_D_IQR[2,], col = palette[3], lty = 1)
    
    box()
    axis(side = 1)
    axis(side = 2, las = 1)
    if(ee == 1)
    {
      obj <- list(sigma = sigma_het[ss])
      mtext(side = 3, line = 0,  bquote(sigma^2 == .(obj$sigma)), font = 2)
    }
    if(ss == 1)
    {
      obj <- list(eir = EIR_equil[ee])
      mtext(side = 2, line = 3.6, bquote('EIR' == .(obj$eir)), font = 2)
    }
    if(ee == 3)
    {
      mtext(side = 1, line = 2.3, 'Follow-up (d)', cex = 0.8)
    }
    if(ss == 1)
    {
      mtext(side = 2, line = 2.3, 'Proportion Recurrence-Free', cex = 0.8)
    }
    if(ee == 3 & ss == 1)
    {
      legend('right', col = c('#222222', '#222222', rev(palette)), lwd = 1.5, lty = c(1,2,1,1,1), legend = c('Treatment', 'Placebo', 'Clinical',
                                                                                                             'LM-Detectable', 'PCR-Detectable'), bty = 'n')
    }
  }
}
dev.off()
