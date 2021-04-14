# install necessary packages
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# specify path to files 
path_in_all_or_none = '../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/survival/'
path_in_leaky = '../../output/analysis/eir_vs_heterogeneity/output_files/leaky/survival/'

# specify parameter ranges
sigma_het <- seq(from = 0, to = 3, by = 1)
EIR_equil <- c(1, 10, 100)
n_rep <- 200

# construct parameter grid 
param_sweep <- expand.grid(sigma_het = sigma_het,
                           EIR_equil = EIR_equil)

param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each  = n_rep), ]

# loop through and load survival curves 
surv_placebo_PCR_all_or_none <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_treatment_PCR_all_or_none <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_placebo_LM_all_or_none <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_treatment_LM_all_or_none <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_placebo_D_all_or_none <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_treatment_D_all_or_none <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)

surv_placebo_PCR_leaky <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_treatment_PCR_leaky <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_placebo_LM_leaky <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_treatment_LM_leaky <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_placebo_D_leaky <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)
surv_treatment_D_leaky <- matrix(NA, nrow = nrow(param_sweep), ncol = 181)

for(ii in 1:nrow(param_sweep))
{
  print(ii)
  file_surv <- paste(path_in_all_or_none, 'survival_', ii - 1, '.csv', sep = '')
  if(file.exists(file_surv))
  {
    surv <- read.csv(file_surv)
    surv_placebo_PCR_all_or_none[ii,] <- surv$surv_placebo_PCR
    surv_treatment_PCR_all_or_none[ii,] <- surv$surv_treatment_PCR
    
    surv_placebo_LM_all_or_none[ii,] <- surv$surv_placebo_LM
    surv_treatment_LM_all_or_none[ii,] <- surv$surv_treatment_LM
    
    surv_placebo_D_all_or_none[ii,] <- surv$surv_placebo_D
    surv_treatment_D_all_or_none[ii,] <- surv$surv_treatment_D
  }
  file_surv <- paste(path_in_leaky, 'survival_', ii - 1, '.csv', sep = '')
  if(file.exists(file_surv))
  {
    surv <- read.csv(file_surv)
    surv_placebo_PCR_leaky[ii,] <- surv$surv_placebo_PCR
    surv_treatment_PCR_leaky[ii,] <- surv$surv_treatment_PCR
    
    surv_placebo_LM_leaky[ii,] <- surv$surv_placebo_LM
    surv_treatment_LM_leaky[ii,] <- surv$surv_treatment_LM
    
    surv_placebo_D_leaky[ii,] <- surv$surv_placebo_D
    surv_treatment_D_leaky[ii,] <- surv$surv_treatment_D
  }
}

# generate plot 
times <- 0:(ncol(surv_placebo_PCR_all_or_none)-1)

jpeg(filename = '../../output/figs/fig_S10.jpg', width = 10, height = 7.5, units = 'in', res = 500)
par(mar = c(3.3, 5.3,2.1,0.8))
par(mar = c(3.3, 5.3,2.1,0.8))
layout(mat = matrix(1:12, nrow = 3, ncol = 4, byrow = T))
for(ee in 1:length(EIR_equil))
{
  for(ss in 1:length(sigma_het))
  {
    indices <- which(param_sweep$sigma_het == sigma_het[ss] & param_sweep$EIR_equil == EIR_equil[ee])
    plot(NA, NA, xlim = c(0, ncol(surv_placebo_LM_all_or_none) - 1), ylim = c(0,1), axes = F,
         xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
    
    surv_placebo_LM_IQR <- apply(surv_placebo_LM_all_or_none[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_placebo_LM_IQR[1,], rev(surv_placebo_LM_IQR[3,])),
            col = col2alpha('red'), border = NA)
    lines(surv_placebo_LM_IQR[2,], col = 'red', lty = 2)
    surv_treatment_LM_IQR <- apply(surv_treatment_LM_all_or_none[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_treatment_LM_IQR[1,], rev(surv_treatment_LM_IQR[3,])),
            col = col2alpha('red'), border = NA)
    lines(surv_treatment_LM_IQR[2,], col = 'red', lty = 1)
    
    surv_placebo_LM_IQR <- apply(surv_placebo_LM_leaky[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_placebo_LM_IQR[1,], rev(surv_placebo_LM_IQR[3,])),
            col = col2alpha('blue'), border = NA)
    lines(surv_placebo_LM_IQR[2,], col = 'blue', lty = 2)
    surv_treatment_LM_IQR <- apply(surv_treatment_LM_leaky[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_treatment_LM_IQR[1,], rev(surv_treatment_LM_IQR[3,])),
            col = col2alpha('blue'), border = NA)
    lines(surv_treatment_LM_IQR[2,], col = 'blue', lty = 1)
    
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
    if(ss == 1 & ee == 3)
    {
      legend('right', col = c('#222222', '#222222', 'red', 'blue'), lwd = 1.5, lty = c(1,2,1,1,1), legend = c('Treatment', 'Placebo', 'All-or-None', 'Leaky'), bty = 'n')
    }
  }
}
dev.off()