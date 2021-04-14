# install necessary packages
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# specify path to files 
path_in = '../../output/analysis/radical_cure_therapeutic/output_files/all_or_none/survival/'

# specify parameter ranges 
EIR_equil <- c(1 / 365, 10 / 365, 100 / 365)
trial_PQ_proph <- c(28, 45)

# construct parameter grid 
n_rep = 200
param_sweep <- expand.grid(EIR_equil = EIR_equil,
                           trial_PQ_proph = trial_PQ_proph)

param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each = n_rep), ]

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

jpeg(filename = '../../output/figs/fig_S9.jpg', width = 7.5, height = 5, units = 'in', res = 500)
par(mar = c(3.3, 5.3,2.1,0.8))
layout(mat = matrix(1:6, nrow = 2, ncol = 3, byrow = F))
for(ee in 1:length(EIR_equil))
{
  for(pp in 1:length(trial_PQ_proph))
  {
    indices <- which(param_sweep$trial_PQ_proph == trial_PQ_proph[pp] & param_sweep$EIR_equil == EIR_equil[ee])
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
      mtext(side = 2, line = 3.6,  paste('Prophylaxis: ', trial_PQ_proph[pp], 'd', sep = ''), font = 1)
    }
    if(pp == 1)
    {
      obj <- list(eir = EIR_equil[ee] * 365)
      mtext(side = 3, line = 0, bquote('EIR' == .(obj$eir)), font = 2)
    }
    if(pp == 2)
    {
      mtext(side = 1, line = 2.3, 'Follow-up (d)', cex = 0.8)
    }
    if(ee == 1)
    {
      mtext(side = 2, line = 2.3, 'Proportion Recurrence-Free', cex = 0.8)
    }
  }
}
dev.off()
