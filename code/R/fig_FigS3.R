# install necessary packages
if(!require(grDevices)){install.packages('grDevices'); library(grDevices)}
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# specify path to files 
path_in = '../../output/analysis/followup_vs_relapse/output_files/all_or_none/survival/'

# specify parameter ranges 
EIR_equil <- c(1 / 365, 10 / 365, 100 / 365)
trial_duration <- c(90, 180, 365, 730)
ff <- c(1/30, 1/60, 1/90, 1/180)

# construct parameter grid 
n_rep <- 200
param_sweep <- expand.grid(EIR_equil = EIR_equil,
                           trial_duration = trial_duration,
                           ff = ff)
param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each = n_rep), ]

# loop through and load survival curves 
surv_placebo_PCR <- matrix(NA, nrow = nrow(param_sweep), ncol = max(trial_duration)+1)
surv_treatment_PCR <- matrix(NA, nrow = nrow(param_sweep), ncol = max(trial_duration)+1)
surv_placebo_LM <- matrix(NA, nrow = nrow(param_sweep), ncol = max(trial_duration)+1)
surv_treatment_LM <- matrix(NA, nrow = nrow(param_sweep), ncol = max(trial_duration)+1)
surv_placebo_D <- matrix(NA, nrow = nrow(param_sweep), ncol = max(trial_duration)+1)
surv_treatment_D <- matrix(NA, nrow = nrow(param_sweep), ncol = max(trial_duration)+1)

for(ii in 1:nrow(param_sweep))
{
  print(ii)
  file_surv <- paste(path_in, 'survival_', ii - 1, '.csv', sep = '')
  if(file.exists(file_surv))
  {
    surv <- read.csv(file_surv)
    surv_placebo_PCR[ii,1:nrow(surv)] <- surv$surv_placebo_PCR
    surv_treatment_PCR[ii,1:nrow(surv)] <- surv$surv_treatment_PCR
    
    surv_placebo_LM[ii,1:nrow(surv)] <- surv$surv_placebo_LM
    surv_treatment_LM[ii,1:nrow(surv)] <- surv$surv_treatment_LM
    
    surv_placebo_D[ii,1:nrow(surv)] <- surv$surv_placebo_D
    surv_treatment_D[ii,1:nrow(surv)] <- surv$surv_treatment_D
  }
}

# generate plot
palette <- brewer.pal(n = 9, name = 'YlOrRd')
palette <- palette[c(5,7,9)]


jpeg(filename = '../../output/figs/fig_S3.jpg', width = 10, height = 7.5, units = 'in', res = 500)
par(mar = c(3.3, 4.8,1.6,0.8))
layout(mat = matrix(1:12, nrow = 3, ncol = 4, byrow = F))
for(rr in 1:length(ff))
{
  for(ee in 1:length(EIR_equil))
  {
    indices <- which(param_sweep$EIR_equil == EIR_equil[ee] & param_sweep$trial_duration == 730 & param_sweep$ff == ff[rr])
    plot(NA, NA, xlim = c(0, ncol(surv_placebo_LM) - 1), ylim = c(0,1), axes = F,
         xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
    
    surv_placebo_PCR_IQR <- apply(surv_placebo_PCR[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    times <- 0:(sum(!is.na(surv_placebo_PCR_IQR[1,]))-1)
    polygon(c(times, rev(times)), c(surv_placebo_PCR_IQR[1,times+1], rev(surv_placebo_PCR_IQR[3,times+1])),
            col = col2alpha(palette[1]), border = NA)
    lines(surv_placebo_PCR_IQR[2,times+1], col = palette[1], lty = 2)
    surv_treatment_PCR_IQR <- apply(surv_treatment_PCR[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_treatment_PCR_IQR[1,times+1], rev(surv_treatment_PCR_IQR[3,times+1])),
            col = col2alpha(palette[1]), border = NA)
    lines(surv_treatment_PCR_IQR[2,times+1], col = palette[1], lty = 1)
    
    surv_placebo_LM_IQR <- apply(surv_placebo_LM[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    times <- 0:(sum(!is.na(surv_placebo_LM_IQR[1,]))-1)
    polygon(c(times, rev(times)), c(surv_placebo_LM_IQR[1,times+1], rev(surv_placebo_LM_IQR[3,times+1])),
            col = col2alpha(palette[2]), border = NA)
    lines(surv_placebo_LM_IQR[2,times+1], col = palette[2], lty = 2)
    surv_treatment_LM_IQR <- apply(surv_treatment_LM[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_treatment_LM_IQR[1,times+1], rev(surv_treatment_LM_IQR[3,times+1])),
            col = col2alpha(palette[2]), border = NA)
    lines(surv_treatment_LM_IQR[2,times+1], col = palette[2], lty = 1)
    
    surv_placebo_D_IQR <- apply(surv_placebo_D[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    times <- 0:(sum(!is.na(surv_placebo_D_IQR[1,]))-1)
    polygon(c(times, rev(times)), c(surv_placebo_D_IQR[1,times+1], rev(surv_placebo_D_IQR[3,times+1])),
            col = col2alpha(palette[3]), border = NA)
    lines(surv_placebo_D_IQR[2,times+1], col = palette[3], lty = 2)
    surv_treatment_D_IQR <- apply(surv_treatment_D[indices,], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
    polygon(c(times, rev(times)), c(surv_treatment_D_IQR[1,times+1], rev(surv_treatment_D_IQR[3,times+1])),
            col = col2alpha(palette[3]), border = NA)
    lines(surv_treatment_D_IQR[2,times+1], col = palette[3], lty = 1)
    box()
    axis(side = 1)
    axis(side = 2, las = 1)
    if(rr == 1)
    {
      mtext(side = 2, line = 2.3, 'Proportion Recurrence-Free', cex = 0.8)
      obj <- list(eir = EIR_equil[ee] * 365)
      mtext(side = 2, line = 3.6, bquote('EIR' == .(obj$eir)), font = 2)
    }
    if(ee == 1)
    {
      mtext(side = 3, line = 0, paste('Mean Relapse Time: ', 1 / ff[rr], 'd', sep = ''), font = 1, cex = 0.8)
    }
    if(ee == length(EIR_equil))
    {
      mtext(side = 1, line = 2.3, 'Follow-up (d)', cex = 0.8)
    } 
    if(ee == 3 & rr == 1)
    {
      legend('right', col = c('#222222', '#222222', rev(palette)), lwd = 1.5, lty = c(1,2,1,1,1), legend = c('Treatment', 'Placebo', 'Clinical',
                                                                                                                'LM-Detectable', 'PCR-Detectable'), bty = 'n')
    }
  }
}
dev.off()