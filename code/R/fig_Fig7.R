# install necessary packages 
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# specify path to output files 
path_output_all_or_none <- '../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/efficacy/'

# list all of the output files 
files_output_all_or_none <- list.files(path = path_output_all_or_none, full.names = T, pattern = 'efficacy')

# load the files 
output_all_or_none <- lapply(files_output_all_or_none, function(ff){read.csv(ff)})
output_efficacy_all_or_none <- as.data.frame((do.call('rbind', output_all_or_none)))

# get the unique eir and sig_het values 
eir_equil <- sort(unique(output_efficacy_all_or_none$eir_equil))

# get the efficacies 
eff_cph_recurrent_LM_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_cph_recurrent_LM[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                           probs = c(0.25, 0.50, 0.75))})

eff_cph_recurrent_PCR_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_cph_recurrent_PCR[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                            probs = c(0.25, 0.50, 0.75))})

eff_cph_recurrent_D_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_cph_recurrent_D[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                              probs = c(0.25, 0.50, 0.75))})

eff_incid_recurrent_PCR_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_incid_trial_PCR[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                              probs = c(0.25, 0.50, 0.75))})

eff_incid_recurrent_LM_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_incid_trial_LM[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                             probs = c(0.25, 0.50, 0.75))})

eff_incid_recurrent_D_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_incid_trial_D[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                            probs = c(0.25, 0.50, 0.75))})

eff_risk_recurrent_PCR_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_risk_trial_PCR[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                             probs = c(0.25, 0.50, 0.75))})

eff_risk_recurrent_LM_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_risk_trial_LM[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                            probs = c(0.25, 0.50, 0.75))})

eff_risk_recurrent_D_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_risk_trial_D[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                           probs = c(0.25, 0.50, 0.75))})

# generate plot 
palette <- brewer.pal(n = 9, name = 'YlOrRd')
palette <- palette[c(5,7,9)]
offset <- seq(from = -0.45, to = 0.45, length.out = 9)
cex = 1.5
rect_width <- c(-0.5, mean(offset[3:4]), mean(offset[6:7]), 0.5)
eff_label <- c('Cox', 'Incidence', 'Risk')

jpeg(filename = '../../output/figs/fig_7.jpg', width = 8, height = 5, units = 'in', res = 500)
par(mar = c(3.6, 3.6, 1.2, 0.8))
plot(NA, NA, type = 'n', xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0,1),
     xaxs = 'i', yaxs = 'i', axes = F, xlab = '', ylab = '')
for(ee in 1:length(eir_equil))
{
  for(ss in 1:(length(rect_width) - 1))
  {
    rect(xleft = ee + rect_width[ss], xright = ee + rect_width[ss+1],
         ybottom = 0, ytop = 1, col = ifelse(ss %% 2 == 0, col2alpha('grey'), 'white'), border = NA)
    text(x = ee + mean(rect_width[ss:(ss+1)]), y = 0.975, eff_label[ss], cex = 0.75)
  }
}
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1))
abline(h = 0.75, lwd = 1, lty = 2, col = '#222222')
for(ee in 1:length(eir_equil))
{
  # cph metrics 
  segments(x0 = ee + offset[1], y0 = eff_cph_recurrent_PCR_all_or_none[1,ee], y1 = eff_cph_recurrent_PCR_all_or_none[3,ee],
           pch = 1, col = palette[1], lwd = cex)
  points(ee + offset[1], eff_cph_recurrent_PCR_all_or_none[2,ee], pch = 16, col = palette[1], cex = cex)
  
  
  segments(x0 = ee + offset[2], y0 = eff_cph_recurrent_LM_all_or_none[1,ee], y1 = eff_cph_recurrent_LM_all_or_none[3,ee],
           pch = 1, col = palette[2], lwd = cex)
  points(ee + offset[2], eff_cph_recurrent_LM_all_or_none[2,ee], pch = 16, col = palette[2], cex = cex)
  
  segments(x0 = ee + offset[3], y0 = eff_cph_recurrent_D_all_or_none[1,ee], y1 = eff_cph_recurrent_D_all_or_none[3,ee],
           pch = 1, col = palette[3], lwd = cex)
  points(ee + offset[3], eff_cph_recurrent_D_all_or_none[2,ee], pch = 16, col = palette[3], cex = cex)
  
  # incidence 
  segments(x0 = ee + offset[4], y0 = eff_incid_recurrent_PCR_all_or_none[1,ee], y1 = eff_incid_recurrent_PCR_all_or_none[3,ee],
           pch = 1, col = palette[1], lwd = cex)
  points(ee + offset[4], eff_incid_recurrent_PCR_all_or_none[2,ee], pch = 16, col = palette[1], cex = cex)
  
  segments(x0 = ee + offset[5], y0 = eff_incid_recurrent_LM_all_or_none[1,ee], y1 = eff_incid_recurrent_LM_all_or_none[3,ee],
           pch = 1, col = palette[2], lwd = cex)
  points(ee + offset[5], eff_incid_recurrent_LM_all_or_none[2,ee], pch = 16, col = palette[2], cex = cex)
  
  segments(x0 = ee + offset[6], y0 = eff_incid_recurrent_D_all_or_none[1,ee], y1 = eff_incid_recurrent_D_all_or_none[3,ee],
           pch = 1, col = palette[3], lwd = cex)
  points(ee + offset[6], eff_incid_recurrent_D_all_or_none[2,ee], pch = 16, col = palette[3], cex = cex)
  
  # risk measures 
  segments(x0 = ee + offset[7], y0 = eff_risk_recurrent_PCR_all_or_none[1,ee], y1 = eff_risk_recurrent_PCR_all_or_none[3,ee],
           pch = 1, col = palette[1], lwd = cex)
  points(ee + offset[7], eff_risk_recurrent_PCR_all_or_none[2,ee], pch = 16, col = palette[1], cex = cex)
  
  segments(x0 = ee + offset[8], y0 = eff_risk_recurrent_LM_all_or_none[1,ee], y1 = eff_risk_recurrent_LM_all_or_none[3,ee],
           pch = 1, col = palette[2], lwd = cex)
  points(ee + offset[8], eff_risk_recurrent_LM_all_or_none[2,ee], pch = 16, col = palette[2], cex = cex)
  
  segments(x0 = ee + offset[9], y0 = eff_risk_recurrent_D_all_or_none[1,ee], y1 = eff_risk_recurrent_D_all_or_none[3,ee],
           pch = 1, col = palette[3], lwd = cex)
  points(ee + offset[9], eff_risk_recurrent_D_all_or_none[2,ee], pch = 16, col = palette[3], cex = cex)
  
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Efficacy')

par(xpd = T)
legend(x = mean(c(0.5, length(eir_equil)+0.5)), y = 1.05, pch = c(rep(15, length(palette)),NA),
       lwd = c(rep(NA, length(palette)),1), col = c(palette, '#222222'),
       lty = c(rep(NA, length(palette)),2),
       legend = c('PCR-Detectable', 'LM-Detectable', 'Clinical', 'Clearance Probability'),
       bty = 'n', pt.cex = 1, ncol = 4, xjust = 0.5, cex = 0.6,
       text.width = 0.45)

dev.off()


