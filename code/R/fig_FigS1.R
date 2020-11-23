# set working directory
setwd('~/Dropbox/Radical_Cure_MOA/code/')

# clear existing workspace
rm(list = ls())

# install necessary packages 
if(!require(ggsci)){install.packages('ggsci'); library(ggsci)}

# specify path to output files 
path_output_leaky <- '../output/analysis_0930/eir_vs_heterogeneity/output_files/leaky/'
path_output_all_or_none <- '../output/analysis_0930/eir_vs_heterogeneity/output_files/all_or_none/'

# list all of the output files 
files_output_leaky <- list.files(path = path_output_leaky, full.names = T, pattern = 'efficacy')
files_output_all_or_none <- list.files(path = path_output_all_or_none, full.names = T, pattern = 'efficacy')

# load the files 
output_leaky <- lapply(files_output_leaky, function(ff){read.csv(ff)})
output_efficacy_leaky <- as.data.frame((do.call('rbind', output_leaky)))

output_all_or_none <- lapply(files_output_all_or_none, function(ff){read.csv(ff)})
output_efficacy_all_or_none <- as.data.frame((do.call('rbind', output_all_or_none)))

# get the unique eir and sig_het values 
eir_equil <- sort(unique(output_efficacy_all_or_none$eir_equil))

# get the efficacies 
eff_cph_recurrent_LM_leaky <- sapply(eir_equil, function(eir){quantile(output_efficacy_leaky$eff_cph_relapse_LM_all[output_efficacy_leaky$sig_het == 0 & output_efficacy_leaky$eir_equil == eir],
                                                                                    probs = c(0.25, 0.50, 0.75))})
eff_cph_recurrent_LM_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_cph_relapse_LM_all[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                                          probs = c(0.25, 0.50, 0.75))})

eff_cph_recurrent_PCR_leaky <- sapply(eir_equil, function(eir){quantile(output_efficacy_leaky$eff_cph_relapse_PCR_all[output_efficacy_leaky$sig_het == 0 & output_efficacy_leaky$eir_equil == eir],
                                                                       probs = c(0.25, 0.50, 0.75))})
eff_cph_recurrent_PCR_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_cph_relapse_PCR_all[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                             probs = c(0.25, 0.50, 0.75))})

eff_incid_relapse_PCR_leaky <- sapply(eir_equil, function(eir){quantile(output_efficacy_leaky$eff_incid_relapse_PCR[output_efficacy_leaky$sig_het == 0 & output_efficacy_leaky$eir_equil == eir],
                                                                        probs = c(0.25, 0.50, 0.75))})
eff_incid_relapse_PCR_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_incid_relapse_PCR[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                              probs = c(0.25, 0.50, 0.75))})

eff_incid_relapse_LM_leaky <- sapply(eir_equil, function(eir){quantile(output_efficacy_leaky$eff_incid_relapse_LM[output_efficacy_leaky$sig_het == 0 & output_efficacy_leaky$eir_equil == eir],
                                                                            probs = c(0.25, 0.50, 0.75))})
eff_incid_relapse_LM_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_incid_relapse_LM[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                                  probs = c(0.25, 0.50, 0.75))})

eff_incid_relapse_D_leaky <- sapply(eir_equil, function(eir){quantile(output_efficacy_leaky$eff_incid_relapse_D[output_efficacy_leaky$sig_het == 0 & output_efficacy_leaky$eir_equil == eir],
                                                                           probs = c(0.25, 0.50, 0.75))})
eff_incid_relapse_D_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_incid_relapse_D[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                                 probs = c(0.25, 0.50, 0.75))})

eff_risk_relapse_PCR_leaky <- sapply(eir_equil, function(eir){quantile(output_efficacy_leaky$eff_risk_relapse_PCR[output_efficacy_leaky$sig_het == 0 & output_efficacy_leaky$eir_equil == eir],
                                                                          probs = c(0.25, 0.50, 0.75))})
eff_risk_relapse_PCR_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_risk_relapse_PCR[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                                probs = c(0.25, 0.50, 0.75))})

eff_risk_relapse_LM_leaky <- sapply(eir_equil, function(eir){quantile(output_efficacy_leaky$eff_risk_relapse_LM[output_efficacy_leaky$sig_het == 0 & output_efficacy_leaky$eir_equil == eir],
                                                                           probs = c(0.25, 0.50, 0.75))})
eff_risk_relapse_LM_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_risk_relapse_LM[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                                 probs = c(0.25, 0.50, 0.75))})

eff_risk_relapse_D_leaky <- sapply(eir_equil, function(eir){quantile(output_efficacy_leaky$eff_risk_relapse_D[output_efficacy_leaky$sig_het == 0 & output_efficacy_leaky$eir_equil == eir],
                                                                           probs = c(0.25, 0.50, 0.75))})
eff_risk_relapse_D_all_or_none <- sapply(eir_equil, function(eir){quantile(output_efficacy_all_or_none$eff_risk_relapse_D[output_efficacy_all_or_none$sig_het == 0 & output_efficacy_all_or_none$eir_equil == eir],
                                                                                 probs = c(0.25, 0.50, 0.75))})

# generate plot 
palette <- pal_material(palette = 'deep-orange', n = 8)(8)[-1][c(1,2,4)]
offset <- seq(from = -0.45, to = 0.45, length.out = 18)
cex = 1.15
rect_width <- c(-0.5, mean(offset[6:7]), mean(offset[12:13]), 0.5)
eff_label <- c('Cox', 'Incidence', 'Risk')

jpeg(filename = '../output/figs_manuscript/fig_S1.jpg', width = 8, height = 5, units = 'in', res = 500)
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
  points(ee + offset[1], eff_cph_recurrent_PCR_all_or_none[2,ee], pch = 16, col = palette[1], bg = 'white', cex = cex)
  
  segments(x0 = ee + offset[2], y0 = eff_cph_recurrent_PCR_leaky[1,ee], y1 = eff_cph_recurrent_PCR_all_or_none[3,ee],
           pch = 1, col = palette[1], lwd = cex)
  points(ee + offset[2], eff_cph_recurrent_PCR_leaky[2,ee], pch = 17, col = palette[1], bg = 'white', cex = cex)
  
  segments(x0 = ee + offset[3], y0 = eff_cph_recurrent_LM_all_or_none[1,ee], y1 = eff_cph_recurrent_LM_all_or_none[3,ee],
           pch = 1, col = palette[2], lwd = cex)
  points(ee + offset[3], eff_cph_recurrent_LM_all_or_none[2,ee], pch = 16, col = palette[2], bg = 'white', cex = cex)
  
  segments(x0 = ee + offset[4], y0 = eff_cph_recurrent_LM_leaky[1,ee], y1 = eff_cph_recurrent_LM_all_or_none[3,ee],
           pch = 1, col = palette[2], lwd = cex)
  points(ee + offset[4], eff_cph_recurrent_LM_leaky[2,ee], pch = 17, col = palette[2], bg = 'white', cex = cex)
  
  # incidence 
  segments(x0 = ee + offset[7], y0 = eff_incid_relapse_PCR_all_or_none[1,ee], y1 = eff_incid_relapse_PCR_all_or_none[3,ee],
           pch = 1, col = palette[1], lwd = cex)
  points(ee + offset[7], eff_incid_relapse_PCR_all_or_none[2,ee], pch = 16, col = palette[1], cex = cex)
  
  segments(x0 = ee + offset[8], y0 = eff_incid_relapse_PCR_leaky[1,ee], y1 = eff_incid_relapse_PCR_leaky[3,ee],
           pch = 1, col = palette[1], lwd = cex)
  points(ee + offset[8], eff_incid_relapse_PCR_leaky[2,ee], pch = 17, col = palette[1], cex = cex)
  
  segments(x0 = ee + offset[9], y0 = eff_incid_relapse_LM_all_or_none[1,ee], y1 = eff_incid_relapse_LM_all_or_none[3,ee],
           pch = 1, col = palette[2], lwd = cex)
  points(ee + offset[9], eff_incid_relapse_LM_all_or_none[2,ee], pch = 16, col = palette[2], cex = cex)
  
  segments(x0 = ee + offset[10], y0 = eff_incid_relapse_LM_leaky[1,ee], y1 = eff_incid_relapse_LM_leaky[3,ee],
           pch = 1, col = palette[2], lwd = cex)
  points(ee + offset[10], eff_incid_relapse_LM_leaky[2,ee], pch = 17, col = palette[2], cex = cex)
  
  segments(x0 = ee + offset[11], y0 = eff_incid_relapse_D_all_or_none[1,ee], y1 = eff_incid_relapse_D_all_or_none[3,ee],
           pch = 1, col = palette[3], lwd = cex)
  points(ee + offset[11], eff_incid_relapse_D_all_or_none[2,ee], pch = 16, col = palette[3], cex = cex)
  
  segments(x0 = ee + offset[12], y0 = eff_incid_relapse_D_leaky[1,ee], y1 = eff_incid_relapse_D_leaky[3,ee],
           pch = 1, col = palette[3], lwd = cex)
  points(ee + offset[12], eff_incid_relapse_D_leaky[2,ee], pch = 17, col = palette[3], cex = cex)
  
  # risk measures 
  segments(x0 = ee + offset[13], y0 = eff_risk_relapse_PCR_all_or_none[1,ee], y1 = eff_risk_relapse_PCR_all_or_none[3,ee],
           pch = 1, col = palette[1], lwd = cex)
  points(ee + offset[13], eff_risk_relapse_PCR_all_or_none[2,ee], pch = 16, col = palette[1], cex = cex)

  segments(x0 = ee + offset[14], y0 = eff_risk_relapse_PCR_leaky[1,ee], y1 = eff_risk_relapse_PCR_leaky[3,ee],
           pch = 1, col = palette[1], lwd = cex)
  points(ee + offset[14], eff_risk_relapse_PCR_leaky[2,ee], pch = 17, col = palette[1], cex = cex)
  
  segments(x0 = ee + offset[15], y0 = eff_risk_relapse_LM_all_or_none[1,ee], y1 = eff_risk_relapse_LM_all_or_none[3,ee],
           pch = 1, col = palette[2], lwd = cex)
  points(ee + offset[15], eff_risk_relapse_LM_all_or_none[2,ee], pch = 16, col = palette[2], cex = cex)
  
  segments(x0 = ee + offset[16], y0 = eff_risk_relapse_LM_leaky[1,ee], y1 = eff_risk_relapse_LM_leaky[3,ee],
           pch = 1, col = palette[2], lwd = cex)
  points(ee + offset[16], eff_risk_relapse_LM_leaky[2,ee], pch = 17, col = palette[2], cex = cex)
  
  segments(x0 = ee + offset[17], y0 = eff_risk_relapse_D_all_or_none[1,ee], y1 = eff_risk_relapse_D_all_or_none[3,ee],
           pch = 1, col = palette[3], lwd = cex)
  points(ee + offset[17], eff_risk_relapse_D_all_or_none[2,ee], pch = 16, col = palette[3], cex = cex)
  
  segments(x0 = ee + offset[18], y0 = eff_risk_relapse_D_leaky[1,ee], y1 = eff_risk_relapse_D_leaky[3,ee],
           pch = 1, col = palette[3], lwd = cex)
  points(ee + offset[18], eff_risk_relapse_D_leaky[2,ee], pch = 17, col = palette[3], cex = cex)
  
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Efficacy')

par(xpd = T)
legend(x = mean(c(0.5, length(eir_equil)+0.5)) - 0.125, y = 1.05, pch = c(1,2,rep(15, length(palette)),NA),
       lwd = c(NA, NA, rep(NA, length(palette)),1), col = c('#222222', '#222222', palette, '#222222'),
       lty = c(NA, NA, rep(NA, length(palette)),2),
       legend = c('All-or-None', 'Leaky', 'PCR-Detectable', 'LM-Detectable', 'Symptomatic', 'Hypnozoite Clearance Prob'),
       bty = 'n', pt.cex = 1, ncol = 6, xjust = 0.5, cex = 0.6,
       text.width = 0.45)


#legend('bottomleft', pch = c(1,2,0, 15, 22, 15, 15, 15, NA),
#       lwd = c(NA, NA, NA, NA, NA, NA, NA, NA, 1), col = c('#222222', '#222222', 'grey', 'grey', '#222222', palette, '#222222'),
#       lty = c(NA, NA, NA, NA, NA, NA, NA, NA, 2), pt.bg = c(NA, NA, NA, NA, 'grey', NA, NA, NA, NA),
#       legend = c('All-or-None', 'Leaky', 'Cox', 'Incidence', 'Risk', 'PCR-Detectable', 'LM-Detectable', 'Symptomatic', 'Hypnozoite Clearance Prob'),
#       bty = 'n', pt.cex = 1.25, cex = 0.675)
dev.off()

