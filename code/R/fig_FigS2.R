# set working directory
setwd('~/Dropbox/Radical_Cure_MOA/code/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(ggsci)){install.packages('ggsci'); library(ggsci)}

# specify the output path
path_output_leaky <- '../output/analysis_1119/vector_control/output_files/leaky/'
path_output_all_or_none <- '../output/analysis_1119/vector_control/output_files/all_or_none/'

# list all of the output files 
files_output_leaky <- list.files(path = path_output_leaky, full.names = T, pattern = 'efficacy')
files_output_all_or_none <- list.files(path = path_output_all_or_none, full.names = T, pattern = 'efficacy')

# load the files 
output_leaky <- lapply(files_output_leaky, function(ff){read.csv(ff)})
output_efficacy_leaky <- as.data.frame((do.call('rbind', output_leaky)))

output_all_or_none <- lapply(files_output_all_or_none, function(ff){read.csv(ff)})
output_efficacy_all_or_none <- as.data.frame((do.call('rbind', output_all_or_none)))

# specify parameter ranges 
EIR_equil <- c(0.1 / 365, 1 / 365, 10 / 365, 100 / 365)
is_LLIN_distributed <- c(0,1)
is_IRS_administered <- c(0,1)

PSI_indoors <- c(0.9, 0.9, 0.9, 0.1)
PSI_bed <- c(0.9 * 0.5, 0.9 * 0.75, 0.9 * 0.25, 0.1 * 0.5)

# construct parameter grid 
n_rep <- 200
param_sweep <- expand.grid(EIR_equil = EIR_equil,
                           is_LLIN_distributed = is_LLIN_distributed,
                           is_IRS_administered = is_IRS_administered)
param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each = length(PSI_indoors)),]
param_sweep$PSI_indoors <- rep(PSI_indoors, nrow(param_sweep) / length(PSI_indoors))
param_sweep$PSI_bed <- rep(PSI_bed, nrow(param_sweep) / length(PSI_indoors))

param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each = n_rep), ]

# add in the relevant mosquito parameters 
output_efficacy_all_or_none$PSI_indoors <- param_sweep[output_efficacy_all_or_none$param_id + 1, 'PSI_indoors']
output_efficacy_all_or_none$PSI_bed <- param_sweep[output_efficacy_all_or_none$param_id + 1, 'PSI_bed']

output_efficacy_leaky$PSI_indoors <- param_sweep[output_efficacy_leaky$param_id + 1, 'PSI_indoors']
output_efficacy_leaky$PSI_bed <- param_sweep[output_efficacy_leaky$param_id + 1, 'PSI_bed']

# subset to the default scenario 
output_efficacy_all_default <- subset(output_efficacy_all_or_none, PSI_indoors == 0.9 & PSI_bed == 0.45)
output_efficacy_leaky_default <- subset(output_efficacy_leaky, PSI_indoors == 0.9 & PSI_bed == 0.45)

output_efficacy_all_endophagic_bed <- subset(output_efficacy_all_or_none, PSI_indoors == 0.9 & PSI_bed == (0.9 * 0.75))
output_efficacy_leaky_endophagic_bed <- subset(output_efficacy_leaky, PSI_indoors == 0.9 & PSI_bed == (0.9 * 0.75))

output_efficacy_all_endophagic_indoors <- subset(output_efficacy_all_or_none, PSI_indoors == 0.9 & PSI_bed == (0.9 * 0.25))
output_efficacy_leaky_endophagic_indoors <- subset(output_efficacy_leaky, PSI_indoors == 0.9 & PSI_bed == (0.9 * 0.25))

output_efficacy_all_exophagic <- subset(output_efficacy_all_or_none, PSI_indoors == 0.1 & PSI_bed == (0.1 * 0.5))
output_efficacy_leaky_exophagic <- subset(output_efficacy_leaky, PSI_indoors == 0.1 & PSI_bed == (0.1 * 0.5))

# specify vector control scenarios and range of eir values considered 
scenarios_vector_control <- data.frame(is_LLIN_distributed = c(0,1,0,1),
                                       is_IRS_administered = c(0,0,1,1))

eir_equil <- sort(unique(output_efficacy_leaky$eir_equil))

# calculate efficacy
eff_leaky_default <- lapply(1:nrow(scenarios_vector_control), function(vc){sapply(eir_equil, function(eir){quantile(output_efficacy_leaky_default$eff_cph_recurrent_LM[output_efficacy_leaky_default$eir_equil == eir & 
                                                                                                                                                                          output_efficacy_leaky_default$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] & 
                                                                                                                                                                          output_efficacy_leaky_default$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                                                             probs = c(0.25, 0.50, 0.75))})})

eff_all_default <- lapply(1:nrow(scenarios_vector_control), function(vc){sapply(eir_equil, function(eir){quantile(output_efficacy_all_default$eff_cph_recurrent_LM[output_efficacy_all_default$eir_equil == eir & 
                                                                                                                                                                                      output_efficacy_all_default$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] & 
                                                                                                                                                                                      output_efficacy_all_default$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                                                                   probs = c(0.25, 0.50, 0.75))})})

eff_leaky_endophagic_bed <- lapply(1:nrow(scenarios_vector_control), function(vc){sapply(eir_equil, function(eir){quantile(output_efficacy_leaky_endophagic_bed$eff_cph_recurrent_LM[output_efficacy_leaky_endophagic_bed$eir_equil == eir & 
                                                                                                                                                                         output_efficacy_leaky_endophagic_bed$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] & 
                                                                                                                                                                         output_efficacy_leaky_endophagic_bed$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                                                    probs = c(0.25, 0.50, 0.75))})})

eff_all_endophagic_bed <- lapply(1:nrow(scenarios_vector_control), function(vc){sapply(eir_equil, function(eir){quantile(output_efficacy_all_endophagic_bed$eff_cph_recurrent_LM[output_efficacy_all_endophagic_bed$eir_equil == eir & 
                                                                                                                                                                     output_efficacy_all_endophagic_bed$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] & 
                                                                                                                                                                     output_efficacy_all_endophagic_bed$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                                                  probs = c(0.25, 0.50, 0.75))})})

eff_leaky_endophagic_indoors <- lapply(1:nrow(scenarios_vector_control), function(vc){sapply(eir_equil, function(eir){quantile(output_efficacy_leaky_endophagic_indoors$eff_cph_recurrent_LM[output_efficacy_leaky_endophagic_indoors$eir_equil == eir & 
                                                                                                                                                                                       output_efficacy_leaky_endophagic_indoors$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] & 
                                                                                                                                                                                       output_efficacy_leaky_endophagic_indoors$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                                                           probs = c(0.25, 0.50, 0.75))})})

eff_all_endophagic_indoors <- lapply(1:nrow(scenarios_vector_control), function(vc){sapply(eir_equil, function(eir){quantile(output_efficacy_all_endophagic_indoors$eff_cph_recurrent_LM[output_efficacy_all_endophagic_indoors$eir_equil == eir & 
                                                                                                                                                                                   output_efficacy_all_endophagic_indoors$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] & 
                                                                                                                                                                                   output_efficacy_all_endophagic_indoors$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                                                         probs = c(0.25, 0.50, 0.75))})})

eff_leaky_exophagic <- lapply(1:nrow(scenarios_vector_control), function(vc){sapply(eir_equil, function(eir){quantile(output_efficacy_leaky_exophagic$eff_cph_recurrent_LM[output_efficacy_leaky_exophagic$eir_equil == eir & 
                                                                                                                                                                                               output_efficacy_leaky_exophagic$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] & 
                                                                                                                                                                                               output_efficacy_leaky_exophagic$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                                                               probs = c(0.25, 0.50, 0.75))})})

eff_all_exophagic <- lapply(1:nrow(scenarios_vector_control), function(vc){sapply(eir_equil, function(eir){quantile(output_efficacy_all_exophagic$eff_cph_recurrent_LM[output_efficacy_all_exophagic$eir_equil == eir & 
                                                                                                                                                                                           output_efficacy_all_exophagic$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] & 
                                                                                                                                                                                           output_efficacy_all_exophagic$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                                                             probs = c(0.25, 0.50, 0.75))})})
# generate plots 
palette <- head(pal_lancet()(9), n = 4)
offset <- seq(from = -0.4, to = 0.4, length.out = nrow(scenarios_vector_control) * 2)
cex = 1.5

jpeg(filename = '../output/figs_manuscript/fig_S2.jpg', width = 8, height = 11, units = 'in', res = 500)
par(mar = c(3.3,3.6,1.2,1.1))
layout(mat = matrix(1:4, nrow = 4, ncol = 1))
plot(NA, NA, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0, 1.0), axes = F, 
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
abline(h = 0.75, lwd = 1, lty = 2)
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(vv in 1:nrow(scenarios_vector_control))
  {
    segments(x0 = ee + offset[1 + 2 * (vv - 1)], y0 = eff_all_default[[vv]][1,ee], y1 = eff_all_default[[vv]][3,ee],
             lwd = cex, col = palette[vv])
    points(ee + offset[1 + 2 * (vv - 1)], eff_all_default[[vv]][2,ee], pch = 16, cex = cex, col = palette[vv])
    
    segments(x0 = ee + offset[2 + 2 * (vv - 1)], y0 = eff_leaky_default[[vv]][1,ee], y1 = eff_leaky_default[[vv]][3,ee],
             lwd = cex, col = palette[vv])
    points(ee + offset[2 + 2 * (vv - 1)], eff_leaky_default[[vv]][2,ee], pch = 17, cex = cex, col = palette[vv])
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 2, line = 2.3, 'Efficacy')
mtext(side = 3, line = 0, at = 0.525, font = 2, 'A')
mtext(side = 3, line = 0, at = mean(c(0.5, length(EIR_equil) + 0.5)),
      'Endophagic, No Preference in Biting Time')
legend('bottomleft', pch = c(1,2, NA, rep(15,4)), lwd = c(NA, NA, 1, NA, NA, NA, NA),
       lty = c(NA, NA, 2, NA, NA, NA, NA), pt.cex = cex, col = c('#222222', '#222222', '#222222', palette),
       legend = c('All-or-None', 'Leaky', 'Hypnozoite Clearance Prob', 'No Vector Control', 'LLINs only', 'IRS only', 'LLINs and IRS'),
       bty = 'n')

plot(NA, NA, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0, 1.0), axes = F, 
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
abline(h = 0.75, lwd = 1, lty = 2)
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(vv in 1:nrow(scenarios_vector_control))
  {
    segments(x0 = ee + offset[1 + 2 * (vv - 1)], y0 = eff_all_endophagic_bed[[vv]][1,ee], y1 = eff_all_endophagic_bed[[vv]][3,ee],
             lwd = cex, col = palette[vv])
    points(ee + offset[1 + 2 * (vv - 1)], eff_all_endophagic_bed[[vv]][2,ee], pch = 16, cex = cex, col = palette[vv])
    
    segments(x0 = ee + offset[2 + 2 * (vv - 1)], y0 = eff_leaky_endophagic_bed[[vv]][1,ee], y1 = eff_leaky_endophagic_bed[[vv]][3,ee],
             lwd = cex, col = palette[vv])
    points(ee + offset[2 + 2 * (vv - 1)], eff_leaky_endophagic_bed[[vv]][2,ee], pch = 17, cex = cex, col = palette[vv])
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 2, line = 2.3, 'Efficacy')
mtext(side = 3, line = 0, at = 0.525, font = 2, 'B')
mtext(side = 3, line = 0, at = mean(c(0.5, length(EIR_equil) + 0.5)),
      'Endophagic, Nighttime Biting')

plot(NA, NA, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0, 1.0), axes = F, 
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
abline(h = 0.75, lwd = 1, lty = 2)
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(vv in 1:nrow(scenarios_vector_control))
  {
    segments(x0 = ee + offset[1 + 2 * (vv - 1)], y0 = eff_all_endophagic_indoors[[vv]][1,ee], y1 = eff_all_endophagic_indoors[[vv]][3,ee],
             lwd = cex, col = palette[vv])
    points(ee + offset[1 + 2 * (vv - 1)], eff_all_endophagic_indoors[[vv]][2,ee], pch = 16, cex = cex, col = palette[vv])
    
    segments(x0 = ee + offset[2 + 2 * (vv - 1)], y0 = eff_leaky_endophagic_indoors[[vv]][1,ee], y1 = eff_leaky_endophagic_indoors[[vv]][3,ee],
             lwd = cex, col = palette[vv])
    points(ee + offset[2 + 2 * (vv - 1)], eff_leaky_endophagic_indoors[[vv]][2,ee], pch = 17, cex = cex, col = palette[vv])
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 2, line = 2.3, 'Efficacy')
mtext(side = 3, line = 0, at = 0.525, font = 2, 'C')
mtext(side = 3, line = 0, at = mean(c(0.5, length(EIR_equil) + 0.5)),
      'Endophagic, Daytime Biting')

plot(NA, NA, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0, 1.0), axes = F, 
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
abline(h = 0.75, lwd = 1, lty = 2)
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(vv in 1:nrow(scenarios_vector_control))
  {
    segments(x0 = ee + offset[1 + 2 * (vv - 1)], y0 = eff_all_exophagic[[vv]][1,ee], y1 = eff_all_exophagic[[vv]][3,ee],
             lwd = cex, col = palette[vv])
    points(ee + offset[1 + 2 * (vv - 1)], eff_all_exophagic[[vv]][2,ee], pch = 16, cex = cex, col = palette[vv])
    
    segments(x0 = ee + offset[2 + 2 * (vv - 1)], y0 = eff_leaky_exophagic[[vv]][1,ee], y1 = eff_leaky_exophagic[[vv]][3,ee],
             lwd = cex, col = palette[vv])
    points(ee + offset[2 + 2 * (vv - 1)], eff_leaky_exophagic[[vv]][2,ee], pch = 17, cex = cex, col = palette[vv])
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Efficacy')
mtext(side = 3, line = 0, at = 0.525, font = 2, 'D')
mtext(side = 3, line = 0, at = mean(c(0.5, length(EIR_equil) + 0.5)),
      'Exophagic, No Preference in Biting Time')
dev.off()



