# set working directory
setwd('~/Dropbox/Radical_Cure_MOA/code/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(ggsci)){install.packages('ggsci'); library(ggsci)}

# specify path to output files 
path_output_leaky <- '../output/analysis_1119/radical_cure_therapeutic/output_files/leaky/'
path_output_all_or_none <- '../output/analysis_1119/radical_cure_therapeutic/output_files/all_or_none/'

# list all of the output files 
files_output_leaky <- list.files(path = path_output_leaky, full.names = T, pattern = 'efficacy')
files_output_all_or_none <- list.files(path = path_output_all_or_none, full.names = T, pattern = 'efficacy')

# load the files 
output_leaky <- lapply(files_output_leaky, function(ff){read.csv(ff)})
output_efficacy_leaky <- as.data.frame((do.call('rbind', output_leaky)))

output_all_or_none <- lapply(files_output_all_or_none, function(ff){read.csv(ff)})
output_efficacy_all_or_none <- as.data.frame((do.call('rbind', output_all_or_none)))

# get the unique eir and prophylaxis values 
eir_equil <- sort(unique(output_efficacy_all_or_none$eir_equil))
trial_PQ_proph <- sort(unique(output_efficacy_all_or_none$trial_PQ_proph))

# get the efficacies 
eff_leaky <- lapply(eir_equil, function(eir){sapply(trial_PQ_proph, function(proph){quantile(output_efficacy_leaky$eff_cph_recurrent_LM[output_efficacy_leaky$trial_PQ_proph == proph & output_efficacy_leaky$eir_equil == eir],
                                                                                    probs = c(0.25, 0.50, 0.75))})})
eff_all_or_none <- lapply(eir_equil, function(eir){sapply(trial_PQ_proph, function(proph){quantile(output_efficacy_all_or_none$eff_cph_recurrent_LM[output_efficacy_all_or_none$trial_PQ_proph == proph & output_efficacy_all_or_none$eir_equil == eir],
                                                                                          probs = c(0.25, 0.50, 0.75))})})

# generate plot 
palette <- c('#9A31CD', '#79CDCD')
offset <- seq(from = -0.25, to = 0.25, length.out = length(trial_PQ_proph) * 2)

jpeg(filename = '../output/figs_manuscript/fig_5.jpg', width = 8, height = 4, units = 'in', res = 500)
par(mar = c(3.6, 3.6, 0.8, 0.8))
plot(NA, NA, type = 'n', axes = F, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0, 1),
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(h = 0.75, lwd = 1, lty = 2)
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(pp in 1:length(trial_PQ_proph))
  {
    segments(x0 = ee + offset[1 + 2 * (pp-1)], y0 = eff_all_or_none[[ee]][1,pp], y1 = eff_all_or_none[[ee]][3,pp], col = palette[pp], lwd = 1.5)
    points(ee + offset[1 + 2 * (pp - 1)], eff_all_or_none[[ee]][2,pp], pch = 16, cex = 1.5, col = palette[pp])
    segments(x0 = ee + offset[2 + 2 * (pp-1)], y0 = eff_leaky[[ee]][1,pp], y1 = eff_leaky[[ee]][3,pp], col = palette[pp], lwd = 1.5)
    points(ee + offset[2 + 2 * (pp - 1)], eff_leaky[[ee]][2,pp], pch = 17, cex = 1.5, col = palette[pp])
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Efficacy')
legend('bottomleft', pch = c(1,2, 15, 15, NA),
       lwd = c(NA, NA, NA, NA, 1), col = c('#222222', '#222222', palette, '#222222'),
       lty = c(NA, NA, NA, NA, 2),
       legend = c('All-or-None', 'Leaky', 'Primaquine', 'Tafenoquine', 'Hyp Clearance Prob'),
       bty = 'n', pt.cex = 1.5, cex = 0.8)
dev.off()
