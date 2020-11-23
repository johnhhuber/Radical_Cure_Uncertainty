# set working directory
setwd('~/Dropbox/Radical_Cure_MOA/code/')

# clear existing workspace
rm(list = ls())

# install necessary packages 
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}

# specify path to output files 
path_output_leaky <- '../output/analysis_1119/eir_vs_heterogeneity/output_files/leaky/'
path_output_all_or_none <- '../output/analysis_1119/eir_vs_heterogeneity/output_files/all_or_none/'

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
sig_het <- sort(unique(output_efficacy_all_or_none$sig_het))

# get the efficacies 
eff_leaky <- lapply(eir_equil, function(eir){sapply(sig_het, function(sig){quantile(output_efficacy_leaky$eff_cph_recurrent_LM[output_efficacy_leaky$sig_het == sig & output_efficacy_leaky$eir_equil == eir],
                                                                                    probs = c(0.25, 0.50, 0.75))})})
eff_all_or_none <- lapply(eir_equil, function(eir){sapply(sig_het, function(sig){quantile(output_efficacy_all_or_none$eff_cph_recurrent_LM[output_efficacy_all_or_none$sig_het == sig & output_efficacy_all_or_none$eir_equil == eir],
                                                                                          probs = c(0.25, 0.50, 0.75))})})

# get the biting propensities to show the distribution
indiv_0 <- read.csv('../output/analysis_1119/eir_vs_heterogeneity/output_files/leaky/indiv/indiv_0.csv.bz2')
indiv_1 <- read.csv('../output/analysis_1119/eir_vs_heterogeneity/output_files/leaky/indiv/indiv_200.csv.bz2')
indiv_2 <- read.csv('../output/analysis_1119/eir_vs_heterogeneity/output_files/leaky/indiv/indiv_400.csv.bz2')
indiv_3 <- read.csv('../output/analysis_1119/eir_vs_heterogeneity/output_files/leaky/indiv/indiv_600.csv.bz2')

# generate plot 
#palette <- pal_material(palette = 'deep-orange', n = 8)(8)[-1]
#palette <- palette[seq(from = 1, to = length(palette), length.out = length(sig_het))]
palette <- brewer.pal(n = 9, name = 'YlOrRd')
palette <- palette[c(3,5,7,9)]
offset <- seq(from = -0.4, to = 0.4, length.out = length(sig_het) * 2)

jpeg(filename = '../output/figs_manuscript/fig_1.jpg', width = 8, height = 5, units = 'in', res = 500)
par(mar = c(3.6, 3.6, 1.2, 0.8))
layout(mat = matrix(c(rep(1, 8), 2,3,4,5), nrow = 3, byrow = T))
plot(NA, NA, type = 'n', axes = F, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0, 1),
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(h = 0.75, lwd = 1, lty = 2)
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(ss in 1:length(sig_het))
  {
    segments(x0 = ee + offset[1 + 2 * (ss-1)], y0 = eff_all_or_none[[ee]][1,ss], y1 = eff_all_or_none[[ee]][3,ss], col = palette[ss], lwd = 2)
    points(ee + offset[1 + 2 * (ss - 1)], eff_all_or_none[[ee]][2,ss], pch = 16, cex = 2, col = palette[ss])
    segments(x0 = ee + offset[2 + 2 * (ss-1)], y0 = eff_leaky[[ee]][1,ss], y1 = eff_leaky[[ee]][3,ss], col = palette[ss], lwd = 2)
    points(ee + offset[2 + 2 * (ss - 1)], eff_leaky[[ee]][2,ss], pch = 17, cex = 2, col = palette[ss])
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Efficacy')
mtext(side = 3, line = 0, at = 0.5, 'A', adj = 0, font = 2)
legend('bottomleft', pch = c(1,2,NA),
       lwd = c(NA, NA, 1), col = '#222222',
       lty = c(NA, NA, 2),
       legend = c('All-or-None', 'Leaky', 'Hypnozoite Clearance Prob'),
       bty = 'n', pt.cex = 1.8)

par(mar = c(3.6, 3.6, 1.2, 0.8))
hist(indiv_0$Zeta_Het, xlim = c(0,20), breaks = 0:100, axes = F, col = palette[1], border = 'white', freq = F,
     xlab = '', ylab = '', main = '', xaxs = 'i', yaxs = 'i')
box()
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Biting Propensity')
mtext(side = 2, line = 2.3, 'Density')
mtext(side = 3, line = 0, at = 0, 'B', adj = 0, font = 2)

hist(indiv_1$Zeta_Het, xlim = c(0,20), breaks = 0:100, axes = F, col = palette[2], border = 'white', freq = F,
     xlab = '', ylab = '', main = '', xaxs = 'i', yaxs = 'i')
box()
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Biting Propensity')
mtext(side = 3, line = 0, at = 0, 'C', adj = 0, font = 2)

hist(indiv_2$Zeta_Het, xlim = c(0,20), breaks = 0:100, axes = F, col = palette[3], border = 'white', freq = F, 
     xlab = '', ylab = '', main = '', xaxs = 'i', yaxs = 'i')
box()
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Biting Propensity')
mtext(side = 3, line = 0, at = 0, 'D', adj = 0, font = 2)

hist(indiv_3$Zeta_Het, xlim = c(0,20), breaks = 0:100, axes = F, col = palette[4], border = 'white', freq = F,
     xlab = '', ylab = '', main = '', xaxs = 'i', yaxs = 'i')
box()
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Biting Propensity')
mtext(side = 3, line = 0, at = 0, 'E', adj = 0, font = 2)

dev.off()

