# install necessary packages 
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(shape)){install.packages('shape'); library(shape)}

# specify path to output files 
path_output_leaky <- '../../output/analysis/eir_vs_heterogeneity/output_files/leaky/efficacy/'
path_output_all_or_none <- '../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/efficacy/'

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
                                                                                    probs = c(0.25, 0.50, 0.75), na.rm = T)})})
eff_all_or_none <- lapply(eir_equil, function(eir){sapply(sig_het, function(sig){quantile(output_efficacy_all_or_none$eff_cph_recurrent_LM[output_efficacy_all_or_none$sig_het == sig & output_efficacy_all_or_none$eir_equil == eir],
                                                                                          probs = c(0.25, 0.50, 0.75), na.rm = T)})})

# generate plot 
palette <- brewer.pal(n = 9, name = 'YlOrRd')
palette <- palette[c(3,5,7,9)]
offset <- seq(from = -0.4, to = 0.4, length.out = 2 * length(sig_het))

jpeg(filename = '../../output/figs/fig_S6.jpg', width = 8, height = 5, units = 'in', res = 500)
par(mar = c(3.6, 3.6, 0.8, 0.8))
plot(NA, NA, type = 'n', axes = F, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0, 1),
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(h = 0.75, lwd = 1, lty = 2)
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(ss in 1:length(sig_het))
  {
    segments(x0 = ee + offset[1 + 2 * (ss-1)], y0 = eff_all_or_none[[ee]][1,ss], y1 = eff_all_or_none[[ee]][3,ss], col = palette[ss], lwd = 1.5)
    points(ee + offset[1 + 2 * (ss - 1)], eff_all_or_none[[ee]][2,ss], pch = 16, cex = 1.5, col = palette[ss])
    
    segments(x0 = ee + offset[2 + 2 * (ss-1)], y0 = eff_leaky[[ee]][1,ss], y1 = eff_leaky[[ee]][3,ss], col = palette[ss], lwd = 1.5)
    points(ee + offset[2 + 2 * (ss - 1)], eff_leaky[[ee]][2,ss], pch = 17, cex = 1.5, col = palette[ss])
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Efficacy')

legend(x = 0.5, y = 0.375,
       lwd = c(NA, NA, 1,rep(NA, length(palette))), col = c('#222222','#222222','#222222', palette),
       lty = c(NA, NA, 2,rep(NA, length(palette))),
       pch = c(1,2,NA, rep(15, length(palette))),
       legend = c('All-or-None','Leaky','Clearance Prob', '', '', '', ''),
       bty = 'n', pt.cex = 1.5, cex = 0.8)
Arrows(x0 = 0.675, x1 = 0.675, y0 = 0.23, y1 = 0.0975, arr.type="triangle",
       arr.length = 0.15, arr.width = 0.1)
text(x = 0.70, y = (0.23 + 0.0975)/2, 'Greater Heterogeneity',
     pos = 4, offset = 0, cex = 0.8)

dev.off()