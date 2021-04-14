# install necessary packages 
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(shape)){install.packages('shape'); library(shape)}

# specify path to output files 
path_output_all_or_none <- '../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/efficacy/'

# list all of the output files 
files_output_all_or_none <- list.files(path = path_output_all_or_none, full.names = T, pattern = 'efficacy')

# load the files 
output_all_or_none <- lapply(files_output_all_or_none, function(ff){read.csv(ff)})
output_efficacy_all_or_none <- as.data.frame((do.call('rbind', output_all_or_none)))

# get the unique eir and sig_het values 
eir_equil <- sort(unique(output_efficacy_all_or_none$eir_equil))
sig_het <- sort(unique(output_efficacy_all_or_none$sig_het))

# get the efficacies 
eff_all_or_none <- lapply(eir_equil, function(eir){sapply(sig_het, function(sig){quantile(output_efficacy_all_or_none$eff_cph_recurrent_LM[output_efficacy_all_or_none$sig_het == sig & output_efficacy_all_or_none$eir_equil == eir],
                                                                                          probs = c(0.25, 0.50, 0.75), na.rm = T)})})

# get the biting propensities to show the distribution
indiv_0 <- read.csv('../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/indiv/indiv_1600.csv.bz2')
indiv_1 <- read.csv('../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/indiv/indiv_1800.csv.bz2')
indiv_2 <- read.csv('../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/indiv/indiv_2000.csv.bz2')
indiv_3<- read.csv('../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/indiv/indiv_2200.csv.bz2')

# generate plot 
palette <- brewer.pal(n = 9, name = 'YlOrRd')
palette <- palette[c(3,5,7,9)]
offset <- seq(from = -0.35, to = 0.35, length.out = length(sig_het))

jpeg(filename = '../../output/figs/fig_2.jpg', width = 8, height = 5, units = 'in', res = 500)
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
    segments(x0 = ee + offset[1 + (ss-1)], y0 = eff_all_or_none[[ee]][1,ss], y1 = eff_all_or_none[[ee]][3,ss], col = palette[ss], lwd = 2)
    points(ee + offset[1 + (ss - 1)], eff_all_or_none[[ee]][2,ss], pch = 16, cex = 2, col = palette[ss])
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Efficacy')
mtext(side = 3, line = 0, at = 0.5, 'A', adj = 0, font = 2)
legend(x = 0.5, y = 0.375,
       lwd = c(1,rep(NA, length(palette))), col = c('#222222', palette),
       lty = c(2,rep(NA, length(palette))),
       pch = c(NA, rep(15, length(palette))),
       legend = c('Clearance Probability', '', '', '', ''),
       bty = 'n', pt.cex = 1.8)
Arrows(x0 = 0.62, x1 = 0.62, y0 = 0.29, y1 = 0.13, arr.type="triangle",
       arr.length = 0.15, arr.width = 0.1)
text(x = 0.675, y = (0.29 + 0.13)/2, 'Greater Biting Heterogeneity',
     pos = 4, offset = 0)

par(mar = c(3.6, 3.6, 1.2, 0.8))
hist(indiv_0$Zeta_Het, xlim = c(0,20), breaks = 0:100, axes = F, col = palette[1], border = 'white', freq = F,
     xlab = '', ylab = '', main = '', xaxs = 'i', yaxs = 'i', right = F)
box()
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Biting Exposure')
mtext(side = 2, line = 2.3, 'Density')
mtext(side = 3, line = 0, at = 0, 'B', adj = 0, font = 2)

hist(indiv_1$Zeta_Het, xlim = c(0,20), breaks = 0:100, axes = F, col = palette[2], border = 'white', freq = F,
     xlab = '', ylab = '', main = '', xaxs = 'i', yaxs = 'i', right = F)
box()
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Biting Exposure')
mtext(side = 3, line = 0, at = 0, 'C', adj = 0, font = 2)

hist(indiv_2$Zeta_Het, xlim = c(0,20), breaks = 0:100, axes = F, col = palette[3], border = 'white', freq = F, 
     xlab = '', ylab = '', main = '', xaxs = 'i', yaxs = 'i', right = F)
box()
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Biting Exposure')
mtext(side = 3, line = 0, at = 0, 'D', adj = 0, font = 2)

hist(indiv_3$Zeta_Het, xlim = c(0,20), breaks = 0:100, axes = F, col = palette[4], border = 'white', freq = F,
     xlab = '', ylab = '', main = '', xaxs = 'i', yaxs = 'i', right = F)
box()
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Biting Exposure')
mtext(side = 3, line = 0, at = 0, 'E', adj = 0, font = 2)

dev.off()

