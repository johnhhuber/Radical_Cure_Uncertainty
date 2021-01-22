# install necessary packages
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# specify the output path
path_output_all_or_none <- '../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/genotyping/'

# list all of the output files 
files_output_all_or_none <- list.files(path = path_output_all_or_none, full.names = T, pattern = 'genotyping')

# load the files 
output_all_or_none <- lapply(files_output_all_or_none, function(ff){read.csv(ff)})
output_efficacy_all_or_none <- as.data.frame((do.call('rbind', output_all_or_none)))

# subset to only include homogeneous biting scenario 
output_efficacy_all_or_none <- subset(output_efficacy_all_or_none, sig_het == 0)

# get the unique eir
eir_equil <- sort(unique(output_efficacy_all_or_none$eir_equil))
sens <- seq(from = 0.25, to = 1, by = 0.25)
spec <- seq(from = 0.25, to = 1, by = 0.25)

# calculate the efficacy under different genotyping scenarios
cols_genotyping <- grepl('SENS', colnames(output_efficacy_all_or_none))

eff_all_or_none <- lapply(eir_equil, function(eir){apply(output_efficacy_all_or_none[output_efficacy_all_or_none$eir_equil == eir, cols_genotyping], 2, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)})

# plot the results 
palette <- brewer.pal(n = 9, name = 'YlOrRd')
palette <- palette[c(3,5,7,9)]
offset <- seq(from = -0.475, to = 0.475, length.out = 2 * length(sens) * length(spec))
rect_width <- c(-0.5, mean(offset[8:9]), mean(offset[16:17]), mean(offset[24:25]), 0.5)

genotyping_method <- 'time_LM_'

jpeg(filename = '../../output/figs/fig_5.jpg', width = 8, height = 4, units = 'in', res = 500)
par(mar = c(4.3, 3.6, 0.8, 0.8))
plot(NA, NA, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0,1.1), axes = F,
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
for(ee in 1:length(eir_equil))
{
  for(ss in 1:(length(rect_width) - 1))
  {
    rect(xleft = ee + rect_width[ss], xright = ee + rect_width[ss+1],
         ybottom = 0, ytop = 1.0, col = ifelse(ss %% 2 == 0, col2alpha('grey', alpha = 0.325), 'white'), border = NA)
    text(x = ee + mean(rect_width[ss:(ss+1)]), y = 1.025, spec[ss] * 100, cex = 0.8)
    text(x = ee, y = 1.075, 'Specificity (%)', cex = 0.8)
    segments(x0 = ee + rect_width[ss+1], y0 = 1, y1 = 1.05)
  }
}
abline(h = c(1,1.05))
abline(h = 0.75, lty = 2)
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1))
for(ee in 1:length(eir_equil))
{
  for(sp in 1:length(spec))
  {
    for(se in 1:length(sens))
    {
      mean_offset <- 0.5 * (offset[1 + 8 * (sp-1) + 2 * (se-1)] + offset[2 + 8 * (sp-1) + 2 * (se-1)])
      segments(x0 = ee + mean_offset, y0 = eff_all_or_none[[ee]][1, paste(genotyping_method, 'SENS_', sens[se] * 100, '_SPEC_', spec[sp] * 100, sep = '')],
               y1 = eff_all_or_none[[ee]][3, paste(genotyping_method, 'SENS_', sens[se] * 100, '_SPEC_', spec[sp] * 100, sep = '')], col = palette[se], lwd = 1.25)
      points(ee + mean_offset, eff_all_or_none[[ee]][2, paste(genotyping_method, 'SENS_', sens[se] * 100, '_SPEC_', spec[sp] * 100, sep = '')], pch = 16, cex = 1.25, col = palette[se])
    }
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Efficacy')

par(xpd = T)
legend(x = mean(c(0.5, length(eir_equil)+0.5)), y = 1.19,
       lwd = 1, col = '#222222',
       lty = 2,
       legend = 'Clearance Probability',
       bty = 'n', pt.cex = 1.1, cex = 0.8, ncol = 1, xjust = 0.5)

par(xpd = F)

par(xpd = T)
legend(x = 0.5, y = -0.15, pch = 15,
      col = palette, legend = c('25', '50', '75', '100'),
      bty = 'n', pt.cex = 1.25, cex = 0.9, ncol = length(palette), xjust = 0,
      title = 'Sensitivity (%): ', title.adj = 0.5)
par(xpd = F)

dev.off()
