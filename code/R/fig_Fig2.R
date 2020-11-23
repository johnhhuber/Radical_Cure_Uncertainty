# install necessary packages
if(!require(ggsci)){install.packages('ggsci'); library(ggsci)}
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# specify the output path
path_output_leaky <- '../../output/analysis/followup_vs_relapse/output_files/leaky/efficacy/'
path_output_all_or_none <- '../../output/analysis/followup_vs_relapse/output_files/all_or_none/efficacy/'

# list all of the output files 
files_output_leaky <- list.files(path = path_output_leaky, full.names = T, pattern = 'efficacy')
files_output_all_or_none <- list.files(path = path_output_all_or_none, full.names = T, pattern = 'efficacy')

# load the files 
output_leaky <- lapply(files_output_leaky, function(ff){read.csv(ff)})
output_efficacy_leaky <- as.data.frame((do.call('rbind', output_leaky)))

output_all_or_none <- lapply(files_output_all_or_none, function(ff){read.csv(ff)})
output_efficacy_all_or_none <- as.data.frame((do.call('rbind', output_all_or_none)))

# get the unique eirs, relapse rates and trial durations 
eir_equil <- sort(unique(output_efficacy_all_or_none$eir_equil))
ff <- rev(sort(unique(output_efficacy_all_or_none$ff)))
trial_duration <- sort(unique(output_efficacy_all_or_none$trial_duration))

# get the efficacies for the leaky and all or none campaigns 
eff_leaky <- lapply(eir_equil, function(eir){lapply(trial_duration, function(dur){sapply(ff, function(rel){quantile(output_efficacy_leaky$eff_cph_recurrent_LM[output_efficacy_leaky$eir_equil == eir & output_efficacy_leaky$trial_duration == dur &
                                                                                                                                                                 abs(output_efficacy_leaky$ff - rel) < 1e-6], probs = c(0.25, 0.50, 0.75))})})})

eff_all_or_none <- lapply(eir_equil, function(eir){lapply(trial_duration, function(dur){sapply(ff, function(rel){quantile(output_efficacy_all_or_none$eff_cph_recurrent_LM[output_efficacy_all_or_none$eir_equil == eir & output_efficacy_all_or_none$trial_duration == dur &
                                                                                                                                                                 abs(output_efficacy_all_or_none$ff - rel) < 1e-6], probs = c(0.25, 0.50, 0.75))})})})

# specify mosquito-to-human transmission prob
mosquito_to_human_transmission_prob = 0.5

# generate plot 
palette <- brewer.pal(n = 9, name = 'YlOrRd')
palette <- palette[c(9,7,5,3)]
offset <- seq(from = -0.475, to = 0.475, length.out = 2 * length(ff) * length(trial_duration))
rect_width <- c(-0.5, mean(offset[8:9]), mean(offset[16:17]), mean(offset[24:25]), 0.5)

jpeg(filename = '../../output/figs/fig_2.jpg', width = 7, height = 5, units = 'in', res = 500)
layout(mat = matrix(c(rep(1,6),rep(c(2,3),2)), nrow = 5, byrow = T))
par(mar = c(3.3,3.6,1.2,0.8))
plot(NA, NA, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0,1.1), axes = F,
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
for(ee in 1:length(eir_equil))
{
  for(ss in 1:(length(rect_width) - 1))
  {
    rect(xleft = ee + rect_width[ss], xright = ee + rect_width[ss+1],
         ybottom = 0, ytop = 1, col = ifelse(ss %% 2 == 0, col2alpha('grey', alpha = 0.325), 'white'), border = NA)
    text(x = ee + mean(rect_width[ss:(ss+1)]), y = 1.025, trial_duration[ss])
    text(x = ee, y = 1.075, 'Duration of Follow-up (d)', cex = 0.8)
    segments(x0 = ee + rect_width[ss+1], y0 = 1, y1 = 1.05)
  }
}
abline(h = c(1,1.05))
abline(h = 0.75, lwd = 1, lty = 2)
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1))
for(ee in 1:length(eir_equil))
{
  for(dd in 1:length(trial_duration))
  {
    for(rr in 1:length(ff))
    {
      segments(x0 = ee + offset[1 + 8 * (dd-1) + 2 * (rr-1)], y0 = eff_all_or_none[[ee]][[dd]][1,rr],
               y1 = eff_all_or_none[[ee]][[dd]][3,rr], col = palette[rr], lwd = 1)
      points(ee + offset[1 + 8 * (dd-1) + 2 * (rr-1)], eff_all_or_none[[ee]][[dd]][2,rr], pch = 16, cex = 0.9, col = palette[rr])
      
      segments(x0 = ee + offset[2 + 8 * (dd-1) + 2 * (rr-1)], y0 = eff_leaky[[ee]][[dd]][1,rr],
               y1 = eff_leaky[[ee]][[dd]][3,rr], col = palette[rr], lwd = 1)
      points(ee + offset[2 + 8 * (dd-1) + 2 * (rr-1)], eff_leaky[[ee]][[dd]][2,rr], pch = 17, cex = 0.9, col = palette[rr])
    }
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 3, line = 0, at = 0.5, 'A', font = 2, adj = 0, cex = 0.8)
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Efficacy')

par(xpd = T)
# y = 1.1825
legend(x = mean(c(0.5, length(eir_equil)+0.5)), y = 1.2, pch = c(1,2,NA),
       lwd = c(NA, NA, 1), col = '#222222',
       lty = c(NA, NA, 2),
       legend = c('All-or-None', 'Leaky', 'Hypnozoite Clearance Prob'),
       bty = 'n', pt.cex = 1.25, ncol = 3, xjust = 0.5)

par(xpd = F)
par(mar = c(3.6,3.6,1.8,0.8))
palette <- (pal_material('orange', n = 10)(10)[c(10, 8, 6, 4)])
palette <- brewer.pal(n = 9, name = 'YlOrRd')[c(9,7,5,3)]
plot(NA, NA, xlim = c(0,max(trial_duration)), ylim = c(0,1.0),
     xaxs = 'i', yaxs = 'i', axes = F, xlab = '', ylab = '')
for(ii in 1:length(ff))
{
  hazard <- rep(ff[ii], max(trial_duration))
  lines(1:max(trial_duration), (1 - exp(-cumsum(hazard))), col = palette[ii], lwd = 2)
}
box()
axis(side = 1, at = c(0,trial_duration), labels = c(0,trial_duration))
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 1, line = 2.3, 'Duration of Follow-up (d)', cex = 0.8)
mtext(side = 2, line = 2.3, 'Relapsed (%)', cex = 0.8)
mtext(side = 3, line = 0, at = 0, 'B', font = 2, adj = 0, cex = 0.8)

par(xpd = T)
# y = 1.275
legend(x = 0.5 * max(trial_duration), y = 1.21, lwd = 2, col = palette, legend = (round(1/ff)), bty = 'n', ncol = length(palette),
       title = 'Mean Time to Relapse (d)', xjust = 0.5, cex = 0.8)

par(xpd = F)
palette <- rev(pal_material('grey', n = 10)(10)[c(10, 7, 5, 4)])
plot(NA, NA, xlim = c(0, max(trial_duration)), ylim = c(0,1.0),
     xaxs = 'i', yaxs = 'i', axes = F, xlab = '', ylab = '')
for(ii in 1:length(eir_equil))
{  
  hazard <- rep(eir_equil[ii] / 365, max(trial_duration)) * mosquito_to_human_transmission_prob
  lines(1:max(trial_duration), 1 - exp(-cumsum(hazard)), col = palette[ii], lwd = 2)
}
box()
axis(side = 1, at = c(0,trial_duration), labels = c(0,trial_duration))
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 1, line = 2.3, 'Duration of Follow-up (d)', cex = 0.8)
mtext(side = 2, line = 2.3, 'Reinfected (%)', cex = 0.8)
mtext(side = 3, line = 0, at = 0, 'C', font = 2, adj = 0, cex = 0.8)

par(xpd = T)
legend(x = 0.5 * max(trial_duration), y = 1.21, lwd = 2, col = palette, legend = eir_equil, bty = 'n', ncol = length(palette),
       title = 'EIR', xjust = 0.5, cex = 0.8)
par(xpd = F)

dev.off()
