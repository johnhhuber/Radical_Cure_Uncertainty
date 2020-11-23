# install necessary packages
if(!require(ggsci)){install.packages('ggsci'); library(ggsci)}

# specify paths in 
path_output_leaky <- '../../output/analysis/vector_control/output_files/leaky/'
path_output_all_or_none <- '../../output/analysis/vector_control/output_files/all_or_none/'

# list all of the output files
files_output_leaky <- list.files(path = path_output_leaky, full.names = T, pattern = 'recurrent')
files_output_all_or_none <- list.files(path = path_output_all_or_none, full.names = T, pattern = 'recurrent')

# load the files 
output_leaky <- lapply(files_output_leaky, function(ff){read.csv(ff)})
output_leaky <- as.data.frame((do.call('rbind', output_leaky)))

output_all_or_none <- lapply(files_output_all_or_none, function(ff){read.csv(ff)})
output_all_or_none <- as.data.frame((do.call('rbind', output_all_or_none)))

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
output_all_or_none$PSI_indoors <- param_sweep[output_all_or_none$param_id + 1, 'PSI_indoors']
output_all_or_none$PSI_bed <- param_sweep[output_all_or_none$param_id + 1, 'PSI_bed']
output_all_or_none$EIR_equil <- param_sweep[output_all_or_none$param_id + 1, 'EIR_equil'] * 365
output_all_or_none$is_LLIN_distributed <- param_sweep[output_all_or_none$param_id + 1, 'is_LLIN_distributed']
output_all_or_none$is_IRS_administered <- param_sweep[output_all_or_none$param_id + 1, 'is_IRS_administered']

output_leaky$PSI_indoors <- param_sweep[output_leaky$param_id + 1, 'PSI_indoors']
output_leaky$PSI_bed <- param_sweep[output_leaky$param_id + 1, 'PSI_bed']
output_leaky$EIR_equil <- param_sweep[output_leaky$param_id + 1, 'EIR_equil'] * 365
output_leaky$is_LLIN_distributed <- param_sweep[output_leaky$param_id + 1, 'is_LLIN_distributed']
output_leaky$is_IRS_administered <- param_sweep[output_leaky$param_id + 1, 'is_IRS_administered']

# specify vector control scenarios 
scenarios_vector_control <- data.frame(is_LLIN_distributed = c(0,1,0,1),
                                       is_IRS_administered = c(0,0,1,1))

# calculate the quantiles for the proportion of recurrent infections that are relapses
eir_equil <- sort(unique(output_leaky$EIR_equil))           

# subset to the various vector scenarios 
output_all_default <- subset(output_all_or_none, PSI_indoors == 0.9 & PSI_bed == 0.45)
output_leaky_default <- subset(output_leaky, PSI_indoors == 0.9 & PSI_bed == 0.45)

output_all_endophagic_bed <- subset(output_all_or_none, PSI_indoors == 0.9 & PSI_bed == (0.9 * 0.75))
output_leaky_endophagic_bed <- subset(output_leaky, PSI_indoors == 0.9 & PSI_bed == (0.9 * 0.75))

output_all_endophagic_indoors <- subset(output_all_or_none, PSI_indoors == 0.9 & PSI_bed == (0.9 * 0.25))
output_leaky_endophagic_indoors <- subset(output_leaky, PSI_indoors == 0.9 & PSI_bed == (0.9 * 0.25))

output_all_exophagic <- subset(output_all_or_none, PSI_indoors == 0.1 & PSI_bed == (0.1 * 0.5))
output_leaky_exophagic <- subset(output_leaky, PSI_indoors == 0.1 & PSI_bed == (0.1 * 0.5))

# calculate the proportion of recurrent infections that are true treatment failures 
prop_leaky_default <- lapply(1:nrow(scenarios_vector_control), 
                            function(vc){sapply(eir_equil, function(eir){quantile(output_leaky_default$prop_relapse[output_leaky_default$EIR_equil == eir & 
                                                                                                            output_leaky_default$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] &
                                                                                                            output_leaky_default$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                    probs = c(0.25,0.50,0.75), na.rm = T)})})
prop_all_default <- lapply(1:nrow(scenarios_vector_control), 
                             function(vc){sapply(eir_equil, function(eir){quantile(output_all_default$prop_relapse[output_all_default$EIR_equil == eir & 
                                                                                                                       output_all_default$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] &
                                                                                                                       output_all_default$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                   probs = c(0.25,0.50,0.75), na.rm = T)})})

prop_leaky_endophagic_bed <- lapply(1:nrow(scenarios_vector_control), 
                             function(vc){sapply(eir_equil, function(eir){quantile(output_leaky_endophagic_bed$prop_relapse[output_leaky_endophagic_bed$EIR_equil == eir & 
                                                                                                                       output_leaky_endophagic_bed$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] &
                                                                                                                       output_leaky_endophagic_bed$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                   probs = c(0.25,0.50,0.75), na.rm = T)})})
prop_all_endophagic_bed <- lapply(1:nrow(scenarios_vector_control), 
                           function(vc){sapply(eir_equil, function(eir){quantile(output_all_endophagic_bed$prop_relapse[output_all_endophagic_bed$EIR_equil == eir & 
                                                                                                                   output_all_endophagic_bed$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] &
                                                                                                                   output_all_endophagic_bed$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                 probs = c(0.25,0.50,0.75), na.rm = T)})})

prop_leaky_endophagic_indoors <- lapply(1:nrow(scenarios_vector_control), 
                                    function(vc){sapply(eir_equil, function(eir){quantile(output_leaky_endophagic_indoors$prop_relapse[output_leaky_endophagic_indoors$EIR_equil == eir & 
                                                                                                                                     output_leaky_endophagic_indoors$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] &
                                                                                                                                     output_leaky_endophagic_indoors$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                          probs = c(0.25,0.50,0.75), na.rm = T)})})
prop_all_endophagic_indoors <- lapply(1:nrow(scenarios_vector_control), 
                                  function(vc){sapply(eir_equil, function(eir){quantile(output_all_endophagic_indoors$prop_relapse[output_all_endophagic_indoors$EIR_equil == eir & 
                                                                                                                                 output_all_endophagic_indoors$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] &
                                                                                                                                 output_all_endophagic_indoors$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                        probs = c(0.25,0.50,0.75), na.rm = T)})})

prop_leaky_exophagic <- lapply(1:nrow(scenarios_vector_control), 
                                        function(vc){sapply(eir_equil, function(eir){quantile(output_leaky_exophagic$prop_relapse[output_leaky_exophagic$EIR_equil == eir & 
                                                                                                                                             output_leaky_exophagic$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] &
                                                                                                                                             output_leaky_exophagic$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                              probs = c(0.25,0.50,0.75), na.rm = T)})})
prop_all_exophagic <- lapply(1:nrow(scenarios_vector_control), 
                                      function(vc){sapply(eir_equil, function(eir){quantile(output_all_exophagic$prop_relapse[output_all_exophagic$EIR_equil == eir & 
                                                                                                                                         output_all_exophagic$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] &
                                                                                                                                         output_all_exophagic$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                            probs = c(0.25,0.50,0.75), na.rm = T)})})

# generate plot
rectwidth = 0.07
palette <- head(pal_lancet()(9), n = 4)
offset <- seq(from = -0.3, to = 0.3, length.out = nrow(scenarios_vector_control) * 2)

jpeg(filename = '../../output/figs/fig_S3.jpg', width = 8, height = 11, units = 'in', res = 500)
layout(mat = matrix(1:4, nrow = 4, ncol = 1))
par(mar = c(3.3,3.6,1.2,1.1))
plot(NA, NA, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(vv in 1:nrow(scenarios_vector_control))
  {
    rect(xleft = ee + offset[1 + 2 * (vv-1)] - (rectwidth/2), xright = ee + offset[1 + 2 * (vv-1)] + (rectwidth/2),
         ybottom = 0, ytop = prop_all_default[[vv]][2,ee], col = palette[vv])
    segments(x0 = ee + offset[1 + 2 * (vv-1)], y0 = prop_all_default[[vv]][1,ee], y1 = prop_all_default[[vv]][3,ee],
             col = '#222222', lwd = 1.5)
    rect(xleft = ee + offset[2 + 2 * (vv-1)] - (rectwidth/2), xright = ee + offset[2 + 2 * (vv-1)] + (rectwidth/2),
         ybottom = 0, ytop = prop_leaky_default[[vv]][2,ee], col = palette[vv], density = 45, border = '#222222')
    segments(x0 = ee + offset[2 + 2 * (vv-1)], y0 = prop_leaky_default[[vv]][1,ee], y1 = prop_leaky_default[[vv]][3,ee],
             col = '#222222', lwd = 1.5)
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.20), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 2, line = 2.3, 'Treatment Failures (%)')
mtext(side = 3, line = 0, at = 0.525, font = 2, 'A')
mtext(side = 3, line = 0, at = mean(c(0.5, length(EIR_equil) + 0.5)),
      'Endophagic, No Preference in Biting Time')
legend('topright', pch = c(15, 7, rep(16, 4)), col = c('#222222', '#222222', palette),
       legend = c('All-or-None', 'Leaky', 'No Vector Control', 'LLINs only', 'IRS only', 'LLINs and IRS'),
       bty = 'n', pt.cex = 2)

plot(NA, NA, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(vv in 1:nrow(scenarios_vector_control))
  {
    rect(xleft = ee + offset[1 + 2 * (vv-1)] - (rectwidth/2), xright = ee + offset[1 + 2 * (vv-1)] + (rectwidth/2),
         ybottom = 0, ytop = prop_all_endophagic_bed[[vv]][2,ee], col = palette[vv])
    segments(x0 = ee + offset[1 + 2 * (vv-1)], y0 = prop_all_endophagic_bed[[vv]][1,ee], y1 = prop_all_endophagic_bed[[vv]][3,ee],
             col = '#222222', lwd = 1.5)
    rect(xleft = ee + offset[2 + 2 * (vv-1)] - (rectwidth/2), xright = ee + offset[2 + 2 * (vv-1)] + (rectwidth/2),
         ybottom = 0, ytop = prop_leaky_endophagic_bed[[vv]][2,ee], col = palette[vv], density = 45, border = '#222222')
    segments(x0 = ee + offset[2 + 2 * (vv-1)], y0 = prop_leaky_endophagic_bed[[vv]][1,ee], y1 = prop_leaky_endophagic_bed[[vv]][3,ee],
             col = '#222222', lwd = 1.5)
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.20), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 3, line = 0, at = 0.525, font = 2, 'B')
mtext(side = 3, line = 0, at = mean(c(0.5, length(EIR_equil) + 0.5)),
      'Endophagic, Nighttime Biting')
mtext(side = 2, line = 2.3, 'Treatment Failures (%)')

plot(NA, NA, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(vv in 1:nrow(scenarios_vector_control))
  {
    rect(xleft = ee + offset[1 + 2 * (vv-1)] - (rectwidth/2), xright = ee + offset[1 + 2 * (vv-1)] + (rectwidth/2),
         ybottom = 0, ytop = prop_all_endophagic_indoors[[vv]][2,ee], col = palette[vv])
    segments(x0 = ee + offset[1 + 2 * (vv-1)], y0 = prop_all_endophagic_indoors[[vv]][1,ee], y1 = prop_all_endophagic_indoors[[vv]][3,ee],
             col = '#222222', lwd = 1.5)
    rect(xleft = ee + offset[2 + 2 * (vv-1)] - (rectwidth/2), xright = ee + offset[2 + 2 * (vv-1)] + (rectwidth/2),
         ybottom = 0, ytop = prop_leaky_endophagic_indoors[[vv]][2,ee], col = palette[vv], density = 45, border = '#222222')
    segments(x0 = ee + offset[2 + 2 * (vv-1)], y0 = prop_leaky_endophagic_indoors[[vv]][1,ee], y1 = prop_leaky_endophagic_indoors[[vv]][3,ee],
             col = '#222222', lwd = 1.5)
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.20), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 3, line = 0, at = 0.525, font = 2, 'C')
mtext(side = 3, line = 0, at = mean(c(0.5, length(EIR_equil) + 0.5)),
      'Endophagic, Daytime Biting')
mtext(side = 2, line = 2.3, 'Treatment Failures (%)')

plot(NA, NA, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(vv in 1:nrow(scenarios_vector_control))
  {
    rect(xleft = ee + offset[1 + 2 * (vv-1)] - (rectwidth/2), xright = ee + offset[1 + 2 * (vv-1)] + (rectwidth/2),
         ybottom = 0, ytop = prop_all_exophagic[[vv]][2,ee], col = palette[vv])
    segments(x0 = ee + offset[1 + 2 * (vv-1)], y0 = prop_all_exophagic[[vv]][1,ee], y1 = prop_all_exophagic[[vv]][3,ee],
             col = '#222222', lwd = 1.5)
    rect(xleft = ee + offset[2 + 2 * (vv-1)] - (rectwidth/2), xright = ee + offset[2 + 2 * (vv-1)] + (rectwidth/2),
         ybottom = 0, ytop = prop_leaky_exophagic[[vv]][2,ee], col = palette[vv], density = 45, border = '#222222')
    segments(x0 = ee + offset[2 + 2 * (vv-1)], y0 = prop_leaky_exophagic[[vv]][1,ee], y1 = prop_leaky_exophagic[[vv]][3,ee],
             col = '#222222', lwd = 1.5)
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.20), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Treatment Failures (%)')
mtext(side = 3, line = 0, at = 0.525, font = 2, 'D')
mtext(side = 3, line = 0, at = mean(c(0.5, length(EIR_equil) + 0.5)),
      'Exophagic, No Preference in Biting Time')

dev.off()
