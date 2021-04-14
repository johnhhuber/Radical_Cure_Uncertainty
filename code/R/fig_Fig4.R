# install necessary packages
if(!require(ggsci)){install.packages('ggsci'); library(ggsci)}

# specify the output path
path_output_all_or_none <- '../../output/analysis/vector_control/output_files/all_or_none/efficacy/'

# list all of the output files 
files_output_all_or_none <- list.files(path = path_output_all_or_none, full.names = T, pattern = 'efficacy')

# load the files 
output_all_or_none <- lapply(files_output_all_or_none, function(ff){read.csv(ff)})
output_efficacy_all_or_none <- as.data.frame((do.call('rbind', output_all_or_none)))

# specify parameter ranges 
EIR_equil <- c(1 / 365, 10 / 365, 100 / 365)
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

# subset to the default scenario 
output_efficacy_all_default <- subset(output_efficacy_all_or_none, PSI_indoors == 0.9 & PSI_bed == 0.45)

# specify vector control scenarios and range of eir values considered 
scenarios_vector_control <- data.frame(is_LLIN_distributed = c(0,1,0,1),
                                       is_IRS_administered = c(0,0,1,1))

eir_equil <- sort(unique(output_efficacy_all_or_none$eir_equil))

# calculate efficacy
eff_all_default <- lapply(1:nrow(scenarios_vector_control), function(vc){sapply(eir_equil, function(eir){quantile(output_efficacy_all_default$eff_cph_recurrent_LM[output_efficacy_all_default$eir_equil == eir & 
                                                                                                                                                                     output_efficacy_all_default$is_LLIN_distributed == scenarios_vector_control$is_LLIN_distributed[vc] & 
                                                                                                                                                                     output_efficacy_all_default$is_IRS_administered == scenarios_vector_control$is_IRS_administered[vc]],
                                                                                                                  probs = c(0.25, 0.50, 0.75), na.rm = T)})})
# generate plots 
palette <- head(pal_lancet()(9), n = 4)
offset <- seq(from = -0.4, to = 0.4, length.out = nrow(scenarios_vector_control))
cex = 1.5

jpeg(filename = '../../output/figs/fig_4.jpg', width = 8, height = 4, units = 'in', res = 500)
par(mar = c(3.6, 3.6, 0.8, 0.8))
plot(NA, NA, xlim = c(0.5, length(eir_equil) + 0.5), ylim = c(0, 1.0), axes = F, 
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
abline(h = 0.75, lwd = 1, lty = 2)
abline(v = seq(from = 0.5, to = length(eir_equil) + 0.5, by = 1), col = '#222222')
for(ee in 1:length(eir_equil))
{
  for(vv in 1:nrow(scenarios_vector_control))
  {
    segments(x0 = ee + offset[1 + 1 * (vv - 1)], y0 = eff_all_default[[vv]][1,ee], y1 = eff_all_default[[vv]][3,ee],
             lwd = cex, col = palette[vv])
    points(ee + offset[1 + 1 * (vv - 1)], eff_all_default[[vv]][2,ee], pch = 16, cex = cex, col = palette[vv])
  }
}
box()
axis(side = 1, at = 1:length(eir_equil), labels = eir_equil)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Efficacy')
legend('bottomleft', pch = c(NA, rep(15,4)), lwd = c(1, NA, NA, NA, NA),
       lty = c(2, NA, NA, NA, NA), pt.cex = 1.5, col = c('#222222', palette),
       legend = c('Clearance Probability', 'No Vector Control', 'LLINs only', 'IRS only', 'LLINs and IRS'),
       bty = 'n', cex = 0.8)

dev.off()



