# install necessary packages
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}

# load necessary functions
source('functions_sample_sizes.R')

# specify parameters 
alpha_level = 0.05
power <- seq(from = 0.50, to = 0.95, by = 0.05)
efficacy <- 0.75
efficacy_lower <- 0.5
allocation_ratio <- 1

# load file with trial outputs in placebo arm
output <- read.csv('../../output/analysis/eir_vs_heterogeneity/output_files/all_or_none/sample_size_calculation.csv')

# get the range of eirs 
eir_equil <- sort(unique(output$EIR_equil))

# calculate necessary samples sizes
n_proportion <- lapply(eir_equil, function(eir){sapply(output$prop_relapse_desired[output$EIR_equil == eir], function(prop){sapply(power, function(p){calcSampleSizeRisk(alpha_level = alpha_level,
                                                                                                                                                                        power = p,
                                                                                                                                                                        prop_placebo = prop,
                                                                                                                                                                        efficacy = efficacy)['treatment']})})})

n_incidence <- lapply(eir_equil, function(eir){sapply(output$incidence_relapse_desired[output$EIR_equil == eir], function(irp){sapply(power, function(p){calcSampleSizeIncidence(alpha_level = alpha_level,
                                                                                                                                                                                 power = p,
                                                                                                                                                                                 incidence_rate_placebo = irp,
                                                                                                                                                                                 efficacy = efficacy)['treatment']})})})

n_cox <- lapply(eir_equil, function(eir){sapply(output$incidence_relapse_desired[output$EIR_equil == eir], function(irp){sapply(power, function(p){calcSampleSizeCox(alpha_level = alpha_level,
                                                                                                                                                                           power = p,
                                                                                                                                                                           allocation_ratio = allocation_ratio,   
                                                                                                                                                                           incidence_rate_placebo = irp,
                                                                                                                                                                           efficacy = efficacy)['treatment']})})})

# calculate the summaries 
n_proportion_summary <- lapply(n_proportion, function(df){apply(df, 1, summary)})
n_incidence_summary <- lapply(n_incidence, function(df){apply(df, 1, summary)})
n_cox_summary <- lapply(n_cox, function(df){apply(df, 1, summary)})

# generate plots 
palette <- brewer.pal(n = 9, name = 'YlOrRd')
palette <- palette[c(5,7,9)]
cex = 1.5

jpeg(filename = '../../output/figs/fig_S7.jpg', width = 8.5, height = 8.5/3, units = 'in', res = 500)
layout(mat = matrix(1:3, nrow = 1))
par(mar = c(3.3,3.6,1.3,0.8))
n_cox_max <- max(sapply(n_cox_summary, function(df){max(df['3rd Qu.',])}))
plot(NA, NA, xlim = c(min(power) - 0.025, max(power) + 0.025), ylim = c(0, 1.1 * n_cox_max), axes = F,
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '', las = 1)
for(ee in 1:length(n_cox_summary))
{
  polygon(x = c(power, rev(power)), y = c(n_cox_summary[[ee]]['1st Qu.',], rev(n_cox_summary[[ee]]['3rd Qu.',])),
          border = NA, col = col2alpha(palette[ee], alpha = 0.50))
  lines(power, n_cox_summary[[ee]]['Median',], lwd = cex, col = palette[ee])
  points(x = power, y = n_cox_summary[[ee]]['Median',], pch = 16, col = palette[ee], cex = cex)
}
box()
axis(side = 1, at = power, labels = power * 100)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Power (%)')
mtext(side = 2, line = 2.3, 'Sample Size (per arm)')
mtext(side = 3, line = 0, adj = 0, 'A', font = 2)
legend('topleft', pch = c(NA,NA,NA,NA,16,15), lwd = c(cex, cex, cex, cex, NA, NA),
       col = c(palette, '#222222', col2alpha('#222222', alpha = 0.75)), pt.cex = cex,
       legend = c('EIR = 1', 'EIR = 10', 'EIR = 100', 'Median', 'IQR'), bty = 'n')

n_incidence_max <- max(sapply(n_incidence_summary, function(df){max(df['3rd Qu.',])}))
plot(NA, NA, xlim = c(min(power) - 0.025, max(power) + 0.025), ylim = c(0, 1.1 * n_incidence_max), axes = F,
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '', las = 1)
for(ee in 1:length(n_incidence_summary))
{
  polygon(x = c(power, rev(power)), y = c(n_incidence_summary[[ee]]['1st Qu.',], rev(n_incidence_summary[[ee]]['3rd Qu.',])),
          border = NA, col = col2alpha(palette[ee], alpha = 0.50))
  lines(power, n_incidence_summary[[ee]]['Median',], lwd = cex, col = palette[ee])
  points(x = power, y = n_incidence_summary[[ee]]['Median',], pch = 16, col = palette[ee], cex = cex)
}
box()
axis(side = 1, at = power, labels = power * 100)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Power (%)')
mtext(side = 3, line = 0, adj = 0, 'B', font = 2)

n_proportion_max <- max(sapply(n_proportion_summary, function(df){max(df['3rd Qu.',])}))
plot(NA, NA, xlim = c(min(power) - 0.025, max(power) + 0.025), ylim = c(0, 1.1 * n_proportion_max), axes = F,
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '', las = 1)
for(ee in 1:length(n_proportion_summary))
{
  polygon(x = c(power, rev(power)), y = c(n_proportion_summary[[ee]]['1st Qu.',], rev(n_proportion_summary[[ee]]['3rd Qu.',])),
          border = NA, col = col2alpha(palette[ee], alpha = 0.50))
  lines(power, n_proportion_summary[[ee]]['Median',], lwd = cex, col = palette[ee])
  points(x = power, y = n_proportion_summary[[ee]]['Median',], pch = 16, col = palette[ee], cex = cex)
}
box()
axis(side = 1, at = power, labels = power * 100)
axis(side = 2, las = 1)
mtext(side = 1, line = 2.3, 'Power (%)')
mtext(side = 3, line = 0, adj = 0, 'C', font = 2)

dev.off()