# set working directory
setwd('~/Dropbox/Radical_Cure_MOA/code/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(ggsci)){install.packages('ggsci'); library(ggsci)}

# load necessary functions
source('functions_sample_sizes.R')

# specify parameters 
alpha_level = 0.05
power <- seq(from = 0.50, to = 0.95, by = 0.05)
efficacy <- 0.75
efficacy_lower <- 0.5
allocation_ratio <- 1
incidence_rate_placebo <- c(0.05, 0.1, 0.5, 1,10)

# generate the sample size calculation 
n_incidence <- lapply(incidence_rate_placebo, function(irp){sapply(power, function(p){as.numeric(calcSampleSizeIncidence(alpha_level = alpha_level,
                                                                                                                         power = p,
                                                                                                                         incidence_rate_placebo = irp,
                                                                                                                         efficacy = efficacy)['treatment'])})})

n_cox <- lapply(incidence_rate_placebo, function(irp){sapply(power, function(p){as.numeric(calcSampleSizeCox(alpha_level = alpha_level,
                                                                                                             power = p,
                                                                                                             allocation_ratio = allocation_ratio,
                                                                                                             incidence_rate_placebo = irp,
                                                                                                             efficacy = efficacy)['treatment'])})})


# generate plot 
palette <- pal_material(palette = 'orange', n = 10)(10)[seq(from = 2, to = 10, length.out = length(incidence_rate_placebo))]

jpeg(filename = '../output/figs_manuscript/fig_S5.jpg', width = 6.5, height = 6.5, units = 'in', res = 500)

layout(mat = matrix(1:2, nrow = 2))
par(mar = c(3.3, 4.3, 0.8, 0.8))
plot(NA, NA, xlim = c(0.45,1), ylim = c(0, 800), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
for(ii in 1:length(n_incidence))
{
  lines(power, n_incidence[[ii]], lwd = 1.5, col = palette[ii])
  points(power, n_incidence[[ii]], pch = 21, col = 'white', bg = palette[ii], cex = 1.5)
}
box()
axis(side = 1, at = power, labels = power * 100)
axis(side = 2, las = 1)
mtext(side = 2, line = 3.3, 'Sample Size in Each Arm')
mtext(side = 3, adj = 0, font = 2, 'A')
legend('topleft', lwd = 2, col = palette, legend = incidence_rate_placebo, bty = 'n', title = 'Incidence Rate')

plot(NA, NA, xlim = c(0.45,1), ylim = c(0, 800), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
for(ii in 1:length(n_cox))
{
  lines(power, n_cox[[ii]], lwd = 1.5, col = palette[ii])
  points(power, n_cox[[ii]], pch = 21, col = 'white', bg = palette[ii], cex = 1.5)
}
box()
axis(side = 1, at = power, labels = power * 100)
axis(side = 2, las = 1)
mtext(side = 2, line = 3.3, 'Sample Size in Each Arm')
mtext(side = 3, adj = 0, font = 2, 'B')
mtext(side = 1, line = 2.3, 'Power (%)')

dev.off()