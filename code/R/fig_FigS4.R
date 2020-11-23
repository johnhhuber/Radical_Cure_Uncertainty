# install necessary packages
if(!require(ggsci)){install.packages('ggsci'); library(ggsci)}

# specify input file paths 
path_output_leaky <- '../output/analysis_0930/eir_vs_heterogeneity/output_files/leaky/'
path_output_all_or_none <- '../output/analysis_0930/eir_vs_heterogeneity/output_files/all_or_none/'

# specify parameter ranges
sigma_het <- seq(from = 0, to = 3, by = 1)
EIR_equil <- c(0.1, 1, 10, 100)
n_rep <- 200

# construct parameter grid 
param_sweep <- expand.grid(sigma_het = sigma_het,
                           EIR_equil = EIR_equil)

param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each  = n_rep), ]

# get the file indices 
file_indices <- lapply(EIR_equil, function(eir){lapply(sigma_het, function(sig){which(param_sweep$EIR_equil == eir & param_sweep$sigma_het == sig) - 1})})

# get the output files and calculate the summary statistics for the recurrent infections 
stats_reinfection_all <- list()
stats_reinfection_detectable <- list()
for(ii in 1:length(file_indices))
{
  print(ii)

  prop_recurr_all <- list()
  prop_recurr_detectable <- list()
  for(jj in 1:length(file_indices[[ii]]))
  {
    
    print(jj)
    vec_all <- c()
    vec_detectable <- c()
    for(ff in 1:length(file_indices[[ii]][[jj]]))
    {
      print(ff)
      # leaky intervention action 
      file_name_recurr <- paste(path_output_leaky, 'recurr/recurr_', file_indices[[ii]][[jj]][ff], '.csv.bz2', sep = '')
      file_name_indiv <- paste(path_output_leaky, 'indiv/indiv_', file_indices[[ii]][[jj]][ff], '.csv.bz2', sep = '')
      if(file.exists(file_name_recurr))
      {
        # load files
        indiv <- read.csv(file_name_indiv)
        recurr <- read.csv(file_name_recurr)
        
        # get number of participants
        num_participants <- nrow(indiv)
        
        # subset recurrent infection outputs to remove infection that caused enrollment and only to consider reinfection events
        recurr <- recurr[-match(indiv$Participant_ID, recurr$Participant_ID),]
        
        # calculate the proportion of recurrent infections that were detectable 
        recurr <- subset(recurr, grepl('REINFECTION', Infection_Type))
        
        # calculate the proportion reinfected 
        reinfected <- sapply(indiv$Participant_ID, function(id){id %in% recurr$Participant_ID})
        vec_all <- c(vec_all, sum(reinfected) / num_participants)
        
        # calculate the proportion with a detectable reinfection 
        detectable <- subset(recurr, !grepl('PCR', Infection_Type))
        reinfected <- sapply(indiv$Participant_ID, function(id){id %in% detectable$Participant_ID})
        vec_detectable <- c(vec_detectable, sum(reinfected) / num_participants)
      }
      
      # all-or-none intervention action 
      file_name_recurr <- paste(path_output_all_or_none, 'recurr/recurr_', file_indices[[ii]][[jj]][ff], '.csv.bz2', sep = '')
      file_name_indiv <- paste(path_output_all_or_none, 'indiv/indiv_', file_indices[[ii]][[jj]][ff], '.csv.bz2', sep = '')
      if(file.exists(file_name_recurr))
      {
        # load files
        indiv <- read.csv(file_name_indiv)
        recurr <- read.csv(file_name_recurr)
        
        # subset recurrent infection outputs to remove infection that caused enrollment and only to consider reinfection events
        recurr <- recurr[-match(indiv$Participant_ID, recurr$Participant_ID),]
        
        # calculate the proportion of recurrent infections that were detectable 
        recurr <- subset(recurr, grepl('REINFECTION', Infection_Type))
        
        # calculate the proportion reinfected 
        reinfected <- sapply(indiv$Participant_ID, function(id){id %in% recurr$Participant_ID})
        vec_all <- c(vec_all, sum(reinfected) / num_participants)
        
        # calculate the proportion with a detectable reinfection 
        detectable <- subset(recurr, !grepl('PCR', Infection_Type))
        reinfected <- sapply(indiv$Participant_ID, function(id){id %in% detectable$Participant_ID})
        vec_detectable <- c(vec_detectable, sum(reinfected) / num_participants)
      }
    }
    prop_recurr_all[[jj]] = vec_all
    prop_recurr_detectable[[jj]] = vec_detectable
  }
  stats_reinfection_all[[ii]] <- prop_recurr_all
  stats_reinfection_detectable[[ii]] <- prop_recurr_detectable
}

# generate plot 
palette <- pal_material(palette = 'deep-orange', n = 8)(8)[-1]
offset <- seq(from = -0.4, to = 0.4, length.out = length(sigma_het) * 2)

jpeg(filename = '../output/figs_manuscript/fig_S4.jpg', width = 8, height = 5, units = 'in', res = 500)
par(mar = c(3.6, 3.6, 0.8, 0.8))
plot(NA, NA, xlim = c(0.5, length(EIR_equil) + 0.5), ylim = c(-0.02,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
for(ee in 1:length(EIR_equil))
{
  for(ss in 1:length(sigma_het))
  {
    
    quantiles_all = quantile(stats_reinfection_all[[ee]][[ss]], probs = c(0.25, 0.50, 0.75))
    segments(x0 = ee + offset[1 + 2 * (ss-1)], y0 = quantiles_all[1], y1 = quantiles_all[3], col = palette[ss], lwd = 2)
    points(ee + offset[1 + 2 * (ss - 1)], quantiles_all[2], pch = 15, cex = 1.5, col = palette[ss])
    
    quantiles_detect = quantile(stats_reinfection_detectable[[ee]][[ss]], probs = c(0.25, 0.50, 0.75))
    segments(x0 = ee + offset[2 + 2 * (ss-1)], y0 = quantiles_detect[1], y1 = quantiles_detect[3], col = palette[ss], lwd = 2)
    points(ee + offset[2 + 2 * (ss - 1)], quantiles_detect[2], pch = 16, cex = 1.5, col = palette[ss])
  }
}
abline(v = seq(from = 0.5, to = length(EIR_equil) + 0.5, by = 1))
box()
axis(side = 1, at = 1:length(EIR_equil), labels = EIR_equil)
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 1, line = 2.3, 'EIR')
mtext(side = 2, line = 2.3, 'Percentage of Trial Participants')
legend('topleft', pch = c(15,16),
       col = '#222222',
       cex = 0.8,
       legend = c('Any Reinfection', 'Detectable Reinfection'),
       bty = 'n', pt.cex = 1.5)

legend(x = 0.5, y = 0.875, pch = 15,
       ncol = length(palette),
       col = palette,
       legend = sigma_het,
       title = 'Heterogeneity in Biting:',
       xjust = 0,
       cex = 0.8,
       title.adj = 0.05,
       bty = 'n', pt.cex = 1.5)
dev.off()


