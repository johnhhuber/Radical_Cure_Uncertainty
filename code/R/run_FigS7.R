# set working directory
setwd('~/Dropbox/Radical_Cure_MOA/code/')

# clear existing workspace
rm(list = ls())

# specify path to output files 
path_output_all_or_none <- '../output/analysis_0320/eir_vs_heterogeneity/output_files/all_or_none/'

# list all of the output files 
#files_output_all_or_none <- list.files(path = path_output_all_or_none, full.names = T, pattern = 'efficacy')

# load the files 
#output_all_or_none <- lapply(files_output_all_or_none, function(ff){read.csv(ff)})
#output <- as.data.frame((do.call('rbind', output_all_or_none)))

# specify parameter ranges
sigma_het <- seq(from = 0, to = 3, by = 1)
EIR_equil <- c(1, 10, 100)
n_rep <- 200

# construct parameter grid 
output <- expand.grid(sigma_het = sigma_het,
                           EIR_equil = EIR_equil)

output <- output[rep(seq_len(nrow(output)), each  = n_rep), ]
output$param_id <- 1:nrow(output) - 1

# subset to only look at homogeneous biting scenario 
output <- subset(output, sigma_het == 0)

# subset to only include homogeneous biting scenario 
#output <- subset(output, sig_het == 0, select = c('param_id', 'eir_equil'))

# add in the incidence rate of interest 
output$incidence_relapse_desired <- NA
output$incidence_relapse_all <- NA
output$incidence_recurrent <- NA
output$prop_relapse_desired <- NA
output$prop_relapse_all <- NA
output$prop_recurrent <- NA

# specify person years of followup
#person_year_followup <- 180 / 365

# loop through and calculate the incidence rates 
for(ii in 1:nrow(output))
{
  print(ii)
  
  # get param id 
  param_id <- output$param_id[ii] 
  
  if(file.exists(paste(path_output_all_or_none, 'indiv/indiv_', param_id, '.csv.bz2', sep = '')))
  {
    # load the files 
    output_participant <- read.csv(paste(path_output_all_or_none, 'indiv/indiv_', param_id, '.csv.bz2', sep = ''))
    output_recurrent <- read.csv(paste(path_output_all_or_none, 'recurr/recurr_', param_id, '.csv.bz2', sep = ''))
    
    # get the id's of placebo participants and subset recurrent infections to get that 
    id_placebo <- output_participant$Participant_ID[output_participant$Trial_Arm == 'PLACEBO']
    output_recurrent <- subset(output_recurrent, Participant_ID %in% id_placebo)
    output_recurrent <- output_recurrent[-match(id_placebo, output_recurrent$Participant_ID),]
    
    # calculate incidence rates 
    output$incidence_relapse_desired[ii] <- (sum(grepl('RELAPSE_LM_PRE_', output_recurrent$Infection_Type)) + sum(grepl('RELAPSE_D_PRE_', output_recurrent$Infection_Type))) / (length(id_placebo))
    output$incidence_relapse_all[ii] <- sum(grepl('RELAPSE', output_recurrent$Infection_Type)) / (length(id_placebo))
    output$incidence_recurrent[ii] <- nrow(output_recurrent) / (length(id_placebo))
    
    # calculate the proportion at risk 
    output$prop_relapse_desired[ii] <- length(unique(output_recurrent$Participant_ID[grepl('RELAPSE_LM_PRE_', output_recurrent$Infection_Type) | grepl('RELAPSE_D_PRE_', output_recurrent$Infection_Type)])) / length(id_placebo)
    output$prop_relapse_all[ii] <- length(unique(output_recurrent$Participant_ID[grepl('RELAPSE', output_recurrent$Infection_Type)])) / length(id_placebo)
    output$prop_recurrent[ii] <- length(unique(output_recurrent$Participant_ID)) / length(id_placebo)
  }
}

# write output to file 
write.csv(output, file = '../output/analysis_0320/eir_vs_heterogeneity/output_files/all_or_none/sample_size_calculation.csv', row.names = F)
