# load necessary functions
source('functions_efficacy.R')

# specify path in 
path_in = 'wkdir/heterogeneous_protection/'

# list output files 
files_output_trial <- list.files(path = paste(path_in, 'trial', sep = ''), full.names = T)
files_output_recurrent <- list.files(path = paste(path_in, 'recurrent', sep = ''), full.names = T)
files_output_participants <- list.files(path = paste(path_in, 'participants', sep = ''), full.names = T)

# specify the ranges 
PQ_prop_stratum_1 = seq(from = 0, to = 1, by = 0.25)
PQ_eff_stratum_1 <- seq(from = 0, to = 1, by = 0.25)
PQ_eff_stratum_2 <- seq(from = 0, to = 1, by = 0.25)
EIR_equil <- c(1e-5, 1e-4, 1e-3)

# create parameter grid 
param_sweep <- expand.grid(PQ_prop_stratum_1 = PQ_prop_stratum_1,
                           PQ_eff_stratum_1 = PQ_eff_stratum_1,
                           PQ_eff_stratum_2 = PQ_eff_stratum_2,
                           EIR_equil = EIR_equil)

param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each = 10), ]
n_sim <- 1250

# create data frame 
df <- data.frame(param_id = rep(NA, n_sim),
                 eir_equil = rep(NA, n_sim),
                 PQ_prop_stratum_1 = rep(NA, n_sim),
                 PQ_eff_stratum_1 = rep(NA, n_sim),
                 PQ_eff_stratum_2 = rep(NA, n_sim),
                 eff_cph = rep(NA, n_sim),
                 eff_any_PCR = rep(NA, n_sim),
                 eff_any_LM = rep(NA, n_sim),
                 eff_any_D = rep(NA, n_sim),
                 eff_relapse_PCR = rep(NA, n_sim),
                 eff_relapse_LM = rep(NA, n_sim),
                 eff_relapse_D = rep(NA, n_sim))

# go through and load 
for(ii in 1:1250)
{
  # get the id 
  id = ii
  #id <- as.numeric(head(strsplit(tail(strsplit(files_output_trial[ii], '_')[[1]], n = 1), '[.]')[[1]], n = 1))
  df$param_id[ii] <- id 
  
  # get parameter information 
  df$eir_equil <- param_sweep[id + 1, 'EIR_equil']
  df$PQ_prop_stratum_1 <- param_sweep[id + 1, 'PQ_prop_stratum_1']
  df$PQ_eff_stratum_1 <- param_sweep[id + 1, 'PQ_eff_stratum_1']
  df$PQ_eff_stratum_2 <- param_sweep[id + 1, 'PQ_eff_stratum_2']
  
  # load the files
  output_trial <- read.csv(grep(paste('_', id, sep = ''), files_output_trial, value = T))
  output_recurrent <- read.csv(grep(paste('_', id, sep = ''), files_output_recurrent, value = T))
  output_participants <- read.csv(grep(paste('_', id, sep = ''), files_output_participant, value = T))
  
  # organize the trial data 
  trial_data <- createEffDataFrame(output_participant = output_participants,
                                   output_recurrent = output_recurrent,
                                   output_trial = output_trial)
  
  # calculate efficacy metrics 
  df$eff_cph <- calcCoxEff(df)
  incid_eff <- calcIncidEff(df)
  df$eff_any_PCR <- incid_eff['num_any_PCR']
  df$eff_any_LM <- incid_eff['num_any_LM']
  df$eff_any_D <- incid_eff['num_any_D']
  df$eff_relapse_PCR <- incid_eff['num_relapse_PCR']
  df$eff_relapse_LM <- incid_eff['num_relapse_LM']
  df$eff_relapse_D <- incid_eff['num_relapse_D']
}

# write to file 
write.csv(df, file = 'efficacy_heterogeneous_protection.csv', row.names = F)
