# load necessary functions 
source('functions_efficacy.R')

# install necessary packages 
if(!require(doParallel)){install.packages('doParallel'); library(doParallel)}
if(!require(foreach)){install.packages('foreach'); library(foreach)}

# specify path in 
args <- commandArgs(trailingOnly = T)
path_input_model_params = args[1]
path_input_trial_design = args[2]
path_output = args[3]
batch_id = args[4]

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# list input and output data files 
files_input_params <- list.files(path = path_input_model_params, pattern = 'model_params', full.names = T)
files_input_trial <- list.files(path = path_input_trial_design, pattern = 'trial_design', full.names = T)

files_output_participants <- list.files(path = paste(path_output, 'indiv/', sep = ''), pattern = 'indiv', full.names = T)
files_output_recurrent <- list.files(path = paste(path_output, 'recurr/', sep = ''), pattern = 'recurr', full.names = T)
files_output_trial <- list.files(path = paste(path_output, 'trial/', sep = ''), pattern = 'trial', full.names = T)

# create a data frame to store output 
n_sims <- length(files_input_params)
print(n_sims)

batch_size = round(n_sims / 50)
batches = split(1:n_sims, ceiling(seq_along(1:n_sims) / batch_size))

# loop through and create the necessary files 
df_final <- foreach(ss = batches[[batch_id]], .combine = rbind, .packages = c('survival', 'icenReg')) %dopar%{
  # create temporary dataframe 
  df <- data.frame(param_id = rep(NA, 1),
                   eir_equil = rep(NA, 1),
                   sig_het = rep(NA, 1),
                   ff = rep(NA, 1),
                   trial_duration = rep(NA, 1),
                   trial_PQ_eff_stratum_1 = rep(NA, 1),
                   trial_PQ_eff_stratum_2 = rep(NA, 1),
                   trial_PQ_prop_stratum_1 = rep(NA, 1),
                   trial_PQ_prop_stratum_2 = rep(NA, 1),
                   trial_PQ_proph = rep(NA, 1),
                   is_LLIN_distributed = rep(NA, 1),
                   is_IRS_administered = rep(NA, 1),
                   eff_cph_recurrent_LM = rep(NA, 1),
                   eff_cph_recurrent_PCR = rep(NA, 1),
                   eff_cph_recurrent_D = rep(NA, 1),
                   eff_cph_any_relapse_LM = rep(NA, 1),
                   eff_cph_any_relapse_PCR = rep(NA, 1),
                   eff_cph_relapse_LM_trial = rep(NA, 1),
                   eff_cph_relapse_PCR_trial = rep(NA, 1),
                   eff_cph_relapse_D_trial = rep(NA, 1),
                   eff_cph_relapse_LM_all = rep(NA, 1),
                   eff_cph_relapse_PCR_all = rep(NA, 1),
                   eff_cph_relapse_D_all = rep(NA, 1),
                   eff_incid_recurrent_PCR = rep(NA, 1),
                   eff_incid_recurrent_LM = rep(NA, 1),
                   eff_incid_recurrent_D = rep(NA, 1),
                   eff_incid_any_relapse_PCR = rep(NA, 1),
                   eff_incid_any_relapse_LM = rep(NA, 1),
                   eff_incid_any_relapse_D = rep(NA, 1),
                   eff_incid_relapse_PCR = rep(NA, 1),
                   eff_incid_relapse_LM = rep(NA, 1),
                   eff_incid_relapse_D = rep(NA, 1),
                   eff_incid_trial_PCR = rep(NA, 1),
                   eff_incid_trial_LM = rep(NA, 1),
                   eff_incid_trial_D = rep(NA, 1),
                   eff_risk_relapse_PCR = rep(NA, 1),
                   eff_risk_relapse_LM = rep(NA, 1),
                   eff_risk_relapse_D = rep(NA, 1),
                   eff_risk_trial_PCR = rep(NA, 1),
                   eff_risk_trial_LM = rep(NA, 1),
                   eff_risk_trial_D = rep(NA, 1))
  
  # read in the parameters file 
  model_params <- read.table(file = files_input_params[ss])
  
  # get the simulation id 
  param_id <- as.numeric(head(strsplit(tail(strsplit(files_input_params[ss], '_')[[1]], n = 1), '[.]')[[1]], n = 1))
  
  # load the trial design file 
  trial_design <- read.csv(file = grep(pattern = paste('_', param_id, '.txt', sep = ''), files_input_trial, value = T), header = T)
  
  # fill in the parameters 
  df$param_id <- param_id
  df$eir_equil <- model_params[model_params[,1] == 'EIR_equil', 2] * 365
  df$sig_het <- model_params[model_params[,1] == 'sig_het', 2] 
  df$ff <- model_params[model_params[,1] == 'ff', 2]
  df$trial_duration <- trial_design[trial_design[,1] == 'trial_duration', 2]
  df$trial_PQ_eff_stratum_1 <- trial_design[trial_design[,1] == 'trial_PQ_eff_stratum_1', 2]
  df$trial_PQ_eff_stratum_2 <- trial_design[trial_design[,1] == 'trial_PQ_eff_stratum_2', 2]
  df$trial_PQ_prop_stratum_1 <- trial_design[trial_design[,1] == 'trial_PQ_prop_stratum_1', 2]
  df$trial_PQ_prop_stratum_2 <- trial_design[trial_design[,1] == 'trial_PQ_prop_stratum_2', 2]
  df$trial_PQ_proph <- trial_design[trial_design[,1] == 'trial_PQ_proph', 2]
  df$is_LLIN_distributed <- trial_design[trial_design[,1] == 'is_LLIN_distributed', 2]
  df$is_IRS_administered <- trial_design[trial_design[,1] == 'is_IRS_administered', 2]
  
  # load output files
  output_participants <- read.csv(file = grep(pattern = paste('_', param_id, '.csv', sep = ''), files_output_participants, value = T), header = T)
  output_recurrent <- read.csv(file = grep(pattern = paste('_', param_id, '.csv', sep = ''), files_output_recurrent, value = T), header = T)
  output_trial <- read.csv(file = grep(pattern = paste('_', param_id, '.csv', sep = ''), files_output_trial, value = T), header = T)
  
  # organize the trial data 
  trial_data <- createEffDataFrame(output_participant = output_participants,
                                   output_recurrent = output_recurrent,
                                   output_trial = output_trial)
  
  # get the incidence efficacy calculations 
  eff_incid <- calcIncidEff(trial_data)
  
  df$eff_incid_recurrent_PCR <- eff_incid['num_any_PCR']
  df$eff_incid_recurrent_LM <- eff_incid['num_any_LM']
  df$eff_incid_recurrent_D <- eff_incid['num_any_D']
  
  df$eff_incid_any_relapse_PCR <- eff_incid['num_any_relapse_PCR']
  df$eff_incid_any_relapse_LM <- eff_incid['num_any_relapse_LM']
  df$eff_incid_any_relapse_D <- eff_incid['num_any_relapse_D']
  
  df$eff_incid_relapse_PCR <- eff_incid['num_relapse_PCR']
  df$eff_incid_relapse_LM <- eff_incid['num_relapse_LM']
  df$eff_incid_relapse_D <- eff_incid['num_relapse_D']

  df$eff_incid_trial_PCR <- eff_incid['num_trial_PCR']
  df$eff_incid_trial_LM <- eff_incid['num_trial_LM']
  df$eff_incid_trial_D <- eff_incid['num_trial_D']
  
  # get the cph efficacy calcuations 
  eff_cph <- calcCoxEff(trial_data, days_followup = unique(output_trial$Days_Since_Enrollment))
  
  df$eff_cph_recurrent_LM <- eff_cph['cph_recurrent_LM']
  df$eff_cph_recurrent_PCR <- eff_cph['cph_recurrent_PCR']
  df$eff_cph_recurrent_D <- eff_cph['cph_recurrent_D']
  df$eff_cph_any_relapse_LM <- eff_cph['cph_any_relapse_LM']
  df$eff_cph_any_relapse_PCR <- eff_cph['cph_any_relapse_PCR']
  df$eff_cph_relapse_LM_trial <- eff_cph['cph_relapse_LM_trial']
  df$eff_cph_relapse_PCR_trial <- eff_cph['cph_relapse_PCR_trial']
  df$eff_cph_relapse_D_trial <- eff_cph['cph_relapse_D_trial']
  df$eff_cph_relapse_LM_all <- eff_cph['cph_relapse_LM_all']
  df$eff_cph_relapse_PCR_all <- eff_cph['cph_relapse_PCR_all']
  df$eff_cph_relapse_D_all <- eff_cph['cph_relapse_D_all']
  
  # get the risk-based efficacy calculations 
  eff_risk <- calcRiskEff(trial_data)
  
  df$eff_risk_relapse_PCR <- eff_risk['risk_relapse_PCR']
  df$eff_risk_relapse_LM <- eff_risk['risk_relapse_LM']
  df$eff_risk_relapse_D <- eff_risk['risk_relapse_D']

  df$eff_risk_trial_PCR <- eff_risk['risk_trial_PCR']
  df$eff_risk_trial_LM <- eff_risk['risk_trial_LM']
  df$eff_risk_trial_D <- eff_risk['risk_trial_D']
  
  # return 
  df
}

# write to file 
write.csv(df_final, file = paste(path_output, 'output_efficacy_', batch_id, '.csv', sep = ''), row.names = F)
