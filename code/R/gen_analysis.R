# load necessary functions 
source('functions_input_files.R')

##### Analysis 1 - How does estimated efficacy vary with EIR and heterogeneity in biting? #####

# specify path out 
path_out = '../../output/analysis/eir_vs_heterogeneity/input_files/'

# specify parameter ranges
sigma_het <- seq(from = 0, to = 3, by = 1)
EIR_equil <- c(0.1 / 365, 1 / 365, 10 / 365, 100 / 365)
n_rep <- 200

# construct parameter grid 
param_sweep <- expand.grid(sigma_het = sigma_het,
                           EIR_equil = EIR_equil)

param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each  = n_rep), ]

# write parameters to file 
for(ff in 1:nrow(param_sweep))
{
  print(ff)
  
  # trial design file 
  genTrialFile(output_file_path = paste(path_out, 'all_or_none/trial_design_', ff - 1, '.txt', sep = ''),
               treatment_arm_sample_size = 1000,
               placebo_arm_sample_size = 1000,
               recruitment_start_date = 739125,
               trial_duration = 180,
               dropout_rate = 0,
               trial_PQ_eff_stratum_1 = 1,
               trial_PQ_eff_stratum_2 = 0,
               trial_PQ_prop_stratum_1 = 0.75,
               trial_PQ_prop_stratum_2 = 0.25,
               trial_PQ_lowage = 16 * 365,
               trial_PQ_proph = 28,
               trial_observation_period = 3,
               is_LLIN_distributed = 1,
               is_IRS_administered = 0,
               G6PD_activity_threshold = 7,
               followup_dates = c(0,1,2,3,8,15,22,29,60,120,180))
  
  genTrialFile(output_file_path = paste(path_out, 'leaky/trial_design_', ff - 1, '.txt', sep = ''),
               treatment_arm_sample_size = 1000,
               placebo_arm_sample_size = 1000,
               recruitment_start_date = 739125,
               trial_duration = 180,
               dropout_rate = 0,
               trial_PQ_eff_stratum_1 = 0.75,
               trial_PQ_eff_stratum_2 = 0,
               trial_PQ_prop_stratum_1 = 1,
               trial_PQ_prop_stratum_2 = 0,
               trial_PQ_lowage = 16 * 365,
               trial_PQ_proph = 28,
               trial_observation_period = 3,
               is_LLIN_distributed = 1,
               is_IRS_administered = 0,
               G6PD_activity_threshold = 7,
               followup_dates = c(0,1,2,3,8,15,22,29,60,120,180))
  
  # model parameters file 
  genModelParamsFile(output_file_path = paste(path_out, 'all_or_none/model_params_', ff - 1, '.txt', sep = ''),
                     N_part = 100000,
                     EIR_equil = param_sweep[ff, 'EIR_equil'],
                     sig_het = param_sweep[ff, 'sigma_het'],
                     start_time = 2019,
                     end_time = 2040)
  
  genModelParamsFile(output_file_path = paste(path_out, 'leaky/model_params_', ff - 1, '.txt', sep = ''),
                     N_part = 100000,
                     EIR_equil = param_sweep[ff, 'EIR_equil'],
                     sig_het = param_sweep[ff, 'sigma_het'],
                     start_time = 2019,
                     end_time = 2040)
}


##### Analysis 2 - How does the estimated efficacy vary with relapse rate and the duration of follow-up? #####

# specify path out 
path_out = '../../output/analysis/followup_vs_relapse/input_files/'

# specify parameter ranges 
EIR_equil <- c(0.1 / 365, 1 / 365, 10 / 365, 100 / 365)
trial_duration <- c(90, 180, 365, 730)
ff <- c(1/30, 1/60, 1/90, 1/180)


# construct parameter grid 
param_sweep <- expand.grid(EIR_equil = EIR_equil,
                           trial_duration = trial_duration,
                           ff = ff)
param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each = n_rep), ]

# specify the follow up dates for each duration of follow up 
followup_dates <- list()
followup_dates[[1]] <- c(0,1,2,3,8,15,22,29,60, 90)
followup_dates[[2]] <- c(followup_dates[[1]], 120, 180)
followup_dates[[3]] <- c(followup_dates[[2]], 240, 300, 365)
followup_dates[[4]] <- c(followup_dates[[2]], 240, 300, 360, 420, 480, 540, 600, 660, 730)

# write parameters to file 
for(ff in 1:nrow(param_sweep))
{
  print(ff)
  
  # trial design file 
  genTrialFile(output_file_path = paste(path_out, 'all_or_none/trial_design/trial_design_', ff - 1, '.txt', sep = ''),
               treatment_arm_sample_size = 1000,
               placebo_arm_sample_size = 1000,
               recruitment_start_date = 739125,
               trial_duration = param_sweep[ff, 'trial_duration'],
               dropout_rate = 0,
               trial_PQ_eff_stratum_1 = 1,
               trial_PQ_eff_stratum_2 = 0,
               trial_PQ_prop_stratum_1 = 0.75,
               trial_PQ_prop_stratum_2 = 0.25,
               trial_PQ_lowage = 16 * 365,
               trial_PQ_proph = 28,
               trial_observation_period = 3,
               is_LLIN_distributed = 1,
               is_IRS_administered = 0,
               G6PD_activity_threshold = 7,
               followup_dates = followup_dates[[which(trial_duration == param_sweep[ff, 'trial_duration'])]])
  
  genTrialFile(output_file_path = paste(path_out, 'leaky/trial_design/trial_design_', ff - 1, '.txt', sep = ''),
               treatment_arm_sample_size = 1000,
               placebo_arm_sample_size = 1000,
               recruitment_start_date = 739125,
               trial_duration = param_sweep[ff, 'trial_duration'],
               dropout_rate = 0,
               trial_PQ_eff_stratum_1 = 0.75,
               trial_PQ_eff_stratum_2 = 0,
               trial_PQ_prop_stratum_1 = 1,
               trial_PQ_prop_stratum_2 = 0,
               trial_PQ_lowage = 16 * 365,
               trial_PQ_proph = 28,
               trial_observation_period = 3,
               is_LLIN_distributed = 1,
               is_IRS_administered = 0,
               G6PD_activity_threshold = 7,
               followup_dates = followup_dates[[which(trial_duration == param_sweep[ff, 'trial_duration'])]])
  
  # model parmaeters file 
  genModelParamsFile(output_file_path = paste(path_out, 'all_or_none/model_params/model_params_', ff - 1, '.txt', sep = ''),
                     N_part = 100000,
                     EIR_equil = param_sweep[ff, 'EIR_equil'],
                     sig_het = 0,
                     ff = param_sweep[ff, 'ff'],
                     start_time = 2019,
                     end_time = 2040)
  
  genModelParamsFile(output_file_path = paste(path_out, 'leaky/model_params/model_params_', ff - 1, '.txt', sep = ''),
                     N_part = 100000,
                     EIR_equil = param_sweep[ff,'EIR_equil'],
                     sig_het = 0,
                     ff = param_sweep[ff, 'ff'],
                     start_time = 2019,
                     end_time = 2040)
}

##### Analysis 3 - How does vector control help to correct biases? ##### 

# specify path out 
path_out = '../../output/analysis/vector_control/input_files/'

# specify parameter ranges 
EIR_equil <- c(0.1 / 365, 1 / 365, 10 / 365, 100 / 365)
is_LLIN_distributed <- c(0,1)
is_IRS_administered <- c(0,1)

PSI_indoors <- c(0.9, 0.9, 0.9, 0.1)
PSI_bed <- c(0.9 * 0.5, 0.9 * 0.75, 0.9 * 0.25, 0.1 * 0.5)

# construct parameter grid 
param_sweep <- expand.grid(EIR_equil = EIR_equil,
                           is_LLIN_distributed = is_LLIN_distributed,
                           is_IRS_administered = is_IRS_administered)
param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each = length(PSI_indoors)),]
param_sweep$PSI_indoors <- rep(PSI_indoors, nrow(param_sweep) / length(PSI_indoors))
param_sweep$PSI_bed <- rep(PSI_bed, nrow(param_sweep) / length(PSI_indoors))

param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each = n_rep), ]

# write parameters to file 
for(ff in 1:nrow(param_sweep))
{
  print(ff)
  
  # trial design file 
  genTrialFile(output_file_path = paste(path_out, 'all_or_none/trial_design/trial_design_', ff - 1, '.txt', sep = ''),
               treatment_arm_sample_size = 1000,
               placebo_arm_sample_size = 1000,
               recruitment_start_date = 739125,
               trial_duration = 180,
               dropout_rate = 0,
               trial_PQ_eff_stratum_1 = 1,
               trial_PQ_eff_stratum_2 = 0,
               trial_PQ_prop_stratum_1 = 0.75,
               trial_PQ_prop_stratum_2 = 0.25,
               trial_PQ_lowage = 16 * 365,
               trial_PQ_proph = 28,
               trial_observation_period = 3,
               is_LLIN_distributed = param_sweep[ff, 'is_LLIN_distributed'],
               is_IRS_administered = param_sweep[ff, 'is_IRS_administered'],
               G6PD_activity_threshold = 7,
               followup_dates = c(0,1,2,3,8,15,22,29,60,120,180))
  
  
  genTrialFile(output_file_path = paste(path_out, 'leaky/trial_design/trial_design_', ff - 1, '.txt', sep = ''),
               treatment_arm_sample_size = 1000,
               placebo_arm_sample_size = 1000,
               recruitment_start_date = 739125,
               trial_duration = 180,
               dropout_rate = 0,
               trial_PQ_eff_stratum_1 = 0.75,
               trial_PQ_eff_stratum_2 = 0,
               trial_PQ_prop_stratum_1 = 1,
               trial_PQ_prop_stratum_2 = 0,
               trial_PQ_lowage = 16 * 365,
               trial_PQ_proph = 28,
               trial_observation_period = 3,
               is_LLIN_distributed = param_sweep[ff, 'is_LLIN_distributed'],
               is_IRS_administered = param_sweep[ff, 'is_IRS_administered'],
               G6PD_activity_threshold = 7,
               followup_dates = c(0,1,2,3,8,15,22,29,60,120,180))
  
  # model parameters file
  genModelParamsFile(output_file_path = paste(path_out, 'all_or_none/model_params/model_params_', ff - 1, '.txt', sep = ''),
                     N_part = 100000,
                     EIR_equil = param_sweep[ff, 'EIR_equil'],
                     sig_het = 0,
                     start_time = 2019,
                     end_time = 2040)
  
  genModelParamsFile(output_file_path = paste(path_out, 'leaky/model_params/model_params_', ff - 1, '.txt', sep = ''),
                     N_part = 100000,
                     EIR_equil = param_sweep[ff, 'EIR_equil'],
                     sig_het = 0,
                     start_time = 2019,
                     end_time = 2040)
  
  # mosquito file 
  genMosqFile(output_file_path = paste(path_out, 'all_or_none/mosquito/mosquito_', ff - 1, '.txt', sep = ''),
              PSI_indoors = param_sweep[ff,'PSI_indoors'],
              PSI_bed = param_sweep[ff, 'PSI_bed'])
  
  genMosqFile(output_file_path = paste(path_out, 'leaky/mosquito/mosquito_', ff - 1, '.txt', sep = ''),
              PSI_indoors = param_sweep[ff,'PSI_indoors'],
              PSI_bed = param_sweep[ff, 'PSI_bed'])
}

###### Analysis 4 - How do the estimated efficacy vary with radical cure therapeutic? #####

# specify path out 
path_out = '../../output/analysis/radical_cure_therapeutic/input_files/'

# specify parameter ranges 
EIR_equil <- c(0.1 / 365, 1 / 365, 10 / 365, 100 / 365)
trial_PQ_proph <- c(28, 45)

# construct parameter grid 
param_sweep <- expand.grid(EIR_equil = EIR_equil,
                           trial_PQ_proph = trial_PQ_proph)

param_sweep <- param_sweep[rep(seq_len(nrow(param_sweep)), each = n_rep), ]

# write parameters to file 
for(ff in 1:nrow(param_sweep))
{
  print(ff)
  
  # trial design file 
  genTrialFile(output_file_path = paste(path_out, 'all_or_none/trial_design_', ff - 1, '.txt', sep = ''),
               treatment_arm_sample_size = 1000,
               placebo_arm_sample_size = 1000,
               recruitment_start_date = 739125,
               trial_duration = 180,
               dropout_rate = 0,
               trial_PQ_eff_stratum_1 = 1,
               trial_PQ_eff_stratum_2 = 0,
               trial_PQ_prop_stratum_1 = 0.75,
               trial_PQ_prop_stratum_2 = 0.25,
               trial_PQ_lowage = 16 * 365,
               trial_PQ_proph = param_sweep[ff, 'trial_PQ_proph'],
               trial_observation_period = 3,
               is_LLIN_distributed = 1,
               is_IRS_administered = 0,
               G6PD_activity_threshold = 7,
               followup_dates = c(0,1,2,3,8,15,22,29,60,120,180))
  
  genTrialFile(output_file_path = paste(path_out, 'leaky/trial_design_', ff - 1, '.txt', sep = ''),
               treatment_arm_sample_size = 1000,
               placebo_arm_sample_size = 1000,
               recruitment_start_date = 739125,
               trial_duration = 180,
               dropout_rate = 0,
               trial_PQ_eff_stratum_1 = 0.75,
               trial_PQ_eff_stratum_2 = 0,
               trial_PQ_prop_stratum_1 = 1,
               trial_PQ_prop_stratum_2 = 0,
               trial_PQ_lowage = 16 * 365,
               trial_PQ_proph = param_sweep[ff, 'trial_PQ_proph'],
               trial_observation_period = 3,
               is_LLIN_distributed = 1,
               is_IRS_administered = 0,
               G6PD_activity_threshold = 7,
               followup_dates = c(0,1,2,3,8,15,22,29,60,120,180))
  
  # model parameters file 
  genModelParamsFile(output_file_path = paste(path_out, 'all_or_none/model_params_', ff - 1, '.txt', sep = ''),
                     N_part = 100000,
                     EIR_equil = param_sweep[ff, 'EIR_equil'],
                     sig_het = 0,
                     start_time = 2019,
                     end_time = 2040)
  
  genModelParamsFile(output_file_path = paste(path_out, 'leaky/model_params_', ff - 1, '.txt', sep = ''),
                     N_part = 100000,
                     EIR_equil = param_sweep[ff, 'EIR_equil'],
                     sig_het = 0,
                     start_time = 2019,
                     end_time = 2040)
}

