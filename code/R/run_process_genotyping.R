# load necessary functions
source('functions_genotyping.R')

# install necessary packages
if(!require(doParallel)){install.packages('doParallel'); library(doParallel)}
if(!require(foreach)){install.packages('foreach'); library(foreach)}

# specify command line arguments
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

batch_size = round(n_sims / 500)
batches = split(1:n_sims, ceiling(seq_along(1:n_sims) / batch_size))

# loop through and create the necessary files 
df_final <- foreach(ss = batches[[batch_id]], .combine = rbind, .packages = 'survival') %dopar%{
  # create temporary dataframe 
  df <- data.frame(param_id = rep(NA, 1),
                   eir_equil = rep(NA, 1),
                   sig_het = rep(NA, 1),
                   ff = rep(NA, 1),
                   trial_duration = rep(NA, 1),
                   trial_PQ_eff_stratum_1 = rep(NA, 1),
                   trial_PQ_eff_stratum_2 = rep(NA, 1),
                   trial_PQ_prop_stratum_1 = rep(NA, 1),
                   trial_PQ_prop_stratum_2 = rep(NA, 1))
  
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
  
  # load output files
  output_participants <- read.csv(file = grep(pattern = paste('_', param_id, '.csv', sep = ''), files_output_participants, value = T), header = T)
  output_recurrent <- read.csv(file = grep(pattern = paste('_', param_id, '.csv', sep = ''), files_output_recurrent, value = T), header = T)
  output_trial <- read.csv(file = grep(pattern = paste('_', param_id, '.csv', sep = ''), files_output_trial, value = T), header = T)
  
  # specify the range of genotyping sensitivities and specificities considered 
  sensitivity <- seq(from = 0, to = 1, by = 0.25)
  specificity <- seq(from = 0, to = 1, by = 0.25)
  genotyping_sens_spec <- expand.grid(sensitivity = sensitivity, specificity = specificity)
  
  # organize the trial data 
  trial_data <- createDataFrame(output_participant = output_participants,
                                   output_recurrent = output_recurrent,
                                   output_trial = output_trial,
                                   genotyping_sens_spec = genotyping_sens_spec)
  
  # get the incidence efficacy calculations 
  eff_incid <- calcIncidEff(trial_data)
  df[names(eff_incid)] <- eff_incid

  # get the cph efficacy calculations
  eff_cph <- calcCoxEff(trial_data)
  df[names(eff_cph)] <- eff_cph
  
  # return 
  df
}

# write to file
write.csv(df_final, file = paste(path_output, 'output_genotyping_', batch_id, '.csv', sep = ''), row.names = F)
