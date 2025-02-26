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
n_sims <- length(files_output_participants)
print(n_sims)

batch_size = round(n_sims / 1000)
batches = split(1:n_sims, ceiling(seq_along(1:n_sims) / batch_size))

# loop through and create the necessary files 
df_final <- foreach(ss = batches[[batch_id]], .combine = rbind, .packages = 'survival') %dopar%{
  # create temporary dataframe 
  df <- data.frame(param_id = rep(NA, 1),
                   prop_relapse = rep(NA, 1))
  
  param_id <- as.numeric(head(strsplit(tail(strsplit(files_output_participants[ss], '_')[[1]], n = 1), '[.]')[[1]], n = 1))
  
  # load the parameter files 
  model_params <- read.table(file = grep(pattern = paste('_', param_id, '.txt', sep = ''), files_input_params, value = T), header = F)

  # load the trial design file 
  trial_design <- read.csv(file = grep(pattern = paste('_', param_id, '.txt', sep = ''), files_input_trial, value = T), header = T)
  
  # fill in the parameters 
  df$param_id <- param_id

  # specify files 
  file_participants <- grep(pattern = paste('_', param_id, '.csv', sep = ''), files_output_participants, value = T)
  file_recurrent <- grep(pattern = paste('_', param_id, '.csv', sep = ''), files_output_recurrent, value = T)

  # load output files
  if(file.exists(file_participants) & file.exists(file_recurrent))
  {
    output_participants <- read.csv(file = file_participants, header = T)
    output_recurrent <- read.csv(file = file_recurrent, header = T)

    df$prop_relapse <- sum(grepl('PRE_', output_recurrent$Infection_Type)) / nrow(output_recurrent)
  }

  # return 
  df
}

# write to file 
write.csv(df_final, file = paste(path_output, 'output_recurrent_', batch_id, '.csv', sep = ''), row.names = F)
