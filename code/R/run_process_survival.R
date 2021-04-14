# load necessary functions 
source('functions_survival.R')

# specify path in 
args <- commandArgs(trailingOnly = T)
path_input_indiv <- args[1]
path_input_recurr <- args[2]
path_output <- args[3]
sim_id <- as.numeric(args[4])

# decrement sim_id by 1 to account for numbering in condor
sim_id <- sim_id - 1

# load the files 
file_indiv <- paste(path_input_indiv, 'indiv_', sim_id, '.csv.bz2', sep = '')
file_recurr <- paste(path_input_recurr, 'recurr_', sim_id, '.csv.bz2', sep = '')

indiv <- read.csv(file_indiv)
recurr <- read.csv(file_recurr)

# generate survival curve 
df <- genSurvDf(recurr = recurr, indiv = indiv)

# write to file 
write.csv(df, file = paste(path_output, 'survival_', sim_id, '.csv', sep = ''), row.names = F)
