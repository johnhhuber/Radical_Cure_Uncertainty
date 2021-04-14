#!/bin/tcsh
#$ -q long
#$ -N followup
#$ -pe smp 32
#$ -t 1-50

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_results.R data/analysis_1119/vector_control/leaky/model_params/ data/analysis_1119/vector_control/leaky/trial_design/ vector_control_leaky/ $SGE_TASK_ID
