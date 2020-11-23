#!/bin/tcsh
#$ -q long
#$ -N followup
#$ -pe smp 8
#$ -t 1-500

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_results.R data/analysis_1119/vector_control/all_or_none/model_params/ data/analysis_1119/vector_control/all_or_none/trial_design/ vector_control_all_or_none/ $SGE_TASK_ID
