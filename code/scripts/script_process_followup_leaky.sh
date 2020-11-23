#!/bin/tcsh
#$ -q long
#$ -N followup
#$ -pe smp 8
#$ -t 1-500

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_results.R data/analysis_1119/followup_vs_relapse/leaky/model_params/ data/analysis_1119/followup_vs_relapse/leaky/trial_design/ followup_vs_relapse_leaky/ $SGE_TASK_ID
