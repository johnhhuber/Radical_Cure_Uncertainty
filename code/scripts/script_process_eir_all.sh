#!/bin/tcsh
#$ -q long
#$ -N eir
#$ -pe smp 8
#$ -t 1-500

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_results.R data/analysis_1119/eir_vs_heterogeneity/all_or_none/ data/analysis_1119/eir_vs_heterogeneity/all_or_none/ eir_vs_heterogeneity_all_or_none/ $SGE_TASK_ID
