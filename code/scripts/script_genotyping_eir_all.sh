#!/bin/tcsh
#$ -q long
#$ -N eir
#$ -pe smp 32
#$ -t 1-50

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_genotyping.R data/analysis_0320/eir_vs_heterogeneity/all_or_none/ data/analysis_0320/eir_vs_heterogeneity/all_or_none/ eir_vs_heterogeneity_all_or_none/ $SGE_TASK_ID
