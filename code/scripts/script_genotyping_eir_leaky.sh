#!/bin/tcsh
#$ -q long
#$ -N eir
#$ -pe smp 8
#$ -t 1-500

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_genotyping.R data/analysis_1119/eir_vs_heterogeneity/leaky/ data/analysis_1119/eir_vs_heterogeneity/leaky/ eir_vs_heterogeneity_leaky/ $SGE_TASK_ID
