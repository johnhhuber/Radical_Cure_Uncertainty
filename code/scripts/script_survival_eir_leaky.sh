#!/bin/tcsh
#$ -q long
#$ -N eir
#$ -t 1-2400

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_survival.R eir_vs_heterogeneity_leaky/indiv/ eir_vs_heterogeneity_leaky/recurr/ eir_vs_heterogeneity_leaky/survival/ $SGE_TASK_ID
