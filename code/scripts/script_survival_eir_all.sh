#!/bin/tcsh
#$ -q long
#$ -N eir
#$ -t 1-3200

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_survival.R eir_vs_heterogeneity_all_or_none/indiv/ eir_vs_heterogeneity_all_or_none/recurr/ eir_vs_heterogeneity_all_or_none/survival/ $SGE_TASK_ID
