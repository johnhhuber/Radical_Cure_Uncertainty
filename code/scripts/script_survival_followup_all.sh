#!/bin/tcsh
#$ -q long
#$ -N eir
#$ -t 1-9600

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_survival.R followup_vs_relapse_all_or_none/indiv/ followup_vs_relapse_all_or_none/recurr/ followup_vs_relapse_all_or_none/survival/ $SGE_TASK_ID
