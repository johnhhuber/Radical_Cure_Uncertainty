#!/bin/tcsh
#$ -q long
#$ -N vec
#$ -t 1-9600

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_survival.R vector_control_all_or_none/indiv/ vector_control_all_or_none/recurr/ vector_control_all_or_none/survival/ $SGE_TASK_ID
