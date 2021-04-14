#!/bin/tcsh
#$ -q long
#$ -N radical
#$ -t 1-1200

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_survival.R radical_cure_therapeutic_all_or_none/indiv/ radical_cure_therapeutic_all_or_none/recurr/ radical_cure_therapeutic_all_or_none/survival/ $SGE_TASK_ID
