#!/bin/tcsh
#$ -q long
#$ -N recurr
#$ -pe smp 8
#$ -t 1-1000

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_recurrent.R data/analysis_0320/vector_control/all_or_none/ data/analysis_0320/vector_control/all_or_none/ vector_control_all_or_none/ $SGE_TASK_ID
