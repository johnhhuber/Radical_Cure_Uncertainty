#!/bin/tcsh
#$ -q long
#$ -N eir
#$ -pe smp 32
#$ -t 1-50

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_results.R data/analysis_0320/radical_cure_therapeutic/all_or_none/ data/analysis_0320/radical_cure_therapeutic/all_or_none/ radical_cure_therapeutic_all_or_none/ $SGE_TASK_ID
