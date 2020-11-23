#!/bin/tcsh
#$ -q long
#$ -N eir
#$ -pe smp 8
#$ -t 1-500

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_results.R data/analysis_1119/radical_cure_therapeutic/leaky/ data/analysis_1119/radical_cure_therapeutic/leaky/ radical_cure_therapeutic_leaky/ $SGE_TASK_ID
