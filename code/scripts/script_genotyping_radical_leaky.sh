#!/bin/tcsh
#$ -q long
#$ -N eir
#$ -pe smp 32
#$ -t 1-50

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_genotyping.R data/analysis_0320/radical_cure_therapeutic/leaky/ data/analysis_0320/radical_cure_therapeutic/leaky/ radical_cure_therapeutic_leaky/ $SGE_TASK_ID
