#!/bin/tcsh
#$ -q long
#$ -N eir
#$ -pe smp 32
#$ -t 1-50

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_genotyping.R data/analysis_0320/followup_vs_relapse/leaky/ data/analysis_0320/followup_vs_relapse/leaky/ followup_vs_relapse_leaky/ $SGE_TASK_ID
