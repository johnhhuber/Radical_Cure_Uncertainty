#!/bin/tcsh
#$ -q long
#$ -N followup
#$ -pe smp 32
#$ -t 1-50

setenv R_LIBS /afs/crc.nd.edu/user/j/jhuber3/myRlibs

module load R

Rscript run_process_genotyping.R data/analysis_0320/followup_vs_relapse/all_or_none/ data/analysis_0320/followup_vs_relapse/all_or_none/ followup_vs_relapse_all_or_none/ $SGE_TASK_ID
