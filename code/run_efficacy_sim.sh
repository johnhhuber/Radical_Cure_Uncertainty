#!/bin/bash
#$ -q long
#$ -N eff_sim
#$ -t 1-10


for REP in {1..50}
do
	./a.out ../data/efficacy_sim/model_params_$SGE_TASK_ID.txt ../data/domestic_mosquitoes.txt ../data/occupational_mosquitoes.txt ../data/efficacy_sim/intervention_cov_$SGE_TASK_ID.txt ../output/efficacy_sim/output_${SGE_TASK_ID}_${REP}.txt
	bzip2 ../output/efficacy_sim/output_${SGE_TASK_ID}_${REP}.txt
done

