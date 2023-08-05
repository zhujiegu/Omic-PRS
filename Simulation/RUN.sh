#!/bin/sh

for comp in 1 6 11
do 
	for cut_pct in 1 0.1 0.01 0.001
	do 
		for hr_y in 0.2 0.5
		do
			sbatch -J $comp'_'$cut_pct'_'$hr_y simu_run.sh $comp $cut_pct $hr_y
			sbatch -J 's'$comp'_'$cut_pct'_'$hr_y simu_run_snp.sh $comp $cut_pct $hr_y
		done
	done
done

