#!/bin/sh

# gly screened
#for rx in 3 6 9 12 16
#do 
#	sbatch -J 'rx_'$rx rx_sbatch.sh $rx
#done

# metb screened
for rx in 2 5 10 12
do 
	sbatch -J 'rx_'$rx rx_sbatch.sh $rx
done

