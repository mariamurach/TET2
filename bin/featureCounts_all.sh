#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10000
#SBATCH --time=1:00:00
#SBATCH --output=count-%a.out 
#SBATCH -p _______
#SBATCH	-A _______
#SBATCH --cpus-per-task=1
#SBATCH --array=1-6

p=10
mapped=$curdir/data/mapped
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
i=$(ls $mapped/*.bam| awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'FNR == ArrayTaskID {print}')
echo $i
f=`echo $i | awk -F"Aligned.sortedByCoord.out.bam" '{print $1}'`    
counts=$curdir/data/counts
mkdir -p $counts
featureCounts -T $p -s 0 -p \
-a genes.gtf \
-o $counts/$(basename $f)_counts.txt \
$i
