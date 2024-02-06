#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=100000
#SBATCH --time=2:00:00
#SBATCH --output=rscript.out
#SBATCH -p _______
#SBATCH	-A _______

module load gcc/9.2.0 openmpi/3.1.6 R/4.2.1

#i=$(ls ../mapped/*.bam| awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'FNR == ArrayTaskID {print}')
#echo $i
for i in $(ls ../mapped/*.bam)
do
    f=`echo $i | awk -F"Aligned.sortedByCoord.out.bam" '{print $1}'`    
    Rscript bcr_analysis.R --sample $(basename $f)
done
