#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=1000000
#SBATCH --time=2:00:00
#SBATCH --output=mixcr_%a.out
#SBATCH -p _______
#SBATCH	-A _______
#SBATCH --array=1-21

module load java
output=bcr/output
mkdir $output
i=$(ls mapped/*.bam| awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'FNR == ArrayTaskID {print}')
echo $i
f=`echo $i | awk -F"Aligned.sortedByCoord.out.bam" '{print $1}'`    

#cpu=19000
#run-trust4 -t 40 -f bcrtcr.fa --ref bcr/IMGT+C.fa -1 trimmed/A1.R1.fq.gz -2 trimmed/A2.R2.fq.gz -o $output
run-trust4 -t 40 -f bcr/bcrtcr.fa --ref bcr/IMGT+C.fa -b $i -o $output/$(basename $f)
