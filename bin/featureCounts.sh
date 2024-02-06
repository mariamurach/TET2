#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=10000
#SBATCH --time=24:00:00
#SBATCH --output=job_output_counts.out
#SBATCH -p standard
#SBATCH	-A iprime

for FILE in files/*.bam; 
do echo $FILE
out="$(basename $FILE)"
mkdir files/counts
featureCounts -T 39 -s 2 \
-a gencode.vM32.annotation.gtf \
-o files/counts/"$out"_counts.txt \
$FILE; done