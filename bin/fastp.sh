#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10000
#SBATCH --time=2:00:00
#SBATCH --output=trim_%A_%a..out
#SBATCH -p _______
#SBATCH	-A _______
#SBATCH --cpus-per-task=1
#SBATCH --array=1-6


p=20
outdir=$curdir/data/trimmed
mkdir -p $outdir
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
i=$(ls $raw/*_1.fq | awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'FNR == ArrayTaskID {print}')
echo $i
f=`echo $i | awk -F"_1" '{print $1}'`    
f1="${f}_1.fq"
f2="${f}_2.fq"


fastp -i $f1 -I $f2 -w 20  -o $outdir/$(basename $f).R1.fq.gz -O $outdir/$(basename $f).R2.fq.gz
