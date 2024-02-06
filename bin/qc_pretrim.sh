#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10000
#SBATCH --time=1:00:00
#SBATCH --output=qc_pretrim.out
#SBATCH -p largemem
#SBATCH	-A iprime


module load fastqc 

outdir=$curdir/data/fastqc/pretrim
mkdir -p $outdir
fastqc -o $outdir -t 20 $curdir/raw/*fq* 
