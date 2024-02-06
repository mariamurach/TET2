#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=39
#SBATCH --mem-per-cpu=19000
#SBATCH --time=24:00:00
#SBATCH --output=job_output_qc_paired.out
#SBATCH -p standard
#SBATCH	-A iprime


module load fastqc 
fastqc -o fastqc/ -t 39 trimmed/*fq*
