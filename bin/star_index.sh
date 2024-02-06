#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --mem-per-cpu=90000
#SBATCH --time=2:00:00
#SBATCH --output=index2.out
#SBATCH -p _______
#SBATCH	-A _______



module load star

STAR --runThreadN 30 \
--runMode genomeGenerate \
--genomeDir $curdir/genomes/temp  \
--genomeFastaFiles genome.fa \
--sjdbGTFfile genes.gtf
