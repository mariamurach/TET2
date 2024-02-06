#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --mem-per-cpu=10000
#SBATCH --time=3:00:00
#SBATCH --output=map-%a.out 
#SBATCH -p _______
#SBATCH	-A _______
#SBATCH --cpus-per-task=1
#SBATCH --array=1-6

module load star

p=20

trimmed=$curdir/data/trimmed
outdir=$curdir/data/mapped
mkdir -p $outdir
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
i=$(ls $trimmed/*.R1.fq.gz | awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'FNR == ArrayTaskID {print}')
echo $i

f=`echo $i | awk -F".R1.fq.gz" '{print $1}'`    
f1="${f}.R1.fq.gz"
f2="${f}.R2.fq.gz"
STAR --genomeDir $curdir/genomes/ \
--runMode alignReads \
--sjdbGTFfile genes.gtf \
--readFilesIn $f1 $f2 \
--runThreadN $p \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $outdir/$(basename $f) \
--quantMode GeneCounts \
--readFilesCommand pigz -dc -p $p -k