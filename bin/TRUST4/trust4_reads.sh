#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=29
#SBATCH --mem=100000
#SBATCH --time=3:00:00
#SBATCH --output=trust4_%a.out
#SBATCH -p ______
#SBATCH	-A ______
#SBATCH --array=1-21

module load java

output=output_reads/
mkdir $output
i=$(ls trimmed/*.R1.fq.gz | awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'FNR == ArrayTaskID {print}')
echo $i
f=`echo $i | awk -F".R1.fq.gz" '{print $1}'`    
f1="${f}.R1.fq.gz"
f2="${f}.R2.fq.gz"

#cpu=19000
#run-trust4 -t 40 -f bcrtcr.fa --ref bcr/IMGT+C.fa -1 trimmed/A1.R1.fq.gz -2 trimmed/A2.R2.fq.gz -o $output
run-trust4 -t 40 \
-f GRCm38_bcrtcr.fa \
--ref mouse_IMGT+C.fa \
-1 $f1 -2 $f2 \
--outputReadAssignment \ 
-o $output/$(basename $f)
