#!/bin/bash
#SBATCH -A snic2019-3-67
#SBATCH -M snowy
#SBATCH -p core -n 6
#SBATCH -t 1-00:00:00
#SBATCH -J arr_trim_HOM_pool
#SBATCH -e trim_HOM_pool_%J_%A_%a.err    # File to which STDERR will be written
#SBATCH -o trim_HOM_pool_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-9

# Load libraries.
module load bioinfo-tools
module load trimmomatic/0.36
module load FastQC/0.11.8

echo This is array job: $SLURM_ARRAY_TASK_ID

R1_FILE='./R1_files.txt'
R2_FILE='./R2_files.txt'
OUT_DIR='/proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/02-clean-reads'
OUT_SUBDIR='./sample_subDir.txt'
ID_FILE='./sample_IDs.txt'
ADAPTERS_FILE='/proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/02-clean-reads/adapters.fa'

R1_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $R1_FILE)
R2_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $R2_FILE)
ID_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $ID_FILE)
SUBDIR_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $OUT_SUBDIR)

echo -e This is the R1 target: ${R1_target}
echo -e This is the R2 target: ${R2_target}
echo -e This is the ID target: ${ID_target}
echo -e This is the R1 output file name: ${OUT_DIR}/${SUBDIR_target}/${ID_target}_trimmed_R1_paired.fq.gz
echo -e This is the R2 output file name: ${OUT_DIR}/${SUBDIR_target}/${ID_target}_trimmed_R2_paired.fq.gz

# Create a sub-folder for each sample
mkdir ${OUT_DIR}/${SUBDIR_target}
mkdir ${OUT_DIR}/${SUBDIR_target}/unpaired

# Run trimmomatic on raw sequence reads - Paired-End
java -Xmx38g -jar /sw/apps/bioinfo/trimmomatic/0.36/rackham/trimmomatic-0.36.jar PE -threads 6 -phred33 \
${R1_target} ${R2_target} \
${OUT_DIR}/${SUBDIR_target}/${ID_target}_trimmed_R1_paired.fq.gz ${OUT_DIR}/${SUBDIR_target}/unpaired/${ID_target}_trimmed_R1_unpaired.fq.gz ${OUT_DIR}/${SUBDIR_target}/${ID_target}_trimmed_R2_paired.fq.gz ${OUT_DIR}/${SUBDIR_target}/unpaired/${ID_target}_trimmed_R2_unpaired.fq.gz \
ILLUMINACLIP:${ADAPTERS_FILE}:2:40:15:8:true SLIDINGWINDOW:4:15 LEADING:15 TRAILING:15 MINLEN:36

# Run fastQC on the trimmed reads file
fastqc ${OUT_DIR}/${SUBDIR_target}/${ID_target}_trimmed_R1_paired.fq.gz ${OUT_DIR}/${SUBDIR_target}/${ID_target}_trimmed_R2_paired.fq.gz -o ${OUT_DIR}/FastQC_results -t 6
