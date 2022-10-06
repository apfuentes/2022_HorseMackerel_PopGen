#!/bin/bash
#SBATCH -A snic2020-5-36
#SBATCH -p core -n 4
#SBATCH -t 4-00:00:00
#SBATCH -J HOM_Tpi
#SBATCH -e HOM_Tpi_%J_%A_%a.err
#SBATCH -o HOM_Tpi_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-12

# Load required software.
module load bioinfo-tools
module load popoolation/1.2.2
module load perl/5.26.2

# Directory where the BAM files will be stored.
cd /proj/snic2020-16-14/private/HorseMackerel/analysis/07-popgen-stats

# Set path to required files and directories as environment variables.
ID_FILE='./HOM_12pops_sampleIDs.txt'
POOLSIZE_FILE='./HOM_12pops_poolSize.txt'

# For a given job in the array, set the correspondent sample files.
ID_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $ID_FILE)
POOLSIZE_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $POOLSIZE_FILE)

# Print current file info to stdout for future reference.
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the ID target: ${ID_target}
echo -e This is the POOLSIZE target: ${POOLSIZE_target}

# Calculate Tajima’s Pi.
perl /sw/apps/bioinfo/popoolation/1.2.2/rackham/Variance-sliding.pl --input ./pileup-files/${ID_target}.pileup --output ./Tajimas-Pi/${ID_target}.pileup.WS-SS10K.minCnt.4.minCov.10.maxCov.885.pi --measure pi --window-size 10000 --step-size 10000 --min-count 4 --min-coverage 10 --max-coverage 885 --min-qual 20 --pool-size ${POOLSIZE_target} --fastq-type sanger
