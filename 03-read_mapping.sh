#!/bin/bash
#SBATCH -A snic2019-3-67
#SBATCH -M snowy
#SBATCH -p core -n 8
#SBATCH -t 7-00:00:00
#SBATCH -J HOM_map
#SBATCH -e HOM_map_%J_%A_%a.err
#SBATCH -o HOM_map_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-12

# Load required software.
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.20.4
module load QualiMap/2.2.1

# Directory where the BAM files will be stored.
cd /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/03-bam-files

# Set path to required files and directories as environment variables.
RES_DIR='/proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/03-bam-files'
R1_FILE='./R1_files_bwa.txt'
R2_FILE='./R2_files_bwa.txt'
ID_FILE='./sample_IDs_bwa.txt'
REF_FILE='/proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/00-genome/v.1.0/fTraTra1.PB.asm1.purge1.scaff1.fa'

# For a given job in the array, set the correspondent sample files.
R1_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $R1_FILE)
R2_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $R2_FILE)
ID_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $ID_FILE)

# Print current file info to stdout for future reference.
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the R1 target: ${R1_target}
echo -e This is the R2 target: ${R2_target}
echo -e This is the ID target: ${ID_target}
echo -e This is the REF genome: ${REF_FILE}

# Create directory for temporal files.
if [ -d $SNIC_TMP/tmp ]; then echo "tmp/ exists in SNIC_TMP"; else mkdir $SNIC_TMP/tmp; fi

echo -e $(date -u) ": Read mapping began..."

# Map reads to reference genome using several threads [-t 8] and mark split alignment [-M]; then sort reads and generate the bam file index.
bwa mem -M -t 8 -R '@RG\tID:${ID_target}\tSM:${ID_target}' ${REF_FILE} ${R1_target} ${R2_target} | samtools view -@ 8 -b -S - > $SNIC_TMP/${ID_target}.bam
samtools sort -@ 8 -T $SNIC_TMP/tmp -o $SNIC_TMP/${ID_target}.sort.bam $SNIC_TMP/${ID_target}.bam && rm $SNIC_TMP/${ID_target}.bam
samtools index -@ 8 $SNIC_TMP/${ID_target}.sort.bam

echo -e $(date -u) ": Read mapping and bam file indexing ended..."

# Mark duplicate reads.
if [ -d $RES_DIR/MarkDup_metrics ]; then echo "MarkDup_metrics/ exists"; else mkdir $RES_DIR/MarkDup_metrics; fi
echo -e $(date -u) ": Mark duplicates began..."

java -Xmx60G -jar /sw/apps/bioinfo/picard/2.20.4/rackham/picard.jar MarkDuplicates \
I=$SNIC_TMP/${ID_target}.sort.bam O=$SNIC_TMP/${ID_target}.sort.MarkDup.bam M=$RES_DIR/MarkDup_metrics/${ID_target}.MarkDup.txt \
ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=$SNIC_TMP/tmp && rm $SNIC_TMP/${ID_target}.sort.bam

echo -e $(date -u) ": Mark duplicates ended..."

# Add read groups and generate bam file index.
java -Xmx60G -jar /sw/apps/bioinfo/picard/2.20.4/rackham/picard.jar AddOrReplaceReadGroups \
I=$SNIC_TMP/${ID_target}.sort.MarkDup.bam O=$SNIC_TMP/${ID_target}.sort.MarkDup.RG.bam \
RGID=${ID_target} RGLB=${ID_target} RGPL=illumina RGPU=${ID_target} RGSM=${ID_target} CREATE_INDEX=TRUE && rm $SNIC_TMP/${ID_target}.sort.MarkDup.bam

# Copy final BAM file and its index to the Results directory.
cp $SNIC_TMP/${ID_target}.sort.MarkDup.RG.bam $RES_DIR
cp $SNIC_TMP/${ID_target}.sort.MarkDup.RG.bam.bai $RES_DIR

# Obtain mapping quality summary statistics.
if [ -d $RES_DIR/Qualimap_results ]; then echo "Qualimap_results/ exists"; else mkdir $RES_DIR/Qualimap_results; fi

echo -e $(date -u) ": Qualimap began..."

qualimap --java-mem-size=60G bamqc -bam $SNIC_TMP/${ID_target}.sort.MarkDup.RG.bam -ip -sd -sdmode 0 -outfile $RES_DIR/Qualimap_results/${ID_target}

echo -e $(date -u) ": Qualimap ended..."
