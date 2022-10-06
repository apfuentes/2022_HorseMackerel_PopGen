#!/bin/bash

#SBATCH -A snic2019-3-67
#SBATCH -p core -n 5
#SBATCH -t 10-00:00:00
#SBATCH -J HOM_genoytpeGVCF
#SBATCH -e HOM_genoytpeGVCF_%J_%A_%a.err
#SBATCH -o HOM_genoytpeGVCF_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-10

# Load required software.
module load bioinfo-tools
module load GATK/3.8-0

# Go to the directory where the gVCF files will be stored.
cd /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/05-VCF-files

# Set path to required files and directories as environment variables.
echo This is array job: $SLURM_ARRAY_TASK_ID

INTERVAL_FILE='./intervals_genotypeGVCF.txt'
INTERVAL_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $INTERVAL_FILE)
echo -e This is the Interval target: ${INTERVAL_target}

# Create directory for temporal files.
if [ -d $SNIC_TMP/tmp ]; then echo "tmp/ exists in SNIC_TMP"; else mkdir $SNIC_TMP/tmp; fi

# Run GATK HaplotypeCaller to generate a single gVCF file for each pool sample and scaffold interval.
java -Djava.io.tmpdir=$SNIC_TMP/tmp -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx25G -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/00-genome/v.1.0/fTraTra1.PB.asm1.purge1.scaff1.fa \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/1a_West_Ireland_2016/1a_West_Ireland_2016.${INTERVAL_target}.g.vcf \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/1b_Southwest_Ireland_2016/1b_Southwest_Ireland_2016.${INTERVAL_target}.g.vcf \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/2_West-Southwest_Ireland_2017/2_West-Southwest_Ireland_2017.${INTERVAL_target}.g.vcf \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/3_Southern_NorthSea_2016/3_Southern_NorthSea_2016.${INTERVAL_target}.g.vcf \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/4_Southern_NorthSea_2017/4_Southern_NorthSea_2017.${INTERVAL_target}.g.vcf \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/5a_Northern_Portugal_2016/5a_Northern_Portugal_2016.${INTERVAL_target}.g.vcf \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/5b_Southern_Portugal_2016/5b_Southern_Portugal_2016.${INTERVAL_target}.g.vcf \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/6a_Northern_Portugal_2017/6a_Northern_Portugal_2017.${INTERVAL_target}.g.vcf \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/6b_Southern_Portugal_2017/6b_Southern_Portugal_2017.${INTERVAL_target}.g.vcf \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/7_NorthAfrica_Mauritania_2016/7_NorthAfrica_Mauritania_2016.${INTERVAL_target}.g.vcf \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/8_Northern_SpanishShelf_2016/8_Northern_SpanishShelf_2016.${INTERVAL_target}.g.vcf \
--variant /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/04-gVCF-files/9_Mediterranean_AlboranSea_2018/9_Mediterranean_AlboranSea_2018.${INTERVAL_target}.g.vcf \
-o ./HOM_12_pools.${INTERVAL_target}.raw.vcf.gz
