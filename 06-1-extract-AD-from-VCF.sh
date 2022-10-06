#!/bin/bash

#SBATCH -A snic2018-8-64
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 2-00:00:00
#SBATCH -J HOM_extract-AD
#SBATCH -e HOM_extract-AD_%J_%A_%a.err
#SBATCH -o HOM_extract-AD_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH --mail-type=ALL

# Load required programs.
module load bioinfo-tools
module load GATK/3.8-0

# Create directory for temporal files.
if [ -d $SNIC_TMP/tmp ]; then echo "tmp/ exists in SNIC_TMP"; else mkdir $SNIC_TMP/tmp; fi

java -Djava.io.tmpdir=$SNIC_TMP/tmp -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx20G -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/00-genome/v.1.0/fTraTra1.PB.asm1.purge1.scaff1.fa \
-V /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/05-VCF-files/HOM_12_pools.SNPs.clean.recode.vcf \
-F CHROM -F POS -F REF -F ALT --genotypeFields AD \
-o /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/analysis/06-calculate-AF-apply-Neff/HOM_12_pools.SNPs.clean.AD.csv
