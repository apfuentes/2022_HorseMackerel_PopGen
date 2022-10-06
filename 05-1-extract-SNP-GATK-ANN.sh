#!/bin/bash

#SBATCH -A snic2018-8-64
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 2-00:00:00
#SBATCH -J HOM_extract-SNPs-GATK-ANN
#SBATCH -e HOM_extract-SNPs-GATK-ANN_%J_%A_%a.err
#SBATCH -o HOM_extract-SNPs-GATK-ANN_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH --mail-type=ALL

# Load required programs.
module load bioinfo-tools
module load GATK/3.8-0
module load vcftools zlib/1.2.11

# Directory where the gVCF files will be stored.
cd /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/05-VCF-files

# Create directory for temporal files.
if [ -d $SNIC_TMP/tmp ]; then echo "tmp/ exists in SNIC_TMP"; else mkdir $SNIC_TMP/tmp; fi

# Extract biallelic SNPs and exclude non variants positions.
java -Djava.io.tmpdir=$SNIC_TMP/tmp -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx20G -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/00-genome/v.1.0/fTraTra1.PB.asm1.purge1.scaff1.fa \
-V HOM_12_pools.raw.variants.vcf.gz \
--selectTypeToInclude SNP \
--restrictAllelesTo BIALLELIC \
--excludeNonVariants
-o HOM_12_pools.raw.SNPs.vcf.gz

# Obtain SNP annotations from the VCF file to make density plots and establish cutoff values for QC filtering.
java -Djava.io.tmpdir=$SNIC_TMP/tmp -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx20G -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/00-genome/v.1.0/fTraTra1.PB.asm1.purge1.scaff1.fa \
-V HOM_12_pools.raw.SNPs.vcf.gz \
-F ID -F CHROM -F POS -F FILTER -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum \
--showFiltered \
-o HOM_12_pools.raw.SNPs.ANN.table
