#!/bin/bash

#SBATCH -A snic2019-3-67
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1-00:00:00
#SBATCH -J HOM_apply-filters
#SBATCH -e HOM_apply-filters_%J_%A_%a.err
#SBATCH -o HOM_apply-filters_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH --mail-type=ALL

# Load required programs.
module load bioinfo-tools
module load GATK/3.8-0
module load bcftools/1.10
module load zlib/1.2.11
module load samtools/1.10

# Set working directory.
cd /proj/snic2020-16-14/private/HorseMackerel/data/05-VCF-files

# Create directory for temporal files.
if [ -d $SNIC_TMP/tmp ]; then echo "tmp/ exists in SNIC_TMP"; else mkdir $SNIC_TMP/tmp; fi

# 1. Apply GATK hard filters. Flag SNPs that do not pass GATK hard-filters.
java -Djava.io.tmpdir=$SNIC_TMP/tmp -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx28G -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /proj/snic2020-16-14/private/HorseMackerel/data/00-genome/v.1.0/fTraTra1.PB.asm1.purge1.scaff1.fa \
-V HOM_12_pools.SNPs.raw.vcf.gz \
--filterExpression "FS > 60.0" --filterName "FS" \
--filterExpression "SOR > 3.0" --filterName "SOR" \
--filterExpression "MQ < 40.0" --filterName "MQ" \
--filterExpression "MQRankSum < -12.5" --filterName "MQRankSum" \
--filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum" \
-o HOM_12_pools.SNPs.hf-flag.vcf

# Retain only variants that passed these filters (source: https://evodify.com/gatk-in-non-model-organism/).
grep -E '^#|PASS' HOM_12_pools.SNPs.hf-flag.vcf > HOM_12_pools.SNPs.hf.vcf

# 2. Apply filters to the FORMAT (genotype) field.
java -Djava.io.tmpdir=$SNIC_TMP/tmp -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx28G -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /proj/snic2020-16-14/private/HorseMackerel/data/00-genome/v.1.0/fTraTra1.PB.asm1.purge1.scaff1.fa \
-V HOM_12_pools.SNPs.hf.vcf \
--genotypeFilterExpression "DP < 10 || DP > 885 || GQ < 10" \
--genotypeFilterName "filterGT" \
--setFilteredGtToNocall \
-o HOM_12_pools.SNPs.hf.DP10-885.GQ10.vcf.gz

# 3. Apply monomorphic, missingness and MAF filters.
# Remove monomorphic loci.
bcftools filter -e 'AC==0 || AC==AN' -O z -o HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.vcf.gz HOM_12_pools.SNPs.hf.DP10-885.GQ10.vcf.gz

# obtain missingness per pool sample.
bcftools stats -s - HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.vcf.gz | grep -E ^PSC | cut -f3,14 > HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.imiss

# generate a sample missingness plot.
#R --vanilla < plot_imiss_pools.R

# Filter for missing rate and MAF.
bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' -O z -o HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.miss20.maf0.05.vcf.gz HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.vcf.gz

# Generate an index for the final VCF file.
tabix /proj/snic2020-16-14/private/HorseMackerel/data/05-VCF-files/HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.miss20.maf0.05.vcf.gz