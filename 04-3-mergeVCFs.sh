#!/bin/bash

#SBATCH -A snic2018-8-64
#SBATCH -M snowy
#SBATCH -p core -n 4
#SBATCH -t 3-00:00:00
#SBATCH -J HOM_mergeVCF
#SBATCH -e HOM_mergeVCF_%J_%A_%a.err
#SBATCH -o HOM_mergeVCF_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com

# Load required software.
module load bioinfo-tools
module load picard/2.20.4
#module load samtools/1.9

# Go to the directory where the VCF file will be stored.
cd /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/05-VCF-files

# Create directory for temporal files.
if [ -d $SNIC_TMP/tmp ]; then echo "tmp/ exists in SNIC_TMP"; else mkdir $SNIC_TMP/tmp; fi

java -Djava.io.tmpdir=$SNIC_TMP/tmp -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx20g -jar /sw/apps/bioinfo/picard/2.20.4/rackham/picard.jar \
MergeVcfs \
INPUT=12_pools_HOM_10intervals_sorted.list \
OUTPUT=HOM_12_pools.raw.variants.vcf.gz
