#!/bin/bash
#SBATCH -A snic2018-8-64
#SBATCH -M snowy
#SBATCH -p core -n 8
#SBATCH -t 7:00:00
#SBATCH -J prep_genome
#SBATCH -e prep_genome_%J_%A_%a.err
#SBATCH -o prep_genome_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com

# Load required software.
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.20.4

cd /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/00-genome/v.1.0/

REF_FILE='/proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/data/00-genome/v.1.0/fTraTra1.PB.asm1.purge1.scaff1.fa'

# Generate BWA file index.
# where -a bwtsw specifies that we want to use the indexing algorithm that is capable of handling the whole human genome.
# Expected Result: This creates a collection of files used by BWA to perform the alignment.
bwa index -a bwtsw ${REF_FILE}

echo "################# bwa index done" ;

# Generate fasta file index.
# This creates a file called reference.fa.fai, with one record per line for each of the contigs in the FASTA reference file. Each record is composed of the contig name, size, location, basesPerLine and bytesPerLine.
samtools faidx ${REF_FILE}

echo "################# samtools faidx done" ;

# Generate sequence dictionary.
# Note that this is the new syntax for use with the latest version of Picard. Older versions used a slightly different syntax because all the tools were in separate jars, so you'd call e.g. java -jar CreateSequenceDictionary.jar directly.
# This creates a file called reference.dict formatted like a SAM header, describing the contents of your reference FASTA file.
java -Xmx60g -jar /sw/apps/bioinfo/picard/2.20.4/rackham/picard.jar CreateSequenceDictionary \
REFERENCE=${REF_FILE} \
OUTPUT=${REF_FILE}.dict

echo "################# picard dictionary done" ;
