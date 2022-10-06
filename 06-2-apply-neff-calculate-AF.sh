#!/bin/bash

#SBATCH -A snic2018-8-64
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1-00:00:00
#SBATCH -J HOM_neff-calc-AF
#SBATCH -e HOM_neff-calc-AF_%J_%A_%a.err
#SBATCH -o HOM_neff-calc-AF_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH --mail-type=ALL

# Load required programs.
module load bioinfo-tools
module load python3/3.7.2

# Apply Neff correction to the raw read counts per allele (AD).
python3 /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/code/utility-code/apply-Neff-correction_to_AD.py \
-i /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/analysis/06-apply-Neff-calc-AF/HOM_12_pools.SNPs.clean.AD.csv \
-o /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/analysis/06-apply-Neff-calc-AF/HOM_12_pools.SNPs.clean.AD.Neff.csv

# Calculate population-level allele frequencies for each SNP.
python3 /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/code/utility-code/calculate_allelefreq_from_AD.py \
-i /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/analysis/06-apply-Neff-calc-AF/HOM_12_pools.SNPs.clean.AD.Neff.csv \
-o /proj/uppstore2017191/nobackup/private/UserDirectories/angela/HorseMackerel/analysis/06-apply-Neff-calc-AF/HOM_12_pools.SNPs.clean.AF.Neff.txt
