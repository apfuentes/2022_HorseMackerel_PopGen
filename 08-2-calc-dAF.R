#####################################################################
#                                                                   #
# Calculate the absolute delta allele frequency between pairs of    #
# super pools and generate a Manhattan plot                         #
#                                                                   #
# E-mail: apfuentesp@gmail.com, Uppsala University                  #
# Date: 2020-03-25                                                  #
#                                                                   #
#####################################################################

## Clean environment space.
rm(list=ls())

## ------------------------------------------------------------------
## Function definition
## ------------------------------------------------------------------

##################### Function to calculate dAF #####################
# Returns a dataframe with SNPs and dAF for a given paired comparison.

calculate.dAF <- function(poolData, comparison, group1, group2){
  ## For debugging.
  #comparison <- 'WIreland_1b_vs_WIreland_1a.2'
  #group1 <- c('X1b_Southwest_Ireland_2016')
  #group2 <- c('X1a_West_Ireland_2016', 'X2_West.Southwest_Ireland_2017')
  #comparison
  #group1
  #group2
  
  ## Subset dataframe to keep only the columns (samples) corresponding to the two superpools, and remove NAs.
  poolData_subset <- na.omit(poolData[ ,c(group1,group2)])
  head(poolData_subset)
  
  ## Calculate the row mean for Group 1.
  meanAF_group1 <- rowSums(as.data.frame(poolData_subset[,c(group1)]))/length(group1)
  head(meanAF_group1)
  length(meanAF_group1)
  
  ## Calculate the row mean for Group 2.
  meanAF_group2 <- rowSums(as.data.frame(poolData_subset[,c(group2)]))/length(group2)
  head(meanAF_group2)
  length(meanAF_group2)
  
  ## Calculate the absolute delta allele frequency (dAF) between the two superpools.
  dAF_df = data.frame(SNP = rownames(poolData_subset), 
                      dAF = round(abs(meanAF_group1 - meanAF_group2),4), stringsAsFactors=FALSE)
  head(dAF_df)
  str(dAF_df)
  
  ## Clean up some memory space.
  rm(meanAF_group1)
  rm(meanAF_group2)
  
  return(dAF_df)
  
}
########################## End of function ##########################


## ------------------------------------------------------------------
## Main body of the code 
## ------------------------------------------------------------------

## Clean environment space.
#rm(list=ls())

## Receive arguments from the bash command.
args <- commandArgs(trailingOnly=TRUE)
cat(args, sep='\n')

## Assign these arguments to variables.
working_dir <- args[1]  # /proj/snic2020-16-14/private/HorseMackerel/analysis/07-popgen-stats/pairwiseFST
input_file <- args[2]   # /proj/snic2020-16-14/private/HorseMackerel/analysis/06-apply-Neff-calc-AF/HOM_12_pools.SNPs.hf.strictDP.mono.AD.Neff.AF.txt
superpools_file <- args[3]  # /proj/snic2020-16-14/private/HorseMackerel/analysis/08-outlier-loci-detection/dAF/List_paired_comparisons_HOM_2020-02-21.csv

## For debugging.
working_dir <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF'
#input_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/06-apply-Neff-calc-AF/HOM_12_pools.SNPs.hf.strictDP.exc-6a-6b-7.mono.AD.Neff.AF.txt'
input_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF/2020-03-25/HOM_12_pools.SNPs.hf.strictDP.exc-6a-6b-7.mono.AD.Neff.AF.txt'
#superpools_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF/List_paired_comparisons_HOM_2020-04-28.csv'
superpools_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF/List_paired_comparisons_HOM_2020-05-04.csv'

print('# ------------- Set R objects from BASH command line: ------------- #')
working_dir
input_file
superpools_file

## Set working directory.
setwd(working_dir)

## ------------------------------------------------------------------
## Data loading and preprocessing 
## ------------------------------------------------------------------

## Load the pool allele frequency data.
path_to_file <- input_file

## Save filename.
library(tools)
filename <- file_path_sans_ext(basename(path_to_file))
filename

library(data.table)
poolData <- fread(path_to_file, data.table=FALSE, header=TRUE, sep='\t', stringsAsFactors=FALSE)

head(poolData)  # Explore first lines of the dataframe.
row.names(poolData) <- poolData[, 1]  # Assign first column with pool names to rownames.
#poolData <- poolData[, -1]  # Remove such first column.

head(poolData)  # Verify modification.
dim(poolData)
str(poolData)
class(poolData)

## Fix column names.
colnames(poolData)
tidy_colnames <- make.names(colnames(poolData), unique=TRUE)
colnames(poolData) <- tidy_colnames
colnames(poolData)
#[1] "CHROM.POS"                        "X1a_West_Ireland_2016"            "X1b_Southwest_Ireland_2016"      
#[4] "X2_West.Southwest_Ireland_2017"   "X3_Southern_NorthSea_2016"        "X4_Southern_NorthSea_2017"       
#[7] "X5a_Northern_Portugal_2016"       "X5b_Southern_Portugal_2016"       "X6a_Northern_Portugal_2017"      
#[10] "X6b_Southern_Portugal_2017"       "X7_NorthAfrica_Mauritania_2016"   "X8_Northern_SpanishShelf_2016"   
#[13] "X9_Mediterranean_AlboranSea_2018"

## Sort poolData df in ascending CHR and BP order.
#poolData <- poolData[order(poolData$CHROM.POS),]  
#head(poolData)

## Read in file listing paired comparisons (super pools).
comparisonPairs_df <- read.csv(superpools_file, sep=';', stringsAsFactors = FALSE)
print('# ---------------- File listing the paired comparisons ---------------- #')
comparisonPairs_df

## ------------------------------------------------------------------
## Calculate dAF for each paired comparison
## ------------------------------------------------------------------

## Create a dataframe that will store the dAF for each pairwise comparison.
summary_dAF_df = data.frame(SNP = rownames(poolData), stringsAsFactors=FALSE)
head(summary_dAF_df)

## Split loci names into CHROM and BP.
library(tidyr)
summary_dAF_df <- separate(data=summary_dAF_df, col=SNP, into=c("CHROM", "BP"), sep="\\-", remove=FALSE)
head(summary_dAF_df)

## Replace scaffold names (characters) by numbers.
summary_dAF_df$CHR <- summary_dAF_df$CHROM  # Copy the CHROM column to a new column called CHR
#summary_dAF_df <- summary_dAF_df[order(summary_dAF_df$CHROM),]
head(summary_dAF_df)

## Remove end part of the long scaffold names.
summary_dAF_df$CHR <- gsub('_arrow_ctg1', '', summary_dAF_df$CHR)
head(summary_dAF_df)
#sort(unique(summary_dAF_df$CHR))

## As there is a chromosome and scaffold with the same number - 20, I will replace scaffold_20 for scaffold_30
summary_dAF_df$CHR <- gsub('scaffold_20', 'scaffold_30', summary_dAF_df$CHR)
head(summary_dAF_df)
#sort(unique(summary_dAF_df$CHR))

## Replace the original scaffold name by a unique number (chromosome for plotting).
## For the HOM samples (when not all the scaffolds are represented and in ascending order).
library(dplyr)
library(tidyr)
num <- summary_dAF_df %>%
  separate(CHR, c("Contig", "scaffoldNum"), "_")
head(num)

## Assign the number of the chromosome (SUPER_) or scaffold to the CHR column for plotting.
summary_dAF_df$CHR <- num$scaffoldNum
head(summary_dAF_df)
#sort(as.numeric(unique(summary_dAF_df$CHR)))
#summary_dAF_df[summary_dAF_df$CHROM == 'scaffold_20_arrow_ctg1', ]
#summary_dAF_df[summary_dAF_df$CHROM == 'SUPER_20', ]

## When all the chromosomes are represented and in ascending order.
#uq <- unique(summary_dAF_df$CHR) # Obtain a list of unique scaffold names.
#num <- seq(1, length(uq)) # Create a vector with a sequence of numbers in ascending order from 1 to the number of unique scaffolds.
#for (i in seq_along(num)) summary_dAF_df$CHR[summary_dAF_df$CHR == uq[i]] <- num[i]

## Save in a TXT file the scaffold - chromosome number equivalence.
#write.table(cbind.data.frame(uq, num), file = paste0('Scaffold-CHR_equivalence_for_dAFplot_',paste(group1, collapse = '-'),'_', paste(group2, collapse = '-'), '.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

## Convert columns CHR and BP from character to numeric values.
str(summary_dAF_df)
summary_dAF_df$CHR <- as.numeric(summary_dAF_df$CHR)
summary_dAF_df$BP <- as.numeric(summary_dAF_df$BP)
str(summary_dAF_df)

## For loop that goes row by row in the df to compute dAF for each pairwise comparison.
for (i in 1:length(rownames(comparisonPairs_df))){
  ## Set group pairs for comparison.
  #i=2
  #i=1
  comparison <- strsplit(comparisonPairs_df[i,1], "\\s+")[[1]]  # strsplit() split character string into character vector by a separator, space(s) in this case
  group1 <- strsplit(comparisonPairs_df[i,2], "\\s+")[[1]]
  group2 <- strsplit(comparisonPairs_df[i,3], "\\s+")[[1]]
  
  ## Print current comparison and groups.
  print(comparison)
  print(group1)
  print(group2)
  
  ## Calculate dAF for the current pair group.
  dAF_tmp <- calculate.dAF(poolData=poolData, 
                           comparison=comparison, 
                           group1=group1, 
                           group2=group2)
  print(head(dAF_tmp))
  
  ## Add the dAF of the current pair comparison to the growing summary table of dAF. 
  summary_dAF_df[ ,paste0('dAF_',comparison)] <- dAF_tmp$dAF[match(summary_dAF_df$SNP, dAF_tmp$SNP, nomatch=NA)] # First larger df, second smaller df.
  head(summary_dAF_df)
  summary(summary_dAF_df[ ,paste0('dAF_',comparison)])
  #summary_dAF_df[is.na(summary_dAF_df$dAF_WIreland_1b_vs_WIreland_1a.2), ]  # Explore which rows are NA = no match.
  #dAF_tmp[dAF_tmp$SNP == 'tarseq_5-3819508', ]  # abscent
  #summary_dAF_df[summary_dAF_df$SNP == 'tarseq_5-3819508', ] # present
  #Example when loci match.
  #dAF_tmp[dAF_tmp$SNP == 'tarseq_85-6299535', ]  # present
  #summary_dAF_df[summary_dAF_df$SNP == 'tarseq_85-6299535', ]  # present
  
  ## Generate Manhattan plot.
  #make.Manhattanplot(summary_dAF_df=summary_dAF_df, 
  #                   filename=filename, 
  #                   comparison=comparison)
  
}

print(head(summary_dAF_df))

## Save the dataframe with all dAF values in a TXT file.
library(tools)
superpools_filename <- file_path_sans_ext(basename(superpools_file))
filename
superpools_filename

write.table(summary_dAF_df, file=paste0('Summary_dAF.',superpools_filename,'.',filename,'.txt'), 
            row.names=FALSE, sep='\t', quote=FALSE)
