#####################################################################
#                                                                   #
# Filter SNPs by per sample DP interval (Mean +- 1SD)               #
# This new version enables exclusion of samples                     #
#                                                                   #
# Created by: Angela Fuentes-Pardo, inspired by code from           #
# Mats Petterson and Dmytro Kryvokhyzha                             #
# E-mail: apfuentesp@gmail.com, Uppsala University                  #
# Date: 2020-02-27, update: 2020-05-18                              #
#                                                                   #
#####################################################################

## Clean environment space.
rm(list=ls())

## Receive arguments from the bash command line.
#args <- commandArgs(trailingOnly=TRUE)
#cat(args, sep='\n')

## Assign these arguments to R objects.
#working_dir <- args[1]  # /proj/snic2020-16-14/private/HorseMackerel/data/05-VCF-files
#input_file <- args[2]   # /proj/snic2020-16-14/private/HorseMackerel/data/05-VCF-files/HOM_12_pools.SNPs.hf.vcf.DP.table
#outputName <- args[3]  # minus6a-6b-7
#samples_to_remove <- strsplit(args[4],',')[[1]]  # 6a_Northern_Portugal_2017,6b_Southern_Portugal_2017,7_NorthAfrica_Mauritania_2016
#type <- args[5]  # mode or mean

## For debugging.
working_dir <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/05-snp-filtering-UG'
#input_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/05-snp-filtering/per_pool_DP_10lines.txt'
input_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/05-snp-filtering-UG/HOM_12_pools.UG.SNPs.hf.DP20-300.GQ20.mono.miss20.mac3.DP.table'
library(tools)
filename <- file_path_sans_ext(basename(input_file))
samples_to_remove <- strsplit('6a_Northern_Portugal_2017,6b_Southern_Portugal_2017,7_NorthAfrica_Mauritania_2016', split=',')[[1]]
#samples_to_remove <- strsplit('7_NorthAfrica_Mauritania_2016', split=',')[[1]]
#samples_to_remove <- strsplit('NA', split=',')[[1]]
outputName <- 'minus6a-6b-7'
#outputName <- 'exc-6a-6b-7'
#outputName <- 'exc-7'
#outputName <- 'exc-none'
type <- 'mode'
#type <- 'mean'

print('# ------------- Set R objects from BASH command line: ------------- #')
working_dir
input_file
filename
outputName
samples_to_remove
type

## Set working directory.
setwd(working_dir)

## ------------------------------------------------------------------
## Data loading and preprocessing 
## ------------------------------------------------------------------

## Load file with DP per locus and sample (extracted with GATK-VariantsToTable function).
path_to_file <- input_file

library(data.table)
DPdata <- fread(path_to_file, data.table=FALSE, header=TRUE, sep='\t', stringsAsFactors=FALSE)

head(DPdata)  # Explore first lines of the dataframe.
DPdata$SNP <- paste0(DPdata$CHROM,'-',DPdata$POS)  # Concatenate CHROM and POS, connected with a dash, to obtain a column for the SNP position.

## Keep columns that contain '.DP' and assign SNP column to rownames.
library(dplyr)

DPdata_subset <- DPdata %>% select(matches('.DP'))
rownames(DPdata_subset) <- DPdata$SNP

## Tidy column names.
colnames(DPdata_subset)
tidy_colnames <- sub('.DP*', '', colnames(DPdata_subset))  # Remove '.DP' from the end of the sample names.
colnames(DPdata_subset) <- tidy_colnames
colnames(DPdata_subset)

## Remove samples
if(!is.null(samples_to_remove)) {  # Check if object is not NA.
  to_exclude <- samples_to_remove
  DPdata_subset <- DPdata_subset[ ,!(names(DPdata_subset) %in% to_exclude)]
  print(paste0('# ------------ Colnames after removing: ',paste0(to_exclude,collapse=', '),' are: ------------ #'))
  print(colnames(DPdata_subset))
}

## Explore the resulting df.
dim(DPdata_subset)
str(DPdata_subset)
class(DPdata_subset)

## Set NA (missing data) to 0.
DPdata_subset[is.na(DPdata_subset)] <- as.integer(0)
str(DPdata_subset)

rm(DPdata)  # Free up space

## ------------------------------------------------------------------
## Generate a DP density plot with central tendency cutoff values
## ------------------------------------------------------------------

## (Inspired by: https://evodify.com/gatk-in-non-model-organism/)

## Create a function to calculate the mode.
get.mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## Define the file names list (9 HOM samples here)
nameList <- colnames(DPdata_subset)

cutoffs_list <- matrix(nrow = ncol(DPdata_subset), ncol = 12) # define number of samples (ncol = number of cutoff values)
cutoffs_list <- data.frame(cutoffs_list, row.names = nameList)
colnames(cutoffs_list) <- c('5%', '10%', '99%', 'minCov20', 'Mean', 'SD', 'MeanMinus1SD', 'MeanPlus1SD', 
                     'Mean3x','Mode', 'ModeMinusHalf', 'ModePlusHalf')

## Initiallize pdf graphic device.
pdf(paste0('DPplot-cutoff-values',filename,'-',outputName,'.pdf'))
#jpeg(paste0('DP-distribution-and-cutoffs_',outputName,'.jpg'), height=1600, width=1200)

## Define plot layout.
nrow <- 3
ncol <- ceiling(ncol(DPdata_subset)/3)
par(mar=c(5, 3, 3, 2), cex=0.8, mfrow=c(nrow,ncol))

## Loop over samples, calculate central tendency values, and plot a DP density distribution.
for (i in 1:ncol(DPdata_subset)) {
  #i=1
  DP <- DPdata_subset[ ,c(nameList[i])]
  cutoffs_list[i,c(1,2,3)] <- quantile(DP, c(0.05, 0.1, 0.99), na.rm=T)  # 5% and 95% quantiles.
  cutoffs_list[i,c(4)] <- 20  # Min coverage (Fan).
  cutoffs_list[i,c(5)] <- round(mean(DP, na.rm=T), 1)  # Mean.
  cutoffs_list[i,c(6)] <- round(sd(DP, na.rm=T), 1)  # Standard deviation.
  cutoffs_list[i,c(7)] <- round(cutoffs_list[i,'Mean'] - cutoffs_list[i,'SD'], 1)  # Mean - 1SD.
  cutoffs_list[i,c(8)] <- round(cutoffs_list[i,'Mean'] + cutoffs_list[i,'SD'], 1)  # Mean + 1SD.
  cutoffs_list[i,c(9)] <- round(cutoffs_list[i,'Mean'] * 3)  # Mean * 3 (Carneiro).
  cutoffs_list[i,c(10)] <- round(get.mode(DP), 1)  # Mode.
  cutoffs_list[i,c(11)] <- round(cutoffs_list[i,'Mode'] - 1/2*cutoffs_list[i,'Mode'])  # Mode - 1/2 mode.
  cutoffs_list[i,c(12)] <- round(cutoffs_list[i,'Mode'] + 1/2*cutoffs_list[i,'Mode'])  # Mode + 1/2 mode.
  
  d <- density(DP, from=0, to=300, bw=1, na.rm =T)  # Generate the density distribution of DP.
  
  # Create the plot.
  plot(d, xlim=c(0,300), ylim=c(0,0.03), main=nameList[i], cex.main=1, col='blue4', lwd=1, xaxt='n',
       xlab=paste0('Depth of coverage, N = ', length(DP)), ylab='Fraction of genome')
  polygon(d, col='blue4', border='blue4')
  axis(1, at=seq(0, 300, by=10), las=2)
  abline(v=20, col='green3', lty=2, lwd=1.2)
  abline(v=300, col='green3', lty=2, lwd=1.2)  
  abline(v=cutoffs_list[i,c('Mean')], col='gray60', lty=3, lwd=1.2)
  abline(v=cutoffs_list[i,c('MeanMinus1SD')], col='magenta', lty=3, lwd=1.2)
  abline(v=cutoffs_list[i,c('MeanPlus1SD')], col='magenta', lty=3, lwd=1.2)
  abline(v=cutoffs_list[i,c('Mode')], col='gray40', lty=1, lwd=1.2)
  abline(v=cutoffs_list[i,c('ModeMinusHalf')], col='red', lty=1, lwd=1.2)
  abline(v=cutoffs_list[i,c('ModePlusHalf')], col='red', lty=1, lwd=1.2)
  abline(v=cutoffs_list[i,c(1,3)], col=' gold3', lty=4, lwd=1.2)
  legend('topright', legend=c('MinCov20, MaxCov300', 'Mean', 'Mean +- 1SD', 'Mode', 'Mode +- 1/2Mode', 'Quantile 5%, 95%'), 
                     col=c('green3', 'gray60', 'magenta', 'gray40', 'red', ' gold3'), lty=c(2,3,3,1,1,4), 
         lwd=1.2, cex=0.5, bg = 'white')
}
dev.off()

# Save table with cutoff values.
write.table(cutoffs_list, paste0('DP-central-tendency_',filename,'-',outputName,'.txt'), 
            quote = FALSE, col.names = NA, row.names = TRUE)


## ------------------------------------------------------------------
## Find loci that fulfil the DP requirements for all samples
## ------------------------------------------------------------------

if(type == 'mean'){
  ## Calculate mean depth of coverage per sample.
  Mean <- round(colMeans(DPdata_subset, na.rm = FALSE, dims = 1), 1)  # It returns a vector of the mean of the columns.
  Mean
  
  ## Calculate standard deviation (SD) of depth of coverage per sample. 
  library(matrixStats)
  set.seed(42)
  SD <- round(colSds(as.matrix(DPdata_subset), na.rm = FALSE), 1)  # It returns a vector of the standard deviation of the columns.
  SD
  
  ## Calculate DP interval (Mean +- 1SD) per sample.
  MeanMinus1SD <- Mean-SD
  MeanMinus1SD
  MeanPlus1SD <- Mean+SD
  MeanPlus1SD
  
  ## Generate a unified df with all the info (mean, SD, -1SD, +1SD).
  DPintervals <- rbind(Mean,SD,MeanMinus1SD,MeanPlus1SD)
  DPintervals
  class(DPintervals)
  
  # Save this df as a file.   HOM_12_pools.SNPs.hf.vcf
  write.table(DPintervals, file=paste0('Intervals_DPfilter_Mean1SD.',filename,'.',outputName,'.txt'), quote=FALSE, row.names=TRUE, col.names=NA)
}

if(type == 'mode'){
  
  ## Initiallize dataframe that will have the DP intervals.
  DPintervals <- data.frame(matrix(NA, nrow = 4, ncol = length(colnames(DPdata_subset))))
  colnames(DPintervals) <- colnames(DPdata_subset)
  rownames(DPintervals) <- c('Mode','Mode1half','ModeMinus1half','ModePlus1half')
  DPintervals
  
  ## Calculate the mode of the depth of coverage per sample.
  for (i in 1:ncol(DPdata_subset)) {
    #i=1
    DPintervals[c(1),i] <- round(get.mode(DPdata_subset[,i]), 1)  # Mode
    DPintervals[c(2),i] <- round(1/2*DPintervals['Mode',i], 1)  # Mode
    DPintervals[c(3),i] <- round(DPintervals['Mode',i] - DPintervals['Mode1half',i],1)  # Mode - 1/2 mode
    DPintervals[c(4),i] <- round(DPintervals['Mode',i] + DPintervals['Mode1half',i],1)  # Mode + 1/2 mode
  }
  DPintervals
  
  # Save this df as a file.   HOM_12_pools.SNPs.hf.vcf
  write.table(DPintervals, file=paste0('Intervals_DPfilter_Mode1half.',filename,'.',outputName,'.txt'), quote=FALSE, row.names=TRUE, col.names=NA)
  
}

## Assign TRUE to a given locus/sample when its depth is within the per samlpe DP interval, otherwise assign FALSE.
if(all.equal(colnames(DPintervals), colnames(DPdata_subset))){  # If samples are in the same order between the two dataframes (DP data and DP intervals).
  print(head(DPdata_subset))
  
  for(i in colnames(DPdata_subset)) {  # For loop that goes column by column (sample names) in both dataframes (DP data and intervals).
    #i='1a_West_Ireland_2016'
    #print(paste0('Interval: [ ', DPintervals[3,i],' - ', DPintervals[4,i], ' ]'))  # Print to screen the current interval under review
    #print(head(DPdata_subset[,i]))
    DPdata_subset[ ,paste0('filter_',i)] <- DPdata_subset[,i] >= DPintervals[3,i] & DPdata_subset[,i] <= DPintervals[4,i]  # If DP is within the interval assign TRUE, otherwise assign FALSE.
    #print(head(DPdata_subset))
  }
  
  print(head(DPdata_subset))
  
  ## Generate a list of SNP names for which their depth falls within the set DP intervals for all samples.
  SNPs_eval <- rowSums(DPdata_subset[ ,grep('filter_',colnames(DPdata_subset))]) == length(grep('filter',colnames(DPdata_subset)))  # Vector with TRUE/FALSE values for each rowname (locus).
  filtered_SNPs <- rownames(DPdata_subset[SNPs_eval, ])
  print(head(filtered_SNPs))
  print(paste0('# ----- Number of SNPs that passed the per-sample DPinterval filter, type (',type,'): [ ',length(filtered_SNPs),' ] ----- #'))
  
  ## Save list of SNPs to a file (each locus is in a new line).
  if(type == 'mean') {
    write.table(filtered_SNPs, file=paste0('List_SNPs_passed_DPfilter.Mean1SD.',filename,'.',outputName,'.txt'), quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
  if(type == 'mode') {
    write.table(filtered_SNPs, file=paste0('List_SNPs_passed_DPfilter.Mode1halfMode.',filename,'.',outputName,'.txt'), quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
} 
