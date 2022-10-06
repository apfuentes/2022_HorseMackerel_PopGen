#####################################################################
#                                                                   #
# Calculate the absolute delta allele frequency between pairs of    #
# super pools and generate a Manhattan plot                         #
#                                                                   #
# E-mail: apfuentesp@gmail.com, Uppsala University                  #
# Date: 2020-02-27                                                  #
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


############### Function to generate a Manhattan plot ###############

make.Manhattanplot <- function(summary_dAF_df, filename, comparison){
  ## For debugging.
  #comparison <- 'WIreland_1b_vs_WIreland_1a.2'
  
  ## Subset dataframe for the current pairwise comparison.
  dAFpair_df <- summary_dAF_df[ ,c(1,2,3,4,grep(paste0('dAF_',comparison),colnames(summary_dAF_df)))]
  head(dAFpair_df)
  
  # Sort data frame in ascending CHR and BP order.
  #dAFpair_df <- dAFpair_df[order(dAFpair_df$CHR,dAFpair_df$BP),]  
  #head(dAFpair_df)
  
  ## Open image device.
  png(paste0('Manhattan_plot_dAF_',filename,'_',comparison,'.png'), 
      #png(paste0('/proj/snic2020-16-14/private/HorseMackerel/analysis/08-outlier-loci-detection/dAF/',dir_file,'/Manhattan_plot_dAF.',dir_file,'.',comparison,'.png'), 
      #width=1900, height=1000
      height = 12, width = 24, units='cm', res=150
      #width=950, height=500
      #width=10, height=5, units='in', res=150
  )
  
  ## ------------------------------------------------------------------
  ## Make Manhattan plot using qqman 
  ## ------------------------------------------------------------------
  
  ## Load library.
  library(qqman)
  
  #manhattan(na.omit(dAFpair_df), chr='CHR', bp='BP', snp='SNP', p=comparison, 
  #          logp = FALSE, suggestiveline = FALSE, genomewideline = FALSE, 
  #          xlab = 'Scaffold', ylab = 'dAF', 
  #          main = paste0('Manhattan plot of dAF \n ', comparison)
  #main = paste0('Manhattan plot of dAF \n ', paste(group1, collapse = '-'),'\n vs. \n', paste(group2, collapse = '-')), cex.main=0.8)
  #          ,cex = 0.3)
  
  # --------------- Generate a list of outlier loci to highlight --------------- #
  ## Load library.
  library(tidyverse)
  
  ## Filter loci that fall within the top 0.995 quantile.
  #outliers_0995 <- filter(dAFpair_df, dAFpair_df[,comparison] >= quantile(dAFpair_df[,comparison], 0.995, na.rm=TRUE))
  #dim(outliers_0995)
  #head(outliers_0995)
  
  ## Filter loci that fall within the top 0.99 quantile.
  #outliers_0999 <- filter(dAFpair_df, dAFpair_df[,comparison] >= quantile(dAFpair_df[,comparison], 0.999, na.rm=TRUE))
  #dim(outliers_0999)
  #head(outliers_0999)
  
  
  ## ------------------------------------------------------------------
  ## Make Manhattan plot using ggplot 
  ## ------------------------------------------------------------------
  
  #Source: https://www.r-graph-gallery.com/101_Manhattan_plot.html
  
  ## Load library.
  library(tidyverse)
  
  ## Compute the cumulative position of SNP.
  don <- dAFpair_df %>% 
    
    # compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # add this info to the initial dataset
    left_join(dAFpair_df, ., by=c("CHR"="CHR")) %>%
    
    # add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate(BPcum=BP+tot) #%>%
  
  # ------- Assign highlight and annotation information to correspondent samples ------- #
  #mutate(is_highlight0995=ifelse(SNP %in% outliers_0995$SNP, "yes", "no")) %>%
  #mutate(is_highlight0999=ifelse(SNP %in% outliers_0999$SNP, "yes", "no")) %>%
  #mutate(is_annotate=ifelse(dAFpair_df[,comparison]>0.9, "yes", "no")) 
  
  head(don)
  str(don)
  sum(is.na(dAFpair_df[,5]))
  sum(is.na(don[,5]))
  
  ## Prepare the X axis. 
  # We do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
  axisdf = don %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum))/2)
  
  ## Calculate rollmean and add it to the dataframe.
  ## Load library.
  library(zoo)
  
  ## Set NA (missing data) to 0 as the rollmean cannot handle well NAs.
  don[,5][is.na(don[,5])] <- as.numeric(0)
  sum(is.na(don[,5]))
  str(don)
  
  #don$rollmean <- zoo::rollmean(don[,5],50,fill=NA)
  don$rollmean <- zoo::rollmean(don[,5],100,fill=NA)  # Rollmean of 100 SNP loci.
  head(don)
  
  ## Generate plot using ggplot2.
  library(ggplot2)
  library(ggrepel)
  
  print(paste0('Current comparison: ',comparison))
  y_axis <- paste0('dAF_',comparison)
  
  main_plot <- #don %>%
    ggplot(don, aes(x=BPcum, y=!!ensym(y_axis), col=as.factor(CHR))) + # Never use [ or $ inside aes(). ensym creates a symbol 
    #from the string contained in a variable (so we first have to create those variables at the start of the function), then !! unquotes it, 
    #which means it will work as if you had fed the function raw names.
    #aes(x=row,y=fst,col=chr_ordered)
    
    # Show all points and rollmean line.
    #geom_point(aes(color=as.factor(CHR)), size=0.9, shape=21, stroke=0.2) +
    geom_point(size=0.9, shape=21, stroke=0.2) +
    scale_color_manual(values=rep(c("grey","skyblue"), sum(unique(don$CHR)))) + #("grey30","grey80")
    geom_line(aes(y=rollmean), lwd=0.15, col="black") +
    
    # Customize X axis.
    scale_x_continuous(label=axisdf$CHR, breaks=axisdf$center) +
    #scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
    scale_y_continuous(expand=c(0,0),limits=c(0,1), breaks=seq(0,1,0.1), minor_breaks=NULL) +
    
    # Add title.
    ggtitle(paste0('Manhattan plot of dAF [',nrow(poolData),' SNPs] \n', comparison)) +
    
    # Customize the theme.
    theme_bw() +
    theme(legend.position="none",
          panel.border = element_blank(),
          plot.title = element_text(color="black", size=16, face="bold", hjust = 0.5),
          axis.text.x = element_text(angle=90, color="black"),
          panel.grid = element_blank(),
          panel.grid.major.y=element_line(color="grey60",size=0.2),
          panel.grid.minor.y=element_line(color="grey60",size=0.1),
          axis.title.y = element_text(size=16),
          axis.text = element_text(size=10),
          #axis.ticks.x=element_blank(),
          axis.ticks.x=element_line(size=0.2),
          axis.ticks.y=element_line(size=0.2)) +
    labs(x="Scaffold", y="dAF (100 SNP loci rollmean)")
  
  print(main_plot)
  
  ## --------- Make Manhattan plot with outlier loci highligthed --------- #
  #plot_out <- main_plot +
  #geom_point(data=subset(don, is_highlight0995=="yes"), shape=21, stroke=0.4, size=0.9, col="orange") +
  #  geom_point(data=subset(don, is_highlight0999=="yes"), shape=21, stroke=0.4, size=0.9, col="red") #+
  #geom_label_repel(data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2)  # Add label using ggrepel to avoid overlapping
  
  #print(plot_out)
  
  dev.off()
  
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
#working_dir <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF'
#input_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF/HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.miss20.maf0.05.Neff.AF.500K.txt'
#superpools_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF/List_paired_comparisons_HOM_2020-02-21.csv'

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

## Replace the original scaffold name by a unique number (chromosome for plotting).
## For the HOM samples (when not all the scaffolds are represented and in ascending order).
library(dplyr)
library(tidyr)
num <- summary_dAF_df %>%
  separate(CHR, c("tarseq", "scaffoldNum"), "_")

summary_dAF_df$CHR <- num$scaffoldNum
head(summary_dAF_df)
#summary_dAF_df[summary_dAF_df$CHROM == 'tarseq_85', ]

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
  make.Manhattanplot(summary_dAF_df=summary_dAF_df, 
                     filename=filename, 
                     comparison=comparison)
  
}

print(head(summary_dAF_df))

## Save the dataframe with all dAF values in a TXT file.
library(tools)
superpools_filename <- file_path_sans_ext(basename(superpools_file))
filename
superpools_filename

write.table(summary_dAF_df, file=paste0('Summary_dAF.',superpools_filename,'.',filename,'.txt'), 
            row.names=FALSE, sep='\t', quote=FALSE)
