# Obtain density plots for the annotations extracted from the VCF file. To run in R

library(dplyr) # include line to install the package if not available
library(ggplot2) # include line to install the package if not available

# Set working directory.
setwd('~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/05-snp-filtering/')

# Load data table with GATK annotations.
data <- read.csv('HOM_12_pools.SNPs.raw.ANN.table',sep='\t') # this might take a while (1.5G)

# Explore the column names of this table.
names(data)
#[1] "ID"             "CHROM"          "POS"            "FILTER"         "QD"             "FS"             "SOR"           
#[8] "MQ"             "MQRankSum"      "ReadPosRankSum"

# ------------ plots for normal X axis, one for each annotation of interest ---------------------------

# Save all the pdf files.
pdf('summary_GATK_QC_ANN_HOM_12pools_raw_SNPs.pdf')
#pdf('summary_GATK_QC_ANN_HOM_12pools_raw_SNPs.pdf', width = 7, height = 3)

ggplot(data, aes(x=QD)) + geom_density(fill=rgb(0,0,.5,.5))+
  geom_vline(xintercept = c(1.0)) + xlim(0, 40)
#scale_x_continuous(limits = c(-4,4))
#ggsave(paste('summary_QD_biallelic_SNPs.pdf',sep=''), width = 7,height = 3)

ggplot(data, aes(x=FS)) + geom_density(fill=rgb(0,0,.5,.5))+
  geom_vline(xintercept = c(60.0)) + xlim(0, 200) + scale_x_log10()
#scale_x_continuous(limits = c(-4,4))
#ggsave(paste('summary_FS_biallelic_SNPs_xLog10.pdf',sep=''), width = 7,height = 3)

ggplot(data, aes(x=SOR)) + geom_density(fill=rgb(0,0,.5,.5))+
  geom_vline(xintercept = c(3.0)) + xlim(0, 10)
#scale_x_continuous(limits = c(-4,4))
#ggsave(paste('summary_SOR_biallelic_SNPs.pdf',sep=''), width = 7,height = 3)

ggplot(data, aes(x=MQ)) + geom_density(fill=rgb(0,0,.5,.5))+
  geom_vline(xintercept = c(40.0)) + xlim(20, 70)
#scale_x_continuous(limits = c(-4,4))
#ggsave(paste('summary_MQ_biallelic_SNPs.pdf',sep=''), width = 7,height = 3)

ggplot(data, aes(x=MQRankSum)) + geom_density(fill=rgb(0,0,.5,.5))+
  geom_vline(xintercept = c(-12.5)) + xlim(-15, 15)
#scale_x_continuous(limits = c(-4,4))
#ggsave(paste('summary_MQRankSum_biallelic_SNPs.pdf',sep=''), width = 7,height = 3)

ggplot(data, aes(x=ReadPosRankSum)) + geom_density(fill=rgb(0,0,.5,.5))+
  geom_vline(xintercept = c(-8.0)) + xlim(-10, 10)
#scale_x_continuous(limits = c(-4,4))
#ggsave(paste('summary_ReadPosRankSum_biallelic_SNPs.pdf',sep=''), width = 7,height = 3)

dev.off()
