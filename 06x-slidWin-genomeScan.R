################################
# Sliding windowed-genome scan #
################################

# Clean environment space.
rm(list=ls())

# Set working directory.
#dir.create('/Users/angfu103/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-outlier-loci-detection/HOM-12pools-Neff-miss20-maf001')
setwd('/Users/angfu103/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/06x-slidWin-genomeScan')

######################
# Data preprocessing #
######################

# Load input data.
library(data.table) 
path_to_file <- '/Users/angfu103/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/06x-slidWin-genomeScan/test_data.txt'

SnpFST.pairwise.round <- fread(path_to_file, data.table=FALSE, header=TRUE, sep='\t', stringsAsFactors=FALSE)
head(SnpFST.pairwise.round)  # Explore first lines of the dataframe.

# Create a dataframe with CHROM and POS info plus the FST by SNP of all pairwise comparisons. 
SnpFST_pairwise_df <- cbind.data.frame(CHR = SnpFST.pairwise.round$CHROM, 
                                       CHROM = SnpFST.pairwise.round$CHROM, 
                                       POS = SnpFST.pairwise.round$POS, 
                                       CHROM_POS = SnpFST.pairwise.round$CHROM_POS, 
                                       SnpFST.pairwise.round[,4:ncol(SnpFST.pairwise.round)], 
                                       stringsAsFactors = FALSE)

SnpFST_pairwise_df[1:10, 1:7]
dim(SnpFST_pairwise_df)
str(SnpFST_pairwise_df)

# Replace scaffold names (characters) by numbers (for plotting).
# Obtain a list of unique scaffold names.
uq <- unique(SnpFST_pairwise_df$CHR)

# Create a vector with a sequence of numbers in ascending order from 1 to the number of unique scaffolds.
num <- seq(1, length(uq))

# Safe in a TXT file the scaffold = chromosome number equivalence.
library(tools)
filename <- file_path_sans_ext(basename(path_to_file))
write.table(cbind.data.frame(uq, num), file = paste0('Scaffold-number_equivalence_Manhattan_plot.',filename, '.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

# Replace the original scaffold/contig name for a number (chromosome for plotting).
for (i in seq_along(num)) { SnpFST_pairwise_df$CHR[SnpFST_pairwise_df$CHR == uq[i]] <- num[i] }

SnpFST_pairwise_df[1:10, 1:7]

# Convert columns CHR and BP from character to numeric values.
SnpFST_pairwise_df$CHR <- as.numeric(SnpFST_pairwise_df$CHR)
SnpFST_pairwise_df$POS <- as.numeric(SnpFST_pairwise_df$POS)

str(SnpFST_pairwise_df)

# As Variables names can't start with a number, use make.names() to pre-pend an "X" before colnames beginning with a digit.
colnames(SnpFST_pairwise_df) <- make.names(colnames(SnpFST_pairwise_df))

head(SnpFST_pairwise_df)

################################
# Sliding windowed-genome scan #
################################

# Set window size.
window_size <- 5

# Save in a vector the list of comparison pairs.
pairs <- colnames(SnpFST_pairwise_df[5:ncol(SnpFST_pairwise_df)])
#pairs

# Open a PDF device to store the plot.
pdf(file=paste0('Manhattan_plot_Pairwise_FST_by_SNP.',filename,'.pdf'), width = 11, height = 8.5, paper = 'a4r')

maxrows <- length(pairs)
if (maxrows < 4) {
  par(mfrow = c(maxrows, 1))  ## Set the layout to be X rows by 1 column.
} else {
  par(mfrow = c(4, 1))  ## Set the layout to be max. 6 rows by 1 column.
}

for (comparison in pairs) {
  #comparison <- 'X1a_vs_2a'
  #comparison <- 'X1a_West_Ireland_2016_vs_1b_Southwest_Ireland_2016'
  
  new_col <- c()
  
  for (scaff in unique(SnpFST_pairwise_df$CHROM)) {
    #scaff <- 'tarseq_0'
    #df_scaff <- SnpFST_pairwise_df[SnpFST_pairwise_df$CHROM == scaff, ]

    df_scaff <- SnpFST_pairwise_df[SnpFST_pairwise_df$CHROM == scaff, c(comparison), drop = FALSE]
    
    smoothed <- filter(df_scaff, rep(1/window_size, window_size))
    #smoothed <- filter(df_scaff[, c(comparison)], rep(1/window_size, window_size))
    new_col <- append(new_col, smoothed)
    
  }
  
  # Create a unified dataframe, following the specifications of the qqman R package (http://www.gettinggeneticsdone.com/2014/05/qqman-r-package-for-qq-and-manhattan-plots-for-gwas-results.html).
  #SNP CHR BP         P    
  #1 rs1   1  1 0.9148060 
  #2 rs2   1  2 0.9370754 
  #3 rs3   1  3 0.2861395 
  #4 rs4   1  4 0.8304476 
  #5 rs5   1  5 0.6417455 
  #6 rs6   1  6 0.5190959 
  plot_df <- cbind.data.frame(SNP=SnpFST_pairwise_df$CHROM_POS, CHR=SnpFST_pairwise_df$CHR, BP=SnpFST_pairwise_df$POS, VALUE=new_col)
  
  # Assign NAs added by filter() at the edges of the scaffolds to zeroes.
  plot_df[, c('VALUE')][is.na(plot_df[, c('VALUE')])] <- 0
  head(plot_df)
  
  # Load modified R function of the package qqman https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html
  source("~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/code/utility-code/makeManhattanPlot.R")

  # Obtain a Manhattan plot for the values of interest. 
  makeManhattanPlot(plot_df, chr="CHR", bp="BP", val="VALUE", mode='smfst',
            #cex = 0.7, cex.axis = 1.1, 
            col = c("#bdbdbd", "#636363"), #"#636363"
            logp = FALSE, suggestiveline = FALSE, genomewideline = FALSE)
  title(main = paste0('Pairwise-FST by SNP, comparison: ', comparison), adj = 1)

  
}  

dev.off() 
par(mfrow=c(1,1)) # Set layout back to default.  
  
  
  
  

########################
# OLD VERSION

#######################################################################################
# Genome-wide examination of pairwise comparisons following a sliding window approach #
#######################################################################################

# Clean environment space.
rm(list=ls())

# Set working directory.
#dir.create('/Users/angfu103/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-outlier-loci-detection/HOM-12pools-Neff-miss20-maf001')
setwd('/Users/angfu103/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/06x-slidWin-genomeScan')

######################
# Data preprocessing #
######################

# Load input data.
library(data.table) 
path_to_file <- '/Users/angfu103/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/06x-slidWin-genomeScan/test_data.txt'

df <- fread(path_to_file, data.table=FALSE, header=TRUE, sep='\t', stringsAsFactors=FALSE)
head(df)  # Explore first lines of the dataframe.
dim(df)
str(df)

# Split locus name by scaffold and position.
library(tidyr)
df_split <- separate(data = df, col = "CHROM-POS", into = c("CHROM", "POS"), sep = "-")
head(df_split)

# Transpose dataframe, t() function will generate a matrix.
#tmp <- t(df)
#class(tmp)
#tmp[1:5,]

# Subset dataset, convert back to dataframe.
#data_sub <- as.data.frame(tmp[1:100,2])
#dim(df_sub)
#data_sub[1:5,]
#class(df_sub)

# Convert missing data encoded as -99 for NA. >>>>>>>>>>>>> HOW filter() HANDLES NAs? NEED TO REMOVE THEM??
#sum(df_sub=='-99')
#df_sub[df_sub=='-99'] <- NA
#sum(is.na(df_sub))
#str(df_sub)
#str(df_sub)

# Remove locus (columns) with missing values.
#df_clean <- df_sub[, !is.na(colSums(df_sub != 0)) & colSums(df_sub != 0) > 0]
#dim(df_clean)


################################
# Sliding windowed-genome scan #
################################

# Set window size.
window_size <- 5
new_col <- c()
#step <- 10

for (scaff in unique(df_split$CHROM)) {
  #scaff <- 'tarseq_0'
  df_scaff <- df_split[df_split$CHROM == scaff, ]
  
  smoothed <- filter(df_scaff$fst_AvsB, rep(1/window_size, window_size))
  
  new_col <- append(new_col, smoothed)
  
  #filter(df_split[df_split$CHROM == scaff, ], rep(1/wd_size, wd_size), sides = 2)
  # sides = 2 is equivalent to align="center" for the zoo::rollmean or RcppRoll::roll_mean. sides = 1 is equivalent to "right" alignment. 
  #ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}   # 
  
  #win.mean <- sum(x[i:(i+2)], na.rm = TRUE)/win.size
}

df_split$Fst_smoothed <- new_col

df_split$CHROM_POS <- paste(df_split$CHROM, df_split$POS, sep="-")

##########################
# Obtain Manhattan plots #
##########################
library(ggplot2)

ggplot(df_split, aes_string('CHROM_POS', 'fst_AvsB')) +
  geom_point(size = 2.5, shape = 16, colour = 'black') +
  theme(panel.border = element_rect(colour = 'black', fill = NA, size = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text.x = element_text(colour='white'))

ggplot(df_split, aes_string('CHROM_POS', 'Fst_smoothed')) +
  geom_point(size = 2.5, shape = 16, colour = 'red') +
  theme(panel.border = element_rect(colour = 'black', fill = NA, size = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text.x = element_text(colour='white'))


# OR, Obtain Manhattan plots wihtin a loop.
library(ggplot2)

p <- 
  ggplot(SnpFST_pairwise_df, aes_string('CHROM_POS', comparison)) +
  geom_point(size = 2.5, shape = 16, colour = 'black') +
  theme(panel.border = element_rect(colour = 'black', fill = NA, size = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle(paste0('FST by SNP pairwise comparison: ', comparison)) +
  labs(x = 'Scaffold', y = expression(paste('Smoothed-',italic(F)[ST])))

print(p)

q <- 
  ggplot(SnpFST_pairwise_df, aes_string('CHROM_POS', 'Fst_smoothed')) +
  geom_point(size = 2.5, shape = 16, colour = 'red') +
  theme(panel.border = element_rect(colour = 'black', fill = NA, size = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle(paste0('FST by SNP pairwise comparison: ', comparison)) +
  labs(x = 'Scaffold', y = expression(paste('Smoothed-',italic(F)[ST])))

print(q)

