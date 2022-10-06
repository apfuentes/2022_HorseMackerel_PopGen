###########################################################################
# Test for Isolation-by-Distance (IBD) and Isolation-by-Environment (IBE) #
###########################################################################

## Clean environment space.
rm(list=ls())

working_dir <- "/Users/angfu103/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/09-GEA/IBD"
coord_txt_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/09-GEA/data/HOM_12pools_coordinates.txt'
coord_csv_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/09-GEA/coordinates_HOM_9pools.csv'
FST_file <- "/Users/angfu103/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/07-1-popgen-stats/pairwiseFST/PairwiseFST_HOM_12_pools.UG.SNPs.hf.DP20-300.GQ20.mono.miss20.mac3.AD.Neff.AF_Feb12011_renamed.txt"

setwd(working_dir)  # Set working directory.

####################################################
# Obtain pairwise least-coast distance (Km) matrix #
####################################################

## Taken from CartDist, function CartDistFunction.R, lines 1 to 138.

#coordinates - path to a csv file with the coordinates and site names. This must have three columns the first being the site name, second the longitude and third the latitude.
#min.depth   - minimum depth that can be considered passable in the least-cost analysis. Note that depths are negative. Default is NA or 0.
#max.depth   - maximum depth that can be considred passable in the least-cost analysis. Note that depths are negative. Default is NA or NULL meaning the analysis depths will permit movement between the min.depth and the maximum depth of the derived bathymetric layer.
#trans       - is a transition object (marmap) that might have been calcualted in a previous analysis. Note that each run will return the transition object. This object can called from the workspace into the fuction. 
#gridres     - resolution (mins) of the bathymetric grid used for least coast path (default = 2)
#directory   - if specified this is where marmap will save and or look for the output for the bathy object created by marmap. 
#outpath     - this is the filepath for the output from the function a .Rdata file. Note this must be a full file path ending in .RData. If no path is provided a 'Output' .RData file will be created in the current directory.

coordinates=coord_csv_file
min.depth=0
max.depth=NA
trans=NA
gridres=2
#gridres=4
directory=working_dir
outpath=paste0(working_dir,"/LeastCoastDistanceMatrix_HOM_9pools.RData")

#Libraries ----------
#Check to make sure the packages required are there and install missing

writeLines("\nChecking on package availability.\n")
packages <- c("gdistance", "ggplot2", "marmap","vegan")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) { 
  writeLines(paste("Installing missing packages: ",paste(setdiff(packages, rownames(installed.packages())),collapse=", "),sep=" "))
  install.packages(setdiff(packages, rownames(installed.packages())))  
} 

#load each library
require(gdistance)
require(marmap)
require(vegan)
require(ggplot2)

## check the outpath
if(!is.na(outpath) & substring(outpath,nchar(outpath)-5,nchar(outpath))!= ".RData"){
  stop("\nParamter outpath must be a full path ending in .RData. Please fix and try again.\n")
}

# Plotting settings for water bathymetry fill
blues <- c("lightsteelblue4", "lightsteelblue3",
           "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

## read in coordinates and check format
coords <- read.csv(coordinates, header = TRUE, stringsAsFactors = FALSE)

writeLines("\nEnsure that your data is set up as three columns with the first column having the sample location name, the second with the longitude and the third with the latitude. \n")

if (length(colnames(coords))<3 | !is.numeric(coords[,2]) | !is.numeric(coords[,3])){
  stop("\nCheck coordinate input, there is a problem here.\n")
}

# Clean up column names
colnames(coords)[which(sapply(coords,is.character))[1]] <- "Code"
colnames(coords)[which(sapply(coords,is.numeric))[1]] <- "Long"
colnames(coords)[which(sapply(coords,is.numeric))[2]] <- "Lat"

## Set map limits------adding and subtracting 2 degrees to make the lc.dists function work better
Lat.lim = c(min(coords$Lat)-2, max(coords$Lat)+2)
Long.lim = c(min(coords$Long)-2, max(coords$Long)+2)

#Get the bathydata and keep it
holddir <- getwd() # grab working directory

# directory is specified.
if(!is.na(directory)){
  setwd(directory) # switch to the output directory
  writeLines("\nGetting bathymetry data from NOAA database\n")
  bathydata <- marmap::getNOAA.bathy(lon1 = Long.lim[1], lon2 = Long.lim[2], lat1 = Lat.lim[1], lat2 = Lat.lim[2],
                                   resolution = gridres,keep=TRUE)
  setwd(holddir)
} # set back to the working directory

# directory is not specified.
if(is.na(directory)){
  writeLines("\nGetting bathymetry data from NOAA database\n")
  bathydata <- marmap::getNOAA.bathy(lon1 = Long.lim[1], lon2 = Long.lim[2], lat1 = Lat.lim[1], lat2 = Lat.lim[2],
                                   resolution = gridres,keep=FALSE)
}


### check depths to ensure they are all in water (e.g., no positive 'land' depths)

# Get depths and plot. If any depths > 0 we will not proceed
depths <- marmap::get.depth(bathydata, x=coords$Long, y=coords$Lat, locator=F)
colnames(depths) <- c("Long","Lat","depth")
coords <- merge(coords, depths, by=c("Long","Lat"))

#colours to assign those locations which are in water "green" and on land "red". 
coords$col <- "green"
coords[coords$depth>=0,"col"] <- "red"
#coords[is.na(coords$depth),"col"] <- "red"

if(length(which(depths$depth>0))!=0){
  
  writeLines("\nSome of your coordinates appear to have a positive depth. Refer to map and bump coordinates for those points marked as 'red' off of land. Suggest moving points farther off land for this analysis.\n")
  
  marmap::plot.bathy(bathydata,image = TRUE, land = T, lwd = 0.03,
                     bpal = list(c(0, max(bathydata), greys),
                                 c(min(bathydata), 0, blues)),deep=0,shallow=0)
  
  marmap::plot.bathy(bathydata, lwd = 1, deep = 0, shallow = 0, step = 0, add = TRUE)
  
  legend("bottomright",
         legend = c("Water","Land"), 
         col=c("green","red"),
         pch=19,
         pt.cex=1.5,
         bg="white")
  
  points(coords$Long, coords$Lat, pch=19, cex=2, col=coords$col)
  
  print(coords[coords$depth>0,c("Code","Long","Lat","depth")])
  
  stop("\nFix and re-run funciton.\n")
  
}

# min and maximum depth for transition object

if(is.na(min.depth)){ min.depth <- 0
writeLines("\n No min.depth specified, defaulting to 0 m depth.\n")
}

if(is.na(max.depth)){ max.depth <- NULL
writeLines("\n No max.depth specified, defaulting to maximum depth of bathymetry object from marmap.\n")
}

#Make the trans mat object then do the lc dist calculation

if(is.na(trans)){
  writeLines("\nCalculating transition object for least-cost analysis.\n")
  trans <- marmap::trans.mat(bathydata, min.depth = min.depth, max.depth = max.depth) 
}

sites <- coords[, c("Long","Lat")]
rownames(sites) <- coords[, "Code"]
sites <- sites[order(sites$Lat, decreasing = TRUE), ]

writeLines("\nCalculating least cost distances. This will probably can take a few minutes depending on resolution. If insufficient memory error is returned try adjusting the gridres argument (default = 2) to a higher number. gridres refers to the resolution of the bathymetric grid in minutes.\n")
lc.dists <- marmap::lc.dist(trans, sites, res="dist")
lc.dists

## Transform the dist class object into a matrix for subsequent analysis.
lcDist <- as.matrix(lc.dists)
lcDist

# Remove Mediterranean sample
tmp <- as.data.frame(lcDist)
tmp$MED <- NULL
tmp <- tmp[!rownames(tmp) == "MED", ]
tmp <- as.matrix(tmp)
lcDist <- tmp

# Rename rows and columns (I got the pop order from excel using the original coordinate locations).
#rownames(lcDist) <- c("Tri", "NDS", "NDF", "Lab", "Bla", "Ste", "BDO", "NsF", "Mus", "Mir", "NsS", "ScB", "InB", "7iA", "GeB", "M514")
#colnames(lcDist) <- c("Tri", "NDS", "NDF", "Lab", "Bla", "Ste", "BDO", "NsF", "Mus", "Mir", "NsS", "ScB", "InB", "7iA", "GeB", "M514")
lcDist

## Remove teo pops not used in the final analysis.
#Toremove <- c("InB", "NsF")
#lcDist <- lcDist[!(rownames(lcDist) %in% Toremove), !(colnames(lcDist) %in% Toremove)]
#lcDist

## Save the matrix with least coast distances (Km) as a CSV file.
write.table(lcDist, file = "LeastCoastDistanceMatrix_HOM_9pops.csv", row.names = TRUE, col.names = NA, sep = ',', quote = FALSE)


# Load data ---------------------------------------------------------------

# Load libraries
library(data.table)
library(tidyverse)
library(tools)

# Load the previously created lcDist file
lcDist <- read.table("LeastCoastDistanceMatrix_HOM_9pops.csv", header = TRUE, row.names = 1, sep = ',')

# Load pairwise FST matrix
pFst <- read.table(FST_file, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)  # Modify the separator for the one used in your file.

# Exclude some samples to only retain the northern samples
to_include <- c("NOS1", "NOS2", "WIE1", "NES", "NPT1")
#to_include <- c("NOS1", "NOS2", "WIE1", "WIE3", "NES", "NPT1")
#to_include <- c("WIE1", "WIE3", "NES", "NPT1")
pFst <- pFst[rownames(pFst) %in% to_include, ]
pFst <- pFst[, colnames(pFst) %in% to_include]
lcDist <- lcDist[rownames(lcDist) %in% to_include, ]
lcDist <- lcDist[, colnames(lcDist) %in% to_include]

# Reorder matrix to match the rownames and colnames of the distance matrix
pFst <- as.matrix(pFst)
pFst <- pFst[rownames(lcDist), colnames(lcDist)]
pFst

# Verify rownames are euqal between GenDit and GeoDist matrixes
all.equal(rownames(pFst), rownames(lcDist))

# Linearize FST values applying the formula FST/(1-FST)
ln_pFst <- apply(pFst, 2, function(x) x/(1-x))
ln_pFst

## Save the matrix with linerialized FSTs as a CSV file
#write.table(pFst, file="Linearized_pairwiseFST_Herring_14pops_reordered.csv", row.names=TRUE, sep=',', quote=FALSE)


# Run the isolaion-by-distance test ---------------------------------------

# Readings: http://dr-k-lo.blogspot.com/2012/04/how-to-do-mantel-test-in-r.html, 
# https://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/
# The test consists of calculating the correlation of the entries in the matrices, 
# then permuting the matrices and calculating the same test statistic under each permutation and 
# comparing the original test statistic to the distribution of test statistics from the permutations 
# to generate a p-value. The number of permutations defines the precision with which the p-value can be calculated. 
# The function to perform the Mantel test is mantel.rtest and the required arguments are the two distance matrices.  
# The number of permutations can also be specified by the user, but is by default 99.

# Load library
library(ade4)

# Convert the distance matrices to dist class objects
GEOdist <- as.dist(lcDist)
GENdist <- as.dist(ln_pFst)  # linearized FST

# Run the Mantel test
set.seed(777)
mantel.rtest(GEOdist, GENdist, nrepet = 9999)

#Monte-Carlo test
#Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
#Observation: 0.8917953 
#Based on 9999 replicates
#Simulated p-value: 0.0153 
#Alternative hypothesis: greater 
#Std.Obs  Expectation     Variance 
#2.678299239 -0.002933674  0.111600163 
# Accept Ho: there is correlation between genetic distance and geographic distance

# Plot the result
set.seed(777)
ibd <- mantel.rtest(GEOdist, GENdist, nrepet = 9999)
plot(ibd)

# To get the slope you can run a linear model:
summary(lm(GENdist~GEOdist))
lm_2M <- lm(GENdist~GEOdist)   # Save lm for plotting.

# Convert dist objects into matrixes for plotting.
geN <- as.matrix(GENdist)
geO <- as.matrix(GEOdist)

# "Melt" the matrices to a three column format and save them in a dataframe
library(reshape2)
split_geN <- setNames(reshape2::melt(geN), c('rows', 'vars', 'FstDist'))
split_geO <- setNames(reshape2::melt(geO), c('rows', 'vars', 'GeoDist'))
df_geOgeN <- cbind.data.frame(split_geN, split_geO, stringsAsFactors=FALSE)

# Get the equation, R-squared adjusted, and significance values as characters
# Source: https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
lm_eqn <- function(df, y, x){
  formula = as.formula(sprintf('%s ~ %s', y, x))
  m <- lm(formula, data = df);
  # formating the values into a summary string to print out
  # ~ give some space, but equal size and comma need to be quoted
  eq <- substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~italic(pvalue), 
                   #italic(target) == a + b %.% italic(input)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~italic(pvalue),
                   list(#target = y,
                        input = x,
                        a = format(as.vector(coef(m)[1]), digits = 2), 
                        b = format(as.vector(coef(m)[2]), digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3),
                        # getting the pvalue is painful
                        pvalue = format(summary(m)$coefficients[2,'Pr(>|t|)'], digits = 2)
                   )
  )
  as.character(as.expression(eq));                 
}

# Extract and remane columns
df <- df_geOgeN[, c("GeoDist", "FstDist")]
colnames(df) <- c("x", "y")

# Make the plot
p <- ggplot(df, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y~x, linetype = "dashed", size = 0.5, color = "black") +
  geom_text(x = 1450, y = 0, label = lm_eqn(df, 'y','x'), parse = TRUE, size = 6) +
  xlab("Geographic distance (km)") + 
  #ylab("Genetic distance (FST/1-FST)") +
  ylab(expression(paste("Genetic distance (", italic("F")[ST]/1-italic("F")[ST],")"))) +
  theme_bw() +
  theme(legend.position = "right", 
        panel.grid = element_blank(), 
        text = element_text(size = 20, colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
#p

# Make a scatter plot to show the relationship between genetic and geographic distance 
pdf("Isolation_by_distance_test_HOM_Northern_pops_NOS1-NOS2-WIE1-NES-NPT1.pdf", width = 5.95, height = 4.95)
p
dev.off()


# Reading on IBD: https://www.molecularecologist.com/2012/09/04/isolating-isolation-by-distance/
# test for IBD within population clusters identified by AMOVA or a clustering algorithm
# partial Mantel test can be used to test whether geography (and thereby IBD) contributes to apparent clustering—by 
# testing to see whether an association between a matrix of cluster assignments and genetic distances disappears 
# when controlled for geographic distance.

GENdist
geN
GEOdist
geO
euclidean_distance <- geO


library(vegan)
eu_dist <- vegdist(geO, method="euclidean")

mantel.partial(geN, geO, geO, method="pearson", permutations=999)
