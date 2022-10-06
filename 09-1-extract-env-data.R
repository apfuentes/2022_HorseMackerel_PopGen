# Extract oceanographic data layers from Bio-Oracle v.2.1
# https://www.bio-oracle.org/code.php
# Modifications by: Angela Fuentes Pardo, e-mail: apfuentesp@gmail.com
# Date: March 04, 2021

# The data available in Bio-ORACLE are documented in two peer reviewed articles that you should cite:
  
# [1] Tyberghein L, Verbruggen H, Pauly K, Troupin C, Mineur F, De Clerck O (2012) Bio-ORACLE: A global environmental dataset 
# for marine species distribution modelling. Global Ecology and Biogeography, 21, 272–281.
# [2] Assis, J., Tyberghein, L., Bosh, S., Verbruggen, H., Serrão, E. A., & De Clerck, O. (2017). Bio-ORACLE v2.0: Extending 
# marine data layers for bioclimatic modelling. Global Ecology and Biogeography.

# Clean environment space
rm(list=ls())

# Set environment variables
working_dir <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/09-GEA'
coord_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/09-GEA/data/HOM_12pools_coordinates.txt'

# Set working directory
setwd(working_dir)

# Load packages 
#library(sdmpredictors) 
#library(leaflet)

# Explore datasets in the package 
#list_datasets() 

# Explore layers in a dataset 
#list_layers() 

# --------------------------------------------------------------------------------
# Extract environmental information for a set of sites
# --------------------------------------------------------------------------------

# Load packages (leaflet allows to load google maps) 
library(sdmpredictors) 
library(leaflet) 

# List layers avaialble in Bio-ORACLE v2 
layers.bio2 <- list_layers(datasets="Bio-ORACLE") 
#write.table(layers.bio2,"layers.bio2.txt", sep="\t", col.names = TRUE, quote=FALSE) 

# Download environmental data layers (Max. Temperature, Min. Salinity and Min. Nitrates at the sea bottom) 
environment.bottom <- load_layers(layercodes = c("BO21_dissoxmax_bdmean","BO21_dissoxmean_bdmean","BO21_dissoxmin_bdmean","BO21_dissoxrange_bdmean",
                                                 "BO21_salinitymax_bdmean","BO21_salinitymean_bdmean","BO21_salinitymin_bdmean","BO21_salinityrange_bdmean",
                                                 "BO21_tempmax_bdmean","BO21_tempmean_bdmean","BO21_tempmin_bdmean","BO21_temprange_bdmean"),
                                  equalarea=FALSE, rasterstack=TRUE) 

# Download bathymetry 
bathymetry <- load_layers("BO_bathymean") 

# Load data frame with the sites of interest 
my.sites <- read.table(coord_file, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
my.sites 

# Visualise sites of interest in google maps 
m <- leaflet() 
m <- addTiles(m) 
m <- addMarkers(m, lng=my.sites$Lon, lat=my.sites$Lat, popup=my.sites$Name, icon=TRUE) 
m 

# Extract environmental values from layers 
env <- data.frame(Sample=my.sites$Sample,
                  #Code=my.sites$Code,
                  depth=extract(bathymetry, my.sites[ ,c("Lon","Lat")]), 
                  extract(environment.bottom, my.sites[ ,c("Lon","Lat")])
                  ) 
env

# Save env data in a file
write.table(env, "Environmental_data_HOM_12pools.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# --------------------------------------------------------------------------------
# Filter environmental datasets
# --------------------------------------------------------------------------------
library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RDA

# Examine whether there is high correlation among predictors
head(env)
str(env)
pdf('Correlation_env_HOM_12pools.pdf')
pairs.panels(env[,3:ncol(env)], scale=TRUE)
dev.off()

# RDA is a regression-based method, and so can be subject to problems when using highly correlated predictors (Dormann et al., 2013). 
# Generally, the |r| > 0.7 “rule of thumb” is a good guideline for removing correlated predictors. We will also check for 
# multicollinearity using Variance Inflation Factors (VIF), below.

# Given that high correlation was observed between max, min, and mean values, I decided to retain only mean and range values, as they 
# would be better descriptors of environmental conditions at a given site
to_exclude <- c("BO21_dissoxmax_bdmean","BO21_dissoxmin_bdmean",#"BO21_dissoxmean_bdmean","BO21_dissoxrange_bdmean",
                "BO21_salinitymin_bdmean","BO21_salinitymax_bdmean",#"BO21_salinitymean_bdmean","BO21_salinityrange_bdmean",
                "BO21_tempmax_bdmean","BO21_tempmin_bdmean"#,"BO21_tempmean_bdmean","BO21_temprange_bdmean"
                )
env_set <- env[ ,!names(env) %in% to_exclude]

str(env_set)
pdf('Correlation_env_HOM_12pools_filtered.pdf')
pairs.panels(env_set[,3:ncol(env_set)], scale=TRUE)
dev.off()

# Save filtered env data in a file
write.table(env_set, "Environmental_data_HOM_12pools_filtered.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# --------------------------------------------------------------------------------
# Explore the selected environmental data layers
# --------------------------------------------------------------------------------
# Load package 
library(sdmpredictors) 

env_list <- colnames(env_set[,3:ncol(env_set)])
#env_list <- c("BO21_salinitymin_bdmean","BO21_salinitymean_bdmean","BO21_salinitymax_bdmean","BO21_salinityrange_bdmean")
#env_list <- c("BO21_tempmin_bdmean","BO21_tempmean_bdmean","BO21_tempmax_bdmean","BO21_temprange_bdmean")
#env_list <- c("BO21_dissoxmin_bdmean","BO21_dissoxmean_bdmean","BO21_dissoxmax_bdmean","BO21_dissoxrange_bdmean")
env_list

# Make the plots for all the variables of interest
pdf("Maps_environmental_data_HOM_12pools_filtered_sst_pal.pdf",7,10)
#pdf("Maps_environmental_data_HOM_12pools_only_salinity.pdf")
#pdf("Maps_environmental_data_HOM_12pools_only_temperature.pdf")
#pdf("Maps_environmental_data_HOM_12pools_only_O2.pdf")
par(mfrow=c(3,3))
#par(mfrow=c(2,2))

for(i in env_list){
  #i="BO21_salinityrange_bdmean"
  #i="BO21_tempmean_bdmean"
  
  # Easy download of raster file 
  env_layer <- load_layers(i) #BO21_tempmin_bdmean
  
  # Crop raster to fit the area of interest
  atlantic.ext <- extent(-23,20,18,69)  # min Log, max Log, min Lat, max Lat
  env_layer.crop <- crop(env_layer, atlantic.ext) 
  
  # Choose a color palette 
  library(palr)
  #my.colors <- colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))
  sst_pal <- sst_pal(palette=TRUE)
  my.colors <- colorRampPalette(sst_pal$cols) 
  
  # Plot the map 
  plot(env_layer.crop, col=my.colors(1000), axes=FALSE, box=FALSE, main=i) 
  #title(cex.sub = 1.25, sub = i) 
}
dev.off()

