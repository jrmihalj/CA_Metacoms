# CA parasite metacommunities
# Lat_Long metacommunities divided by year
# August 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Subset into latitudinal metacommunities:

latNorth <- subset(allsites_latlong, Long > -121.63)
nrow(latNorth) # 17 sites

latMidNorth <- subset(allsites_latlong, Long < -121.63 & Long > -121.86)
nrow(latMidNorth) # 99 sites

latMidSouth <- subset(allsites_latlong, Long < -121.86 & Long > -122.04)
nrow(latMidSouth) # 81 sites

latSouth <- subset(allsites_latlong, Long < -122.04)
nrow(latSouth) # 41 sites

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# The latNorth doesn't have very many ponds, so I'll just do the two mid-latitudes

#############################################################################################
#############################################################################################

# latMidNorth by year: 2009

latMidNorth$site.yr <- as.character(latMidNorth$site.yr)
latMidNorth_2009 <- subset(latMidNorth, 
                          substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2009")

# Need to remove some columns
latMidNorth_2009 <- latMidNorth_2009[, -c(1, 18:19)]
head(latMidNorth_2009)

# Turn into a matrix
latMidNorth_2009_mat <- as.matrix(latMidNorth_2009)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(latMidNorth_2009_mat)){
  if(sum(latMidNorth_2009_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
latMidNorth_2009_mat <- latMidNorth_2009_mat[, -zeros]

nrow(latMidNorth_2009_mat) # 31 sites
ncol(latMidNorth_2009_mat) # 10 parasites 

#### METACOMMUNITY ANALYSES #####
library(metacom)

# Run the gamet of metacommunity analyses to determine structure
latMidNorth_2009_metacom <- metacommunity(latMidNorth_2009_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
latMidNorth_2009_Comm <- latMidNorth_2009_metacom$Comm
latMidNorth_2009_plot <- Matrix_Plot(latMidNorth_2009_Comm)
quartz(height=6, width=3)
print(latMidNorth_2009_plot) 
# looks somewhat nested

# View Coherence
latMidNorth_2009_metacom$Coherence
# p-val = 0.036, positive coherence

latMidNorth_2009_metacom$Turnover
# p-val = 0.18, no sig. turnover, but less replacements than expected (quasi-nested)

latMidNorth_2009_metacom$Boundary
# 1.65, p<0.0001, clumped losses

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# latMidNorth by year: 2010

latMidNorth$site.yr <- as.character(latMidNorth$site.yr)
latMidNorth_2010 <- subset(latMidNorth, 
                           substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2010")

# Need to remove some columns
latMidNorth_2010 <- latMidNorth_2010[, -c(1, 18:19)]
head(latMidNorth_2010)

# Turn into a matrix
latMidNorth_2010_mat <- as.matrix(latMidNorth_2010)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(latMidNorth_2010_mat)){
  if(sum(latMidNorth_2010_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
latMidNorth_2010_mat <- latMidNorth_2010_mat[, -zeros]

nrow(latMidNorth_2010_mat) # 38 sites
ncol(latMidNorth_2010_mat) # 10 parasites 

#### METACOMMUNITY ANALYSES #####
library(metacom)

# Run the gamet of metacommunity analyses to determine structure
latMidNorth_2010_metacom <- metacommunity(latMidNorth_2010_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
latMidNorth_2010_Comm <- latMidNorth_2010_metacom$Comm
latMidNorth_2010_plot <- Matrix_Plot(latMidNorth_2010_Comm)
quartz(height=6, width=3)
print(latMidNorth_2010_plot) 
# doesn't look like anything

# View Coherence
latMidNorth_2010_metacom$Coherence
# p-val = 0.29, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# latMidNorth by year: 2011

latMidNorth$site.yr <- as.character(latMidNorth$site.yr)
latMidNorth_2011 <- subset(latMidNorth, 
                           substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")

# Need to remove some columns
latMidNorth_2011 <- latMidNorth_2011[, -c(1, 18:19)]
head(latMidNorth_2011)

# Turn into a matrix
latMidNorth_2011_mat <- as.matrix(latMidNorth_2011)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(latMidNorth_2011_mat)){
  if(sum(latMidNorth_2011_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
latMidNorth_2011_mat <- latMidNorth_2011_mat[, -zeros]

nrow(latMidNorth_2011_mat) # 22 sites
ncol(latMidNorth_2011_mat) # 10 parasites 

#### METACOMMUNITY ANALYSES #####
library(metacom)

# Run the gamet of metacommunity analyses to determine structure
latMidNorth_2011_metacom <- metacommunity(latMidNorth_2011_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
latMidNorth_2011_Comm <- latMidNorth_2011_metacom$Comm
latMidNorth_2011_plot <- Matrix_Plot(latMidNorth_2011_Comm)
quartz(height=6, width=3)
print(latMidNorth_2011_plot) 
# doesn't look like anything

# View Coherence
latMidNorth_2011_metacom$Coherence
# p-val = 0.21, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

#### NOT ENOUGH SITES SAMPLED IN 2012 FOR THIS AREA #####


#############################################################################################
#############################################################################################

# latMidSouth by year: 2009

latMidSouth$site.yr <- as.character(latMidSouth$site.yr)
latMidSouth_2009 <- subset(latMidSouth, 
                           substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2009")

# Need to remove some columns
latMidSouth_2009 <- latMidSouth_2009[, -c(1, 18:19)]
head(latMidSouth_2009)

# Turn into a matrix
latMidSouth_2009_mat <- as.matrix(latMidSouth_2009)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(latMidSouth_2009_mat)){
  if(sum(latMidSouth_2009_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
latMidSouth_2009_mat <- latMidSouth_2009_mat[, -zeros]

nrow(latMidSouth_2009_mat) # 11 sites
ncol(latMidSouth_2009_mat) # 6 parasites 

#### METACOMMUNITY ANALYSES #####
library(metacom)

# Run the gamet of metacommunity analyses to determine structure
latMidSouth_2009_metacom <- metacommunity(latMidSouth_2009_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
latMidSouth_2009_Comm <- latMidSouth_2009_metacom$Comm
latMidSouth_2009_plot <- Matrix_Plot(latMidSouth_2009_Comm)
quartz(height=6, width=3)
print(latMidSouth_2009_plot) 
# looks somewhat nested

# View Coherence
latMidSouth_2009_metacom$Coherence
# p-val = 0.91, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# latMidSouth by year: 2010

latMidSouth$site.yr <- as.character(latMidSouth$site.yr)
latMidSouth_2010 <- subset(latMidSouth, 
                           substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2010")

# Need to remove some columns
latMidSouth_2010 <- latMidSouth_2010[, -c(1, 18:19)]
head(latMidSouth_2010)

# Turn into a matrix
latMidSouth_2010_mat <- as.matrix(latMidSouth_2010)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(latMidSouth_2010_mat)){
  if(sum(latMidSouth_2010_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
latMidSouth_2010_mat <- latMidSouth_2010_mat[, -zeros]

nrow(latMidSouth_2010_mat) # 31 sites
ncol(latMidSouth_2010_mat) # 12 parasites 

#### METACOMMUNITY ANALYSES #####
library(metacom)

# Run the gamet of metacommunity analyses to determine structure
latMidSouth_2010_metacom <- metacommunity(latMidSouth_2010_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
latMidSouth_2010_Comm <- latMidSouth_2010_metacom$Comm
latMidSouth_2010_plot <- Matrix_Plot(latMidSouth_2010_Comm)
quartz(height=6, width=3)
print(latMidSouth_2010_plot) 
# looks complicated...at least two compartments?

# View Coherence
latMidSouth_2010_metacom$Coherence
# p-val = 0.17, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# latMidSouth by year: 2011

latMidSouth$site.yr <- as.character(latMidSouth$site.yr)
latMidSouth_2011 <- subset(latMidSouth, 
                           substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")

# Need to remove some columns
latMidSouth_2011 <- latMidSouth_2011[, -c(1, 18:19)]
head(latMidSouth_2011)

# Turn into a matrix
latMidSouth_2011_mat <- as.matrix(latMidSouth_2011)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(latMidSouth_2011_mat)){
  if(sum(latMidSouth_2011_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
latMidSouth_2011_mat <- latMidSouth_2011_mat[, -zeros]

nrow(latMidSouth_2011_mat) # 20 sites
ncol(latMidSouth_2011_mat) # 10 parasites 

#### METACOMMUNITY ANALYSES #####
library(metacom)

# Run the gamet of metacommunity analyses to determine structure
latMidSouth_2011_metacom <- metacommunity(latMidSouth_2011_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
latMidSouth_2011_Comm <- latMidSouth_2011_metacom$Comm
latMidSouth_2011_plot <- Matrix_Plot(latMidSouth_2011_Comm)
quartz(height=6, width=3)
print(latMidSouth_2011_plot) 
# looks nested

# View Coherence
latMidSouth_2011_metacom$Coherence
# p-val = 0.033, positive coherence

latMidSouth_2011_metacom$Turnover
# p-val = 0.071, borderline significant turnover (nested/quasi-nested)

latMidSouth_2011_metacom$Boundary
# 1.78, p-val < 0.0001

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# latMidSouth by year: 2012

latMidSouth$site.yr <- as.character(latMidSouth$site.yr)
latMidSouth_2012 <- subset(latMidSouth, 
                           substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2012")

# Need to remove some columns
latMidSouth_2012 <- latMidSouth_2012[, -c(1, 18:19)]
head(latMidSouth_2012)

# Turn into a matrix
latMidSouth_2012_mat <- as.matrix(latMidSouth_2012)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(latMidSouth_2012_mat)){
  if(sum(latMidSouth_2012_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
latMidSouth_2012_mat <- latMidSouth_2012_mat[, -zeros]

nrow(latMidSouth_2012_mat) # 19 sites
ncol(latMidSouth_2012_mat) # 10 parasites 

#### METACOMMUNITY ANALYSES #####
library(metacom)

# Run the gamet of metacommunity analyses to determine structure
latMidSouth_2012_metacom <- metacommunity(latMidSouth_2012_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
latMidSouth_2012_Comm <- latMidSouth_2012_metacom$Comm
latMidSouth_2012_plot <- Matrix_Plot(latMidSouth_2012_Comm)
quartz(height=6, width=3)
print(latMidSouth_2012_plot) 
# looks maybe nested

# View Coherence
latMidSouth_2012_metacom$Coherence
# p-val = 0.065, borderline positive coherence

latMidSouth_2012_metacom$Turnover
# p-val = 0.10, borderline significant turnover (nested/quasi-nested)

latMidSouth_2012_metacom$Boundary
# 1.90, p-val < 0.0001


