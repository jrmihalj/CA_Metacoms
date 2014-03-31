# California snail host metacommunities

# Separated by year as well
# Includes only sites with >10 PSRE dissected. 
# August 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Import the data set that has all sites and pres/abs of all hosts
# 'Sites_Snails.csv'
snails <- read.csv(file.choose(), header=T)
head(snails)

# Need to pick out the sites that are present in the parasite data set
# Sites with >10 PSRE dissected, essentially

zeros <- NULL
for (i in 1:nrow(snails)){
  out <- as.numeric(all_sites$site.yr %in% paste(snails$site.yr[i]))
  if(sum(out)==0){
    zeros <- c(zeros, i)
  }
}
zeros
snails <- snails[-zeros ,]

# Remove first column, then convert to matrix
snails_mat <- snails[, -1]
snails_mat <- as.matrix(snails_mat)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(snails_mat)){
  if(sum(snails_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # no zeros


# Check if any rows sums = 0 and remove
zeros <- NULL
for (j in 1:nrow(snails_mat)){
  if(sum(snails_mat[j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros 
snails_mat <- snails_mat[-zeros ,]

nrow(snails_mat) # 227 sites
ncol(snails_mat) # 5 snails

#### METACOMMUNITY ANALYSES #####
library(metacom)

# Run the gamet of metacommunity analyses to determine structure
snails_metacom <- metacommunity(snails_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
snails_Comm <- snails_metacom$Comm
snails_plot <- Matrix_Plot(snails_Comm, xlab="Snail Species", ylab="Sites")
quartz(height=6, width=3)
print(snails_plot) # looks like two clear compartments?

# View Coherence
snails_metacom$Coherence
# p-val = 0.032, positive coherence

snails_metacom$Turnover
# p-val = 0.23, less than expected, so quasi-nested

snails_metacom$Boundary
# index=1.34, p < 0.0001, clumped losses

#######################
# First compartment: Rows 1-51 of the first ordinated matrix

snails_compA_metacom <- metacommunity(snails_Comm[1:51 ,], order=T, method="r1", sims=1000)
# Store and view the ordinated occupancy matrix
snails_compA_Comm <- snails_compA_metacom$Comm
snails_compA_plot <- Matrix_Plot(snails_compA_Comm, xlab="Snail Species", ylab="Sites")
quartz(height=3, width=3)
print(snails_compA_plot) # looks like two clear compartments?

# View Coherence
snails_compA_metacom$Coherence
# p-val < 0.0001, positive coherence

snails_compA_metacom$Turnover
# p-val = 0.064, borderline. quasi-nested

snails_compA_metacom$Boundary
# index=1.17, p=0.010, no clumping

#######################
# Second compartment: Rows 52-227 of the first ordinated matrix

snails_compB_metacom <- metacommunity(snails_Comm[52:227 , 1:3], order=T, method="r1", sims=1000)
# Store and view the ordinated occupancy matrix
snails_compB_Comm <- snails_compB_metacom$Comm
snails_compB_plot <- Matrix_Plot(snails_compB_Comm, xlab="Snail Species", ylab="Sites")
quartz(height=3, width=3)
print(snails_compB_plot) # looks nested

# View Coherence
snails_compB_metacom$Coherence
# p-val < 0.0001, positive coherence
snails_compB_metacom$Turnover
# p-val < 0.0001, nested
snails_compB_metacom$Boundary
# index=1, no clumping

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now subset the data, based on year:

# Year 2009
snails$site.yr <- as.character(snails$site.yr)
snails_2009 <- subset(snails, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2009")
head(snails_2009)

# Need to remove headings and first column
snails_2009_mat <- snails_2009[, 2:6]
head(snails_2009_mat)

# Turn into a matrix
snails_2009_mat <- as.matrix(snails_2009_mat)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(snails_2009_mat)){
  if(sum(snails_2009_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # No zeros

# Check rows
zeros <- NULL
for (j in 1:nrow(snails_2009_mat)){
  if(sum(snails_2009_mat[j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros 
snails_2009_mat <- snails_2009_mat[-zeros ,]


nrow(snails_2009_mat) # 69 sites
ncol(snails_2009_mat) # 5 snails 

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
snails_2009_metacom <- metacommunity(snails_2009_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_snails_2009 <- snails_2009_metacom$Comm
ord_snails_2009_plot <- Matrix_Plot(ordinated_snails_2009)
quartz(height=6, width=3)
print(ord_snails_2009_plot)
# Doesn't look like much

snails_2009_metacom$Coherence
#p-val = 0.65, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Year 2010
snails_2010 <- subset(snails, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2010")
head(snails_2010)

# Need to remove headings and first column
snails_2010_mat <- snails_2010[, 2:6]
head(snails_2010_mat)

# Turn into a matrix
snails_2010_mat <- as.matrix(snails_2010_mat)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(snails_2010_mat)){
  if(sum(snails_2010_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # No zeros

# Check rows
zeros <- NULL
for (j in 1:nrow(snails_2010_mat)){
  if(sum(snails_2010_mat[j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros 
snails_2010_mat <- snails_2010_mat[-zeros ,]

nrow(snails_2010_mat) # 79 sites
ncol(snails_2010_mat) # 5 snails 

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
snails_2010_metacom <- metacommunity(snails_2010_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_snails_2010 <- snails_2010_metacom$Comm
ord_snails_2010_plot <- Matrix_Plot(ordinated_snails_2010, xlab="Snail Species", ylab="Sites")
quartz(height=6, width=3)
print(ord_snails_2010_plot)
# Looks nested?

snails_2010_metacom$Coherence
#p-val = 0.0023, positive coherence

snails_2010_metacom$Turnover
#p-val = 0.0015, nested

snails_2010_metacom$Boundary
#index=1.53, p < 0.0001, clumped losses

################################
# Extract the ordination site-scores to correlate with host-level ordination. 
# Save the site scores of the ordinated matrix 
scores_snails_2010 <- OrderMatrix(ordinated_snails_2010, output.scores=T)
sitescores_snails_2010 <- scores_snails_2010$sitescores
sitescores_snails_2010 <- melt(sitescores_snails_2010)
snailsSiteID <- rownames(sitescores_snails_2010)
sitescores_snails_2010 <- data.frame(snailsSiteID, sitescores_snails_2010)
head(sitescores_snails_2010)


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Year 2011
snails_2011 <- subset(snails, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")
head(snails_2011)

# Need to remove headings and first column
snails_2011_mat <- snails_2011[, 2:6]
head(snails_2011_mat)

# Turn into a matrix
snails_2011_mat <- as.matrix(snails_2011_mat)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(snails_2011_mat)){
  if(sum(snails_2011_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
snails_2011_mat <- snails_2011_mat[, -zeros]

# Check rows
zeros <- NULL
for (j in 1:nrow(snails_2011_mat)){
  if(sum(snails_2011_mat[j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros 
snails_2011_mat <- snails_2011_mat[-zeros ,]

nrow(snails_2011_mat) # 51 sites
ncol(snails_2011_mat) # 4 snails 

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
snails_2011_metacom <- metacommunity(snails_2011_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_snails_2011 <- snails_2011_metacom$Comm
ord_snails_2011_plot <- Matrix_Plot(ordinated_snails_2011)
quartz(height=6, width=3)
print(ord_snails_2011_plot)
# two nested compartments?

snails_2011_metacom$Coherence
#p-val = 0.05, negative coherence? checkerboard? 

#################################
##### THE PROGRAM ISN"T WORKING FOR THESE COMPARTMENTS?
# First compartment: rows 1:10
snails_compA2011_metacom <- metacommunity(ordinated_snails_2011[1:10 ,], order=F, method="r1", sims=1000)
# Store and view the ordinated occupancy matrix
snails_compA2011_Comm <- snails_compA2011_metacom$Comm
snails_compA2011_plot <- Matrix_Plot(snails_compA2011_Comm, xlab="Host Species", ylab="Sites")
quartz(height=6, width=3)
print(snails_compA2011_plot) # looks like two clear compartments?

# View Coherence
snails_compA2011_metacom$Coherence


#################################
# Second compartment: rows 11:51
snails_compB2011_metacom <- metacommunity(ordinated_snails_2011[11:51 ,-4], order=F, method="r1", sims=1000)
# Store and view the ordinated occupancy matrix
snails_compB2011_Comm <- snails_compB2011_metacom$Comm
snails_compB2011_plot <- Matrix_Plot(snails_compB2011_Comm, xlab="Host Species", ylab="Sites")
quartz(height=6, width=3)
print(snails_compB2011_plot) # looks like two clear compBrtments?

# View Coherence
snails_compB2011_metacom$Coherence
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Year 2012
snails_2012 <- subset(snails, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2012")
head(snails_2012)

# Need to remove headings and first column
snails_2012_mat <- snails_2012[, 2:6]
head(snails_2012_mat)

# Turn into a matrix
snails_2012_mat <- as.matrix(snails_2012_mat)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(snails_2012_mat)){
  if(sum(snails_2012_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # No zeros

# Check rows
zeros <- NULL
for (j in 1:nrow(snails_2012_mat)){
  if(sum(snails_2012_mat[j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # No zeros

nrow(snails_2012_mat) # 28 sites
ncol(snails_2012_mat) # 5 snails 

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
snails_2012_metacom <- metacommunity(snails_2012_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_snails_2012 <- snails_2012_metacom$Comm
ord_snails_2012_plot <- Matrix_Plot(ordinated_snails_2012, xlab="Snail Species", ylab="Sites")
quartz(height=6, width=3)
print(ord_snails_2012_plot)
# Looks nested

snails_2012_metacom$Coherence
#p-val = 0.0032, positive coherence

snails_2012_metacom$Turnover
# p-val = 0.0016, nested subsets

snails_2012_metacom$Boundary
# 1.51, p=0.0018, clumped losses

################################
# Extract the ordination site-scores to correlate with host-level ordination. 
# Save the site scores of the ordinated matrix 
scores_snails_2012 <- OrderMatrix(ordinated_snails_2012, output.scores=T)
sitescores_snails_2012 <- scores_snails_2012$sitescores
sitescores_snails_2012 <- melt(sitescores_snails_2012)
snailsSiteID <- rownames(sitescores_snails_2012)
sitescores_snails_2012 <- data.frame(snailsSiteID, sitescores_snails_2012)
head(sitescores_snails_2012)

