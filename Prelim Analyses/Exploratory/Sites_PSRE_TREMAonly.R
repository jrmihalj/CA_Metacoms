# California parasite metacommunities
# Analysis of metacommunity structure at the site (pond)-level 
# Including only trematode species
# Includes only sites with >10 PSRE dissected. 
# August 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Choose the file: "Sites_ALL_PSRE_With headers.csv", which contains headers and 
# site-yrs
all_sites <- read.csv(file.choose(), header=T)
head(all_sites)

#Remove the four non-trematode species:
all_trema <- all_sites[, -c(12:14, 17)]
head(all_trema)

#Turn the data.frame into an occupancy matrix:
#Remove first column (site.yr) and other column names
all_trema_mat <- all_trema[, -1]
all_trema_mat <- as.matrix(all_trema_mat)
# Remove rows with sum = 0
zeros <- NULL
for (j in 1:nrow(all_trema_mat)){
  if(sum(all_trema_mat[j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros

all_trema_mat <- all_trema_mat[-zeros ,]

colnames(all_trema_mat) <- NULL
head(all_trema_mat)

##### METACOMMUNITY ANALYSES ####
# for all PSRE sites, only trematodes

library(metacom)
ordinated_all_trema <- OrderMatrix(all_trema_mat)
image(apply(ordinated_all_trema,2,rev), col=c(0,'black'), axes=FALSE, xlab='species', ylab='sites', cex.lab=1.4)
box(lwd=2,col=1)

#Determine coherence. Must be significant to continue with analyses. Otherwise
#the metacommunity is assembled randomly
site_trema_coherence <- Coherence(all_trema_mat, 
                                sims=1000, method="r1", allow.empty=F)
site_trema_coherence$pval
# p-val = 0.12, so no coherence in the full dataset. No more analyses here. 

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for sites with 
# >10 PSRE dissected, only trematodes

# Start with 2009 sampled sites:

trema_2009 <- subset(all_trema, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2009")
head(trema_2009)
length(trema_2009$site.yr)
# 71 sites in 2009

# Need to remove headings and first column
trema_2009 <- trema_2009[, -1]
head(trema_2009)

# Turn into a matrix
trema_2009 <- as.matrix(trema_2009)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(trema_2009)){
  if(sum(trema_2009[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
trema_2009 <- trema_2009[, -zeros]
# Seems like some rows might have zeros as well?
zeros <- NULL
for (j in 1:nrow(trema_2009)){
  if(sum(trema_2009[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
trema_2009 <- trema_2009[-zeros ,]
nrow(trema_2009) # 70 sites
ncol(trema_2009) # 9 parasites


#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
trema_2009_metacom <- metacommunity(trema_2009, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_trema_2009 <- trema_2009_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(ordinated_trema_2009))
# Matrix looks nested

trema_2009_metacom$Coherence
#p-val = 0.0048, significant
#Fewer embedded absences than expected, so positive coherence

trema_2009_metacom$Turnover
#p-val= 0.033
#Fewer replacements than expected, so nested subsets

trema_2009_metacom$Boundary
#index=1.89, p<0.0001
#Clumped species loss

# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for sites with 
# >10 PSRE dissected, only trematodes

# Start with 2010 sampled sites:

trema_2010 <- subset(all_trema, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2010")
head(trema_2010)


# Need to remove headings and first column
trema_2010 <- trema_2010[, -1]
head(trema_2010)

# Turn into a matrix
trema_2010 <- as.matrix(trema_2010)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(trema_2010)){
  if(sum(trema_2010[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
trema_2010 <- trema_2010[, -zeros]
# Seems like some rows might have zeros as well?
zeros <- NULL
for (j in 1:nrow(trema_2010)){
  if(sum(trema_2010[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
trema_2010 <- trema_2010[-zeros ,]

nrow(trema_2010) # 80 sites
ncol(trema_2010) # 10 parasites


#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
trema_2010_metacom <- metacommunity(trema_2010, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_trema_2010 <- trema_2010_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(ordinated_trema_2010))
# Matrix looks random

trema_2010_metacom$Coherence
#p-val = 0.44, not significant coherence

# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for sites with 
# >10 PSRE dissected, only trematodes

# Start with 2011 sampled sites:

trema_2011 <- subset(all_trema, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")
head(trema_2011)


# Need to remove headings and first column
trema_2011 <- trema_2011[, -1]
head(trema_2011)

# Turn into a matrix
trema_2011 <- as.matrix(trema_2011)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(trema_2011)){
  if(sum(trema_2011[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
trema_2011 <- trema_2011[, -zeros]
# Seems like some rows might have zeros as well?
zeros <- NULL
for (j in 1:nrow(trema_2011)){
  if(sum(trema_2011[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
trema_2011 <- trema_2011[-zeros ,]

nrow(trema_2011) # 51 sites
ncol(trema_2011) # 7 parasites

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
trema_2011_metacom <- metacommunity(trema_2011, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_trema_2011 <- trema_2011_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(ordinated_trema_2011))
# Matrix looks random

trema_2011_metacom$Coherence
#p-val = 0.79, not significant coherence

# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for sites with 
# >10 PSRE dissected, only trematodes

# Start with 2012 sampled sites:

trema_2012 <- subset(all_trema, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2012")
head(trema_2012)


# Need to remove headings and first column
trema_2012 <- trema_2012[, -1]
head(trema_2012)

# Turn into a matrix
trema_2012 <- as.matrix(trema_2012)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(trema_2012)){
  if(sum(trema_2012[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
trema_2012 <- trema_2012[, -zeros]
# Seems like some rows might have zeros as well?
zeros <- NULL
for (j in 1:nrow(trema_2012)){
  if(sum(trema_2012[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # No zeros

nrow(trema_2012) # 28 sites
ncol(trema_2012) # 9 parasites

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
trema_2012_metacom <- metacommunity(trema_2012, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_trema_2012 <- trema_2012_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(ordinated_trema_2012))
# Matrix looks nested

trema_2012_metacom$Coherence
#p-val = 0.0078, positive coherence

trema_2012_metacom$Turnover
#p-val = 0.095, but fewer replacements than mean (quasi-nested)

trema_2012_metacom$Boundary
#index=1.74, p<0.0001, clumped losses