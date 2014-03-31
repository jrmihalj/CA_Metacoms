# California parasite metacommunities
# Analysis of metacommunity structure at the site (pond)-level 
# Including all parasite species
# Includes only sites with >10 PSRE dissected. 
# August 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# First, check the coherence of all PSRE sites, where the number of dissected 
# individuals is greater than 10

#Choose the correct file: "Sites_ALL_PSRE.csv"
site_all_mat <- read.csv(file.choose(), header=F)
#Turn the data.frame into an occupancy matrix:
site_all_mat <- as.matrix(site_all_mat)


library(metacom)
ordinated_all <- OrderMatrix(site_all_mat)
image(apply(ordinated_all,2,rev), col=c(0,'black'), axes=FALSE, xlab='species', ylab='sites', cex.lab=1.4)
box(lwd=2,col=1)

#Determine coherence. Must be significant to continue with analyses. Otherwise
#the metacommunity is assembled randomly
site_all_coherence <- Coherence(site_all_mat, 
                                sims=1000, method="r1", allow.empty=F)
site_all_coherence$pval
# 0.14 Means not significantly coherent. No further analyses for this level of
# the dataset 

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for sites with 
# >10 PSRE dissected

# Start with 2009 sampled sites
# I need to subset the dataframe to include only sites sampled from 2009
# This is difficult because the "site-yr" variable is concatinated,
# such as: Sitename_2009. So I need to subset in a tricky manner.

# Choose the file: "Sites_ALL_PSRE_With headers.csv", which contains headers and 
# site-yrs
all_sites <- read.csv(file.choose(), header=T)
head(all_sites)
#Use the function substr() to search for the last 4 characters of the "site.yr"
#variable. Search for ones that have "2009"

# site.yr needs to be a character string
all_sites$site.yr <- as.character(all_sites$site.yr)

# test the method:
site_names <- all_sites$site.yr
is.character(site_names)
# Look for "2009" in the last 4 characters of the site name
substr(site_names, nchar(site_names)-3, nchar(site_names)) %in% "2009"

# Now subset the data, based on this search:
sites_2009 <- subset(all_sites, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2009")
head(sites_2009)

# 71 sites in 2009

# Need to remove headings and first column
sites_2009 <- sites_2009[, 2:17]
head(sites_2009)

# Turn into a matrix
sites_2009 <- as.matrix(sites_2009)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(sites_2009)){
  if(sum(sites_2009[, j]) == 0){
   zeros <- c(zeros, j)
  }
}
zeros
sites_2009 <- sites_2009[, -zeros]

nrow(sites_2009) # 71 sites
ncol(sites_2009) # 12 parasites

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
site_2009_metacom <- metacommunity(sites_2009, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_2009 <- site_2009_metacom$Comm
ord2009_plot <- Matrix_Plot(ordinated_2009)
quartz(height=6, width=3)
print(ord2009_plot)
# Matrix looks nested, somewhat

site_2009_metacom$Coherence
#p-val = 0.066, no coherence (random structure)
# marginally significant

site_2009_metacom$Turnover
# p-val = 0.052, nested

site_2009_metacom$Boundary
#index = 2.10, p < 0.0001 clumped losses

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now for sites sampled in 2010

sites_2010 <- subset(all_sites, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2010")
head(sites_2010)
nrow(sites_2010)
# 87 sites in 2010

# Need to remove headings and first column
sites_2010 <- sites_2010[, 2:17]
head(sites_2010)

# Turn into a matrix
sites_2010 <- as.matrix(sites_2010)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(sites_2010)){
  if(sum(sites_2010[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
sites_2010 <- sites_2010[, -zeros]

nrow(sites_2010) # 87 sites
ncol(sites_2010) # 13 parasites

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
site_2010_metacom <- metacommunity(sites_2010, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_2010 <- site_2010_metacom$Comm
ord2010_plot <- Matrix_Plot(ordinated_2010)
quartz(height=3, width=3)
print(ord2010_plot)
# Doesn't seem like much going on

site_2010_metacom$Coherence
#p-val = 0.14, no coherence (random structure)

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now for sites sampled in 2011

sites_2011 <- subset(all_sites, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")
head(sites_2011)
nrow(sites_2011)
# 54 sites in 2011

# Need to remove headings and first column
sites_2011 <- sites_2011[, 2:17]
head(sites_2011)

# Turn into a matrix
sites_2011 <- as.matrix(sites_2011)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(sites_2011)){
  if(sum(sites_2011[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
sites_2011 <- sites_2011[, -zeros]

nrow(sites_2011) # 54 sites
ncol(sites_2011) # 11 parasites

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
site_2011_metacom <- metacommunity(sites_2011, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_2011 <- site_2011_metacom$Comm
ord2011_plot <- Matrix_Plot(ordinated_2011)
quartz(height=3, width=3)
print(ord2011_plot)
# Doesn't seem like much going on

site_2011_metacom$Coherence
#p-val = 0.85, no coherence (random structure)

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now for sites sampled in 2012

sites_2012 <- subset(all_sites, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2012")
head(sites_2012)
nrow(sites_2012)
# 28 sites in 2012

# Need to remove headings and first column
sites_2012 <- sites_2012[, 2:17]
head(sites_2012)

# Turn into a matrix
sites_2012 <- as.matrix(sites_2012)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(sites_2012)){
  if(sum(sites_2012[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
sites_2012 <- sites_2012[, -zeros]

nrow(sites_2012) # 28 sites
ncol(sites_2012) # 10 parasites

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
site_2012_metacom <- metacommunity(sites_2012, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_2012 <- site_2012_metacom$Comm
ord2012_plot <- Matrix_Plot(ordinated_2012)
quartz(height=3, width=3)
print(ord2012_plot)
# Looks somewhat nested, but depauperate

site_2012_metacom$Coherence
#p-val = 0.18, no coherence (random structure)