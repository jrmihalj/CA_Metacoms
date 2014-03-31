# California amphibian HOST metacommunities

# Separated by year as well
# Includes only sites with >10 PSRE dissected. 
# August 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Import the data set that has all sites and pres/abs of all hosts
# 'Sites_Hosts.csv'
hosts <- read.csv(file.choose(), header=T)
head(hosts)

# Need to pick out the sites that are present in the parasite data set ('all_sites')
# Sites with >10 PSRE dissected, essentially

zeros <- NULL
for (i in 1:nrow(hosts)){
  out <- as.numeric(all_sites$site.yr %in% paste(hosts$site.yr[i]))
  if(sum(out)==0){
    zeros <- c(zeros, i)
  }
}
zeros
hosts <- hosts[-zeros ,]

# Remove first column, then convert to matrix
hosts_mat <- hosts[, -1]
hosts_mat <- as.matrix(hosts)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(hosts_mat)){
  if(sum(hosts_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # no zeros


# Check if any rows sums = 0 and remove
zeros <- NULL
for (j in 1:nrow(hosts_mat)){
  if(sum(hosts_mat[j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # NO zeros

nrow(hosts_mat) # 240 sites
ncol(hosts_mat) # 6 hosts

#### METACOMMUNITY ANALYSES #####
library(metacom)

# Run the gamet of metacommunity analyses to determine structure
hosts_metacom <- metacommunity(hosts_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
hosts_Comm <- hosts_metacom$Comm
hosts_plot <- Matrix_Plot(hosts_Comm, xlab="Host Species", ylab="Sites")
quartz(height=6, width=3)
print(hosts_plot) # looks like two clear compartments?

# View Coherence
hosts_metacom$Coherence
# p-val = 0.46, no coherence, but two compartments

#######################
# First compartment: Rows 1-36 of the first ordinated matrix

hosts_compA_metacom <- metacommunity(hosts_Comm[1:36 ,], order=F, method="r1", sims=1000)
# Store and view the ordinated occupancy matrix
hosts_compA_Comm <- hosts_compA_metacom$Comm
hosts_compA_plot <- Matrix_Plot(hosts_compA_Comm, xlab="Host Species", ylab="Sites")
quartz(height=6, width=3)
print(hosts_compA_plot) 

# View Coherence
hosts_compA_metacom$Coherence
# p-val = 0.31, no coherence

#######################
# Second compartment: Rows 37-240 of the first ordinated matrix

hosts_compB_metacom <- metacommunity(hosts_Comm[37:240 , 1:5], order=F, method="r1", sims=1000)
# Store and view the ordinated occupancy matrix
hosts_compB_Comm <- hosts_compB_metacom$Comm
hosts_compB_plot <- Matrix_Plot(hosts_compB_Comm, xlab="Host Species", ylab="Sites")
quartz(height=6, width=3)
print(hosts_compB_plot) # looks like two clear compBrtments?

# View Coherence
hosts_compB_metacom$Coherence
# p-val = 0.24, no coherence
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now subset the data, based on year:

# Year 2009
hosts$site.yr <- as.character(hosts$site.yr)
hosts_2009 <- subset(hosts, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2009")
head(hosts_2009)

# Need to remove headings and first column
hosts_2009_mat <- hosts_2009[, 2:7]
head(hosts_2009_mat)

# Turn into a matrix
hosts_2009_mat <- as.matrix(hosts_2009_mat)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(hosts_2009_mat)){
  if(sum(hosts_2009_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # No zeros

nrow(hosts_2009_mat) # 71 sites
ncol(hosts_2009_mat) # 6 hosts 

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
hosts_2009_metacom <- metacommunity(hosts_2009_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_hosts_2009 <- hosts_2009_metacom$Comm
ord_hosts_2009_plot <- Matrix_Plot(ordinated_hosts_2009)
quartz(height=6, width=3)
print(ord_hosts_2009_plot)
# Doesn't look like much

hosts_2009_metacom$Coherence
#p-val = 0.46, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Year 2010
hosts_2010 <- subset(hosts, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2010")
head(hosts_2010)

# Need to remove headings and first column
hosts_2010_mat <- hosts_2010[, 2:7]
head(hosts_2010_mat)

# Turn into a matrix
hosts_2010_mat <- as.matrix(hosts_2010_mat)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(hosts_2010_mat)){
  if(sum(hosts_2010_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # No zeros

nrow(hosts_2010_mat) # 87 sites
ncol(hosts_2010_mat) # 6 hosts 

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
hosts_2010_metacom <- metacommunity(hosts_2010_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_hosts_2010 <- hosts_2010_metacom$Comm
ord_hosts_2010_plot <- Matrix_Plot(ordinated_hosts_2010, xlab="Host Species", ylab="Sites")
quartz(height=6, width=3)
print(ord_hosts_2010_plot)
# Interesting, looks like two compartments, just like the 2010 host-parasite dataset

hosts_2010_metacom$Coherence
#p-val = 0.44, no coherence though...

###################
# First compartment
hosts_2010_compA <- metacommunity(ordinated_hosts_2010[1:14 ,], order=F, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
hosts_2010_compA_Comm <- hosts_2010_compA$Comm
hosts_2010_compA_plot <- Matrix_Plot(hosts_2010_compA_Comm, xlab="Host Species", ylab="Sites")
quartz(height=3, width=3)
print(hosts_2010_compA_plot)
# 

hosts_2010_compA$Coherence
# p-val = 0.090 borderline positive coherence

hosts_2010_compA$Turnover
# p-val = 0.025, significant, turnover (nested subsets)

hosts_2010_compA$Boundary
# 0.923, p=0.44, no clumped losses (or at least couldn't find significance)

###################
# Second compartment
hosts_2010_compB <- metacommunity(ordinated_hosts_2010[15:87 , 1:5], order=F, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
hosts_2010_compB_Comm <- hosts_2010_compB$Comm
hosts_2010_compB_plot <- Matrix_Plot(hosts_2010_compB_Comm, xlab="Host Species", ylab="Sites")
quartz(height=3, width=3)
print(hosts_2010_compB_plot)
# 

hosts_2010_compB$Coherence
# p-val = 0.10, borderline positive coherence

hosts_2010_compB$Turnover
# p-val = 0.015, significant turnover (nested subsets)

hosts_2010_compB$Boundary
# 1.22, p=0.001, clumped losses

################################
# Extract the ordination site-scores to correlate with host-level ordination. 
# Save the site scores of the ordinated matrix 
scores_amph_2010 <- OrderMatrix(ordinated_hosts_2010, output.scores=T)
sitescores_amph_2010 <- scores_amph_2010$sitescores
sitescores_amph_2010 <- melt(sitescores_amph_2010)
AmphSiteID <- rownames(sitescores_amph_2010)
sitescores_amph_2010 <- data.frame(AmphSiteID, sitescores_amph_2010)
head(sitescores_amph_2010)

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Year 2011
hosts_2011 <- subset(hosts, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")
head(hosts_2011)

# Need to remove headings and first column
hosts_2011_mat <- hosts_2011[, 2:7]
head(hosts_2011_mat)

# Turn into a matrix
hosts_2011_mat <- as.matrix(hosts_2011_mat)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(hosts_2011_mat)){
  if(sum(hosts_2011_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # No zeros

nrow(hosts_2011_mat) # 54 sites
ncol(hosts_2011_mat) # 6 hosts 

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
hosts_2011_metacom <- metacommunity(hosts_2011_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_hosts_2011 <- hosts_2011_metacom$Comm
ord_hosts_2011_plot <- Matrix_Plot(ordinated_hosts_2011)
quartz(height=6, width=3)
print(ord_hosts_2011_plot)
# Nothing happening here

hosts_2011_metacom$Coherence
#p-val = 0.42, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Year 2012
hosts_2012 <- subset(hosts, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2012")
head(hosts_2012)

# Need to remove headings and first column
hosts_2012_mat <- hosts_2012[, 2:7]
head(hosts_2012_mat)

# Turn into a matrix
hosts_2012_mat <- as.matrix(hosts_2012_mat)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(hosts_2012_mat)){
  if(sum(hosts_2012_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros # No zeros

nrow(hosts_2012_mat) # 28 sites
ncol(hosts_2012_mat) # 6 hosts 

#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
hosts_2012_metacom <- metacommunity(hosts_2012_mat, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_hosts_2012 <- hosts_2012_metacom$Comm
ord_hosts_2012_plot <- Matrix_Plot(ordinated_hosts_2012, xlab="Host Species", ylab="Sites")
quartz(height=6, width=3)
print(ord_hosts_2012_plot)
# maybe something happening here
# Two compartments

hosts_2012_metacom$Coherence
#p-val = 0.025, positive coherence

hosts_2012_metacom$Turnover
# p-val = 0.14, more replacements than expected sort of

hosts_2012_metacom$Boundary
# 1.03, p=0.24, no clumping 

# So quasi-Gleasonian structure...

################################
# Extract the ordination site-scores to correlate with host-level ordination. 
# Save the site scores of the ordinated matrix 
scores_amph_2012 <- OrderMatrix(ordinated_hosts_2012, output.scores=T)
sitescores_amph_2012 <- scores_amph_2012$sitescores
sitescores_amph_2012 <- melt(sitescores_amph_2012)
AmphSiteID <- rownames(sitescores_amph_2012)
sitescores_amph_2012 <- data.frame(AmphSiteID, sitescores_amph_2012)
head(sitescores_amph_2012)




