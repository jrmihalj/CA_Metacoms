# California parasite metacommunities
# Analysis of metacommunity structure at the Site(Pond)-level 
# Including only parasite species acquired within the pond
# Includes only sites with >10 PSRE dissected. 
# Oct 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Here is a list of parasite species that will be included in the study:
# Alar - Trematode
# Clin - Trematode
# Fib - Trematode
# Glob - Trematode
# Echi - Trematode
# Mano - Trematode
# Nyct - Protozoa
# Opal - Protist
# Pinw - Nematode
# Rib - Trematode
# Thic - Trematode

# Gorg, Glyp, Hae, Hali, Mega, Ozwa, Rhab and Spir are all thought to be 
# acquired post-metamorphosis.

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# I will use the yearly PSRE-level datasets to compile site-level data for each year.

# Build a function to construct site-level data from the individ-level data:

Site.Com <- function(dataset){
  data <- dataset[, c(1, 9:ncol(dataset))] #Just site.yr and parasites
  
  site.names <- unique(dataset[, 1]) # Store the unique site names
  parasite.names <- colnames(data)[2:ncol(data)] # Store the parasite names
  
  # Create a data frame to store parasite data
  site.df <- mat.or.vec(nr=length(site.names), nc=length(parasite.names)+1)
  site.df <- data.frame(site.df, stringsAsFactors=F)
  colnames(site.df) <- c("site.yr", parasite.names)
  
  # Sum the parasite columns and store
  for(i in 1:length(site.names)){
    which.site <- site.names[i]
    which.hosts <- which(data$site.yr==which.site) #Which rows have hosts from site[i]
    host.parasites <- data[which.hosts, 2:ncol(data)]
    
    parasite.sums <- NULL
    for(j in 1:ncol(host.parasites)){
      parasite.sums[j] <- sum(host.parasites[, j])
    }
    
    site.df[i ,] <- c(which.site, parasite.sums)  
  }
  
  site.df 
  # Because I want an incidence matrix, turn the sums >=1 to 1
  for(i in 1:nrow(site.df)){
    for(j in 2:ncol(site.df)){
      if(site.df[i, j] >=1) site.df[i, j] <- 1
    }   
  }
   
  # Remove columns (parasite) that sum to zero
  site.mat <- data.matrix(site.df[, 2:ncol(site.df)])
  zeros <- NULL
  for (i in 1:ncol(site.mat)){
    if(sum(site.mat[, i]) == 0){
      zeros <- c(zeros, i)
    }
  }
  if(sum(zeros) >= 1) site.mat <- site.mat[, -zeros]
  
  site.df <- data.frame(site.names, site.mat)
  
  # Remove columns (parasite) that sum to zero
  zeros <- NULL
  for (i in 1:nrow(site.df)){
    if(sum(site.df[i, 2:ncol(site.df) ]) == 0){
      zeros <- c(zeros, i)
    }
  }
  if(sum(zeros) >= 1) site.df <- site.df[-zeros, ]
  
  return(site.df)
}


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Sites: 2009 (compare 2009 to 2010)
# Create an amalgamated site-level data frame 
Sites_9 <- Site.Com(PSRE_9)

# Remove site names and store as a matrix for the metacommunity analyses
Sites_9_mat <- data.matrix(Sites_9[, 2:ncol(Sites_9)])

nrow(Sites_9_mat) #30 sites (PRPND008 removed because no parasites present)
ncol(Sites_9_mat) #9 parasites

###################################################
######## Metacommunity Analyses: 2009 #############
###################################################

library(metacom)
Sites_9_metacom <- metacommunity(Sites_9_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
Sites_9_Comm <- Sites_9_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(Sites_9_Comm, xlab="Parasite Species", ylab="Sites"))
# Two compartments?

# View coherence:
Sites_9_metacom$Coherence
# p-val=0.39

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# Look at all of 2009 sites
# Sites: 2009_all
PSRE_9_all <- subset(PSRE, 
                 substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2009") 

# Create an amalgamated site-level data frame 
Sites_9_all <- Site.Com(PSRE_9_all)

# Remove site names and store as a matrix for the metacommunity analyses
Sites_9_all_mat <- data.matrix(Sites_9_all[, 2:ncol(Sites_9_all)])

nrow(Sites_9_all_mat) #84 sites (2 sites removed - no parasites)
ncol(Sites_9_all_mat) #10 parasites

###################################################
####### Metacommunity Analyses: 2009_all ##########
###################################################

library(metacom)
Sites_9_all_metacom <- metacommunity(Sites_9_all_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
Sites_9_all_Comm <- Sites_9_all_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(Sites_9_all_Comm, xlab="Parasite Species", ylab="Sites"))
# A jumble...

# View coherence:
Sites_9_all_metacom$Coherence
# p-val=0.78
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Sites: 2010(a) (comparing 2010 to 2009)
# Create an amalgamated site-level data frame 
Sites_10a <- Site.Com(PSRE_10)

# Remove site names and store as a matrix for the metacommunity analyses
Sites_10a_mat <- data.matrix(Sites_10a[, 2:ncol(Sites_10a)])

nrow(Sites_10a_mat) #31 sites
ncol(Sites_10a_mat) #8 parasites

###################################################
######## Metacommunity Analyses: 2010(a) ##########
###################################################

library(metacom)
Sites_10a_metacom <- metacommunity(Sites_10a_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
Sites_10a_Comm <- Sites_10a_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(Sites_10a_Comm, xlab="Parasite Species", ylab="Sites"))
# Doesn't look like much at all...

# View coherence:
Sites_10a_metacom$Coherence
# p-val=0.75, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Sites: 2010(b) (comparing 2010 to 2011)
# Create an amalgamated site-level data frame 
Sites_10b <- Site.Com(PSRE_10b)

# Remove site names and store as a matrix for the metacommunity analyses
Sites_10b_mat <- data.matrix(Sites_10b[, 2:ncol(Sites_10b)])

nrow(Sites_10b_mat) #55 sites (PRPND004 removed - no parasites)
ncol(Sites_10b_mat) #9 parasites

###################################################
######## Metacommunity Analyses: 2010(b) ##########
###################################################

library(metacom)
Sites_10b_metacom <- metacommunity(Sites_10b_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
Sites_10b_Comm <- Sites_10b_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(Sites_10b_Comm, xlab="Parasite Species", ylab="Sites"))
# Looks nested (Sub-compartment with sites that contain Fib parasite)

# View coherence:
Sites_10b_metacom$Coherence
# p-val=0.017, positive coherence

# Turnover
Sites_10b_metacom$Turnover
# p-val=0.33, fewer replacements than mean, so quasi-nested

# Boundary Clumping
Sites_10b_metacom$Boundary
# index=1.92, p<0.0001, clumped losses

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# Because, 2010 showed different patterns with different subsets of data, include
# all data here:

# Sites: 2010_all
# Create an amalgamated site-level data frame 
Sites_10_all <- Site.Com(PSRE_10_all)

# Remove site names and store as a matrix for the metacommunity analyses
Sites_10_all_mat <- data.matrix(Sites_10_all[, 2:ncol(Sites_10_all)])

nrow(Sites_10_all_mat) #108 sites (PRPND004 removed for no parasites)
ncol(Sites_10_all_mat) #9 parasites

###################################################
####### Metacommunity Analyses: 2010_all ##########
###################################################

library(metacom)
Sites_10_all_metacom <- metacommunity(Sites_10_all_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
Sites_10_all_Comm <- Sites_10_all_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(Sites_10_all_Comm, xlab="Parasite Species", ylab="Sites"))
# Busy...sort of nested, but hard to say

# View coherence:
Sites_10_all_metacom$Coherence
# p-val=0.67, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# Sites: 2011_a (compared to 2010)
# Create an amalgamated site-level data frame 
Sites_11a <- Site.Com(PSRE_11)

# Remove site names and store as a matrix for the metacommunity analyses
Sites_11a_mat <- data.matrix(Sites_11a[, 2:ncol(Sites_11a)])

nrow(Sites_11a_mat) #56 sites (no sites removed)
ncol(Sites_11a_mat) #11 parasites

###################################################
####### Metacommunity Analyses: 2011_a ############
###################################################

library(metacom)
Sites_11a_metacom <- metacommunity(Sites_11a_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
Sites_11a_Comm <- Sites_11a_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(Sites_11a_Comm, xlab="Parasite Species", ylab="Sites"))
# Looks pretty nested

# View coherence:
Sites_11a_metacom$Coherence
# p-val=0.03, positive coherence

Sites_11a_metacom$Turnover
# p-val=0.57, fewer replacements, so quasi-nested

Sites_11a_metacom$Boundary
# index=1.65, p<0.0001, clumped boundaries

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# Sites: 2011_b (compared to 2012)
# Create an amalgamated site-level data frame 
Sites_11b <- Site.Com(PSRE_11b)

# Remove site names and store as a matrix for the metacommunity analyses
Sites_11b_mat <- data.matrix(Sites_11b[, 2:ncol(Sites_11b)])

nrow(Sites_11b_mat) #31 sites (no sites removed)
ncol(Sites_11b_mat) #10 parasites

###################################################
####### Metacommunity Analyses: 2011_b ############
###################################################

library(metacom)
Sites_11b_metacom <- metacommunity(Sites_11b_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
Sites_11b_Comm <- Sites_11b_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(Sites_11b_Comm, xlab="Parasite Species", ylab="Sites"))
# Looks somewhat nested

# View coherence:
Sites_11b_metacom$Coherence
# p-val=0.025, positive coherence

Sites_11b_metacom$Turnover
# p-val=0.42, fewer replacements, so quasi-nested

Sites_11b_metacom$Boundary
# index=1.83, p<0.0001, clumped boundaries

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# Check all of 2011, as subsets a and b both seem nested
# Sites: 2011_all
# Create an amalgamated site-level data frame 
Sites_11_all <- Site.Com(PSRE_11_all)

# Remove site names and store as a matrix for the metacommunity analyses
Sites_11_all_mat <- data.matrix(Sites_11_all[, 2:ncol(Sites_11_all)])

nrow(Sites_11_all_mat) #69 sites (no sites removed)
ncol(Sites_11_all_mat) #11 parasites

###################################################
####### Metacommunity Analyses: 2011_all ##########
###################################################

library(metacom)
Sites_11_all_metacom <- metacommunity(Sites_11_all_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
Sites_11_all_Comm <- Sites_11_all_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(Sites_11_all_Comm, xlab="Parasite Species", ylab="Sites"))
# Still looks somewhat nested

# View coherence:
Sites_11_all_metacom$Coherence
# p-val=0.029, positive coherence

Sites_11_all_metacom$Turnover
# p-val=0.92, fewer replacements, so quasi-nested

Sites_11_all_metacom$Boundary
# index=1.36, p<0.0001, clumped boundaries

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# Sites: 2012 (compared to 2011)
# Create an amalgamated site-level data frame 
Sites_12 <- Site.Com(PSRE_12)

# Remove site names and store as a matrix for the metacommunity analyses
Sites_12_mat <- data.matrix(Sites_12[, 2:ncol(Sites_12)])

nrow(Sites_12_mat) #28 sites (3 sites removed)
ncol(Sites_12_mat) #7 parasites

###################################################
####### Metacommunity Analyses: 2011_all ##########
###################################################

library(metacom)
Sites_12_metacom <- metacommunity(Sites_12_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
Sites_12_Comm <- Sites_12_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(Sites_12_Comm, xlab="Parasite Species", ylab="Sites"))
# Looks somewhat nested

# View coherence:
Sites_12_metacom$Coherence
# p-val=0.060, positive coherence

Sites_12_metacom$Turnover
# p-val=0.42, fewer replacements, so quasi-nested

Sites_12_metacom$Boundary
# index=1.76, p<0.0001, clumped boundaries
