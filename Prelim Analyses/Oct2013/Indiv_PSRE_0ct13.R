# California parasite metacommunities
# Analysis of metacommunity structure at the inidividual (PSRE)-level 
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

# I'm going to start from the PSRE individuals and then scale up to the site-level.
# I will remove individuals that contain any of the parasites that are acquired post-
# metamorphosis and later determine the site-level pres/abs of each parasite. 

# Import the dataset that has all individual PSRE dissected from sites with >=10 PSRE dissected
# file name: 'Indiv_ALL_Host.csv'
PSRE <- read.csv(file.choose(), header=T)
head(PSRE)
# I need to subset to only have PSRE

PSRE <- subset(PSRE, Spp=="PSRE")
nrow(PSRE) # 3086 PSRE individuals

# Determine which rows (i.e. individuals) have any of the non-pond parasites
# These parasites are in columns:
rem_cols <- c(13:16, 19, 22, 24, 26)
remove <- NULL
for(i in 1:nrow(PSRE)){
  for(j in rem_cols){
    
    if( PSRE[i, j]==1){
      
      remove <- c(remove, i)
      
    }
  }
}
remove # 20 individuals

# Remove the identified individuals
PSRE <- PSRE[-remove, ]
nrow(PSRE) # 3066 individuals

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# 2009 vs. 2010 analysis
# Now I will pull out individuals from sites that were sampled in both 2009 and 2010.

# Need to make 'site.yr' a character string:
PSRE$site.yr <- as.character(PSRE$site.yr)

# Create two data frames that have either sites sampled in 2009 or 2010
PSRE_9 <- subset(PSRE, 
                 substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2009")
PSRE_10 <- subset(PSRE, 
                  substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2010")

# Now figure out which individuals were from sites sampled in 2010, that were also sampled 
# in 2009
remove <- NULL
for(i in 1:nrow(PSRE_10)){
  eval <- NULL
  eval <- substr(PSRE_9$site.yr, start=1, stop = nchar(PSRE_9$site.yr)-3) %in%
    substr(PSRE_10$site.yr[i], start=1, stop = nchar(PSRE_10$site.yr[i])-3)
  eval <- as.numeric(eval)
  if(sum(eval)==0){
    remove <- c(remove, i)
  }
}
length(remove)
PSRE_10 <- PSRE_10[-remove, ]

# Now see which 2009 match up to the 2010:

remove <- NULL
for(i in 1:nrow(PSRE_9)){
  eval <- NULL
  eval <- substr(PSRE_10$site.yr, start=1, stop = nchar(PSRE_10$site.yr)-3) %in%
    substr(PSRE_9$site.yr[i], start=1, stop = nchar(PSRE_9$site.yr[i])-3)
  eval <- as.numeric(eval)
  if(sum(eval)==0){
    remove <- c(remove, i)
  }
}
length(remove)
PSRE_9 <- PSRE_9[-remove, ]

nrow(PSRE_9) # 406 indiv
nrow(PSRE_10) # 309 indiv

# Make sure the same site.yr's are present in each data frame
unique(PSRE_9$site.yr)
unique(PSRE_10$site.yr) #31 sites each

###################################################
######## Metacommunity Analyses: 2009 #############
###################################################

# Remove unneccessary columns, convert to matrix
PSRE_9_mat <- as.matrix(PSRE_9[, 9:27])
head(PSRE_9_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_9_mat)){
  if(sum(PSRE_9_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_9_mat <- PSRE_9_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_9_mat)){
  if(sum(PSRE_9_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_9_mat <- PSRE_9_mat[, -zeros]

nrow(PSRE_9_mat) # 348 indiv.
ncol(PSRE_9_mat) # 9 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_9_metacom <- metacommunity(PSRE_9_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_9_Comm <- PSRE_9_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_9_Comm, xlab="Parasite Species", ylab="Host Indiv."))
# Maybe nested, but hard to say

# View coherence:
PSRE_9_metacom$Coherence
# pval = 0.44, no coherence

###################################################
######## Metacommunity Analyses: 2010 #############
###################################################

# Remove unneccessary columns, convert to matrix
PSRE_10_mat <- as.matrix(PSRE_10[, 9:27])
head(PSRE_10_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_10_mat)){
  if(sum(PSRE_10_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_10_mat <- PSRE_10_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_10_mat)){
  if(sum(PSRE_10_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_10_mat <- PSRE_10_mat[, -zeros]

nrow(PSRE_10_mat) # 302 indiv.
ncol(PSRE_10_mat) # 8 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_10_metacom <- metacommunity(PSRE_10_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_10_Comm <- PSRE_10_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_10_Comm, xlab="Parasite Species", ylab="Host Indiv."))
# Two clear compartments

# View coherence:
PSRE_10_metacom$Coherence
# pval = 0.22, no coherence

#########################
# Compartment A
# rows 1-175

PSRE_10_compA <- PSRE_10_Comm[1:175, ]

PSRE_10_compA_metacom <- metacommunity(PSRE_10_compA, order=F, method="r1", sims=1000)

PSRE_10_compA_Comm <- PSRE_10_compA_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_10_compA_Comm, xlab="Parasite Species", ylab="Host Indiv."))

PSRE_10_compA_metacom$Coherence
# p-val=0.96, no coherence

#########################
# Compartment B
# rows 176-302, columns 

PSRE_10_compB <- PSRE_10_Comm[176:302, 1:5]

PSRE_10_compB_metacom <- metacommunity(PSRE_10_compB, order=F, method="r1", sims=1000)

PSRE_10_compB_Comm <- PSRE_10_compB_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_10_compB_Comm, xlab="Parasite Species", ylab="Host Indiv."))

PSRE_10_compB_metacom$Coherence
# p-val = 0.087

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# 2010 vs. 2011 analysis
# Now I will pull out individuals from sites that were sampled in both 2009 and 2010.

# Create two data frames that have either sites sampled in 2009 or 2010
PSRE_10b <- subset(PSRE, 
                 substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2010")
PSRE_11 <- subset(PSRE, 
                  substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")

# Now figure out which individuals were from sites sampled in 2011, that were also sampled 
# in 2010
remove <- NULL
for(i in 1:nrow(PSRE_11)){
  eval <- NULL
  eval <- substr(PSRE_10b$site.yr, start=1, stop = nchar(PSRE_10b$site.yr)-3) %in%
    substr(PSRE_11$site.yr[i], start=1, stop = nchar(PSRE_11$site.yr[i])-3)
  eval <- as.numeric(eval)
  if(sum(eval)==0){
    remove <- c(remove, i)
  }
}
length(remove)
PSRE_11 <- PSRE_11[-remove, ]

# Now see which 2010 match up to the 2011:

remove <- NULL
for(i in 1:nrow(PSRE_10b)){
  eval <- NULL
  eval <- substr(PSRE_11$site.yr, start=1, stop = nchar(PSRE_11$site.yr)-3) %in%
    substr(PSRE_10b$site.yr[i], start=1, stop = nchar(PSRE_10b$site.yr[i])-3)
  eval <- as.numeric(eval)
  if(sum(eval)==0){
    remove <- c(remove, i)
  }
}
length(remove)
PSRE_10b <- PSRE_10b[-remove, ]

nrow(PSRE_10b) # 546 indiv
nrow(PSRE_11) # 534 indiv

# Make sure the same site.yr's are present in each data frame
unique(PSRE_10b$site.yr)
unique(PSRE_11$site.yr) # 56 sites each

###################################################
######## Metacommunity Analyses: 2010b #############
###################################################

# Remove unneccessary columns, convert to matrix
PSRE_10b_mat <- as.matrix(PSRE_10b[, 9:27])
head(PSRE_10b_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_10b_mat)){
  if(sum(PSRE_10b_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_10b_mat <- PSRE_10b_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_10b_mat)){
  if(sum(PSRE_10b_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_10b_mat <- PSRE_10b_mat[, -zeros]

nrow(PSRE_10b_mat) # 496 indiv.
ncol(PSRE_10b_mat) # 9 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_10b_metacom <- metacommunity(PSRE_10b_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_10b_Comm <- PSRE_10b_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_10b_Comm, xlab="Parasite Species", ylab="Host Indiv."))
# Again, two clear compartments.

# View coherence:
PSRE_10b_metacom$Coherence
# pval = 0.87, no coherence

#########################
# Compartment A
# rows 1-272

PSRE_10b_compA <- PSRE_10b_Comm[1:272, -1]

PSRE_10b_compA_metacom <- metacommunity(PSRE_10b_compA, order=F, method="r1", sims=1000)

PSRE_10b_compA_Comm <- PSRE_10b_compA_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_10b_compA_Comm, xlab="Parasite Species", ylab="Host Indiv."))

PSRE_10b_compA_metacom$Coherence
# p-val=0.43, no coherence

#########################
# Compartment B
# rows 272-496, columns 1-7

PSRE_10b_compB <- PSRE_10b_Comm[272:496, 1:7]

PSRE_10b_compB_metacom <- metacommunity(PSRE_10b_compB, order=F, method="r1", sims=1000)

PSRE_10b_compB_Comm <- PSRE_10b_compB_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_10b_compB_Comm, xlab="Parasite Species", ylab="Host Indiv."))

PSRE_10b_compB_metacom$Coherence
# p-val = 0.083

###########################
# Because 2010 and 2010b showed similar patterns, I'm going to use the entire 2010 dataset

PSRE_10_all <- subset(PSRE, 
                   substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2010")

# Remove unneccessary columns, convert to matrix
PSRE_10_all_mat <- as.matrix(PSRE_10_all[, 9:27])
head(PSRE_10_all_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_10_all_mat)){
  if(sum(PSRE_10_all_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_10_all_mat <- PSRE_10_all_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_10_all_mat)){
  if(sum(PSRE_10_all_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_10_all_mat <- PSRE_10_all_mat[, -zeros]

nrow(PSRE_10_all_mat) # 924 indiv.
ncol(PSRE_10_all_mat) # 9 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_10_all_metacom <- metacommunity(PSRE_10_all_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_10_all_Comm <- PSRE_10_all_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_10_all_Comm, xlab="Parasite Species", ylab="Host Indiv."))
# Again, two clear compartments.

# View coherence:
PSRE_10_all_metacom$Coherence
# pval = 0.40, no coherence, but two compartments

#########################
# Compartment A
# rows 1-403

PSRE_10_all_compA <- PSRE_10_all_Comm[1:272, -c(1:2)]

PSRE_10_all_compA_metacom <- metacommunity(PSRE_10_all_compA, order=T, method="r1", sims=1000)

PSRE_10_all_compA_Comm <- PSRE_10_all_compA_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_10_all_compA_Comm, xlab="Parasite Species", ylab="Host Indiv."))

PSRE_10_all_compA_metacom$Coherence
# p-val=0.037, positive coherence

PSRE_10_all_compA_metacom$Turnover
# p-val=0.21, more replacements...quasi gleasonian?

PSRE_10_all_compA_metacom$Boundary
# index=1.51, p<0.0001, clumped losses

#########################
# Compartment B
# rows 404-924, columns 1-7

PSRE_10_all_compB <- PSRE_10_all_Comm[404:924, 1:8]

PSRE_10_all_compB_metacom <- metacommunity(PSRE_10_all_compB, order=T, method="r1", sims=1000)

PSRE_10_all_compB_Comm <- PSRE_10_all_compB_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_10_all_compB_Comm, xlab="Parasite Species", ylab="Host Indiv."))

PSRE_10_all_compB_metacom$Coherence
# p-val = 0.88, but looks like there's two more compartments?


###################################################
######## Metacommunity Analyses: 2011 #############
###################################################

# Remove unneccessary columns, convert to matrix
PSRE_11_mat <- as.matrix(PSRE_11[, 9:27])
head(PSRE_11_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_11_mat)){
  if(sum(PSRE_11_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11_mat <- PSRE_11_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_11_mat)){
  if(sum(PSRE_11_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11_mat <- PSRE_11_mat[, -zeros]

nrow(PSRE_11_mat) # 517 indiv.
ncol(PSRE_11_mat) # 11 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_11_metacom <- metacommunity(PSRE_11_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_11_Comm <- PSRE_11_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_11_Comm, xlab="Parasite Species", ylab="Host Indiv."))
# Looks like no pattern at all...

# View coherence:
PSRE_11_metacom$Coherence
# pval = 0.81, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# 2011 vs. 2012 analysis
# Now I will pull out individuals from sites that were sampled in both 2011 and 2012.

# Create two data frames that have either sites sampled in 2011 or 2012
PSRE_11b <- subset(PSRE, 
                   substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")
PSRE_12 <- subset(PSRE, 
                  substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2012")

# Now figure out which individuals were from sites sampled in 2012, that were also sampled 
# in 2010
remove <- NULL
for(i in 1:nrow(PSRE_12)){
  eval <- NULL
  eval <- substr(PSRE_11b$site.yr, start=1, stop = nchar(PSRE_11b$site.yr)-3) %in%
    substr(PSRE_12$site.yr[i], start=1, stop = nchar(PSRE_12$site.yr[i])-3)
  eval <- as.numeric(eval)
  if(sum(eval)==0){
    remove <- c(remove, i)
  }
}
length(remove)
PSRE_12 <- PSRE_12[-remove, ]

# Now see which 2010 match up to the 2012:

remove <- NULL
for(i in 1:nrow(PSRE_11b)){
  eval <- NULL
  eval <- substr(PSRE_12$site.yr, start=1, stop = nchar(PSRE_12$site.yr)-3) %in%
    substr(PSRE_11b$site.yr[i], start=1, stop = nchar(PSRE_11b$site.yr[i])-3)
  eval <- as.numeric(eval)
  if(sum(eval)==0){
    remove <- c(remove, i)
  }
}
length(remove)
PSRE_11b <- PSRE_11b[-remove, ]

nrow(PSRE_11b) # 289 indiv
nrow(PSRE_12) # 302 indiv

# Make sure the same site.yr's are present in each data frame
unique(PSRE_11b$site.yr)
unique(PSRE_12$site.yr) # 31 sites each

###################################################
######## Metacommunity Analyses: 2011b #############
###################################################

# Remove unneccessary columns, convert to matrix
PSRE_11b_mat <- as.matrix(PSRE_11b[, 9:27])
head(PSRE_11b_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_11b_mat)){
  if(sum(PSRE_11b_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11b_mat <- PSRE_11b_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_11b_mat)){
  if(sum(PSRE_11b_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11b_mat <- PSRE_11b_mat[, -zeros]

nrow(PSRE_11b_mat) # 279 indiv.
ncol(PSRE_11b_mat) # 10 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_11b_metacom <- metacommunity(PSRE_11b_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_11b_Comm <- PSRE_11b_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_11b_Comm, xlab="Parasite Species", ylab="Host Indiv."))
#Actually looks nested

# View coherence:
PSRE_11b_metacom$Coherence
# pval = 0.039

PSRE_11b_metacom$Turnover
# pval = 0.78, fewer replacements than expected (quasi-nested)

PSRE_11b_metacom$Boundary
# index=1.99, p<0.0001

###########################
# Because 2011 and 2011b showed different patterns, I'm going to use the entire 2011 dataset

PSRE_11_all <- subset(PSRE, 
                      substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")

# Remove unneccessary columns, convert to matrix
PSRE_11_all_mat <- as.matrix(PSRE_11_all[, 9:27])
head(PSRE_11_all_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_11_all_mat)){
  if(sum(PSRE_11_all_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11_all_mat <- PSRE_11_all_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_11_all_mat)){
  if(sum(PSRE_11_all_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11_all_mat <- PSRE_11_all_mat[, -zeros]

nrow(PSRE_11_all_mat) # 620 indiv.
ncol(PSRE_11_all_mat) # 11 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_11_all_metacom <- metacommunity(PSRE_11_all_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_11_all_Comm <- PSRE_11_all_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_11_all_Comm, xlab="Parasite Species", ylab="Host Indiv."))
# Two clear compartments.

# View coherence:
PSRE_11_all_metacom$Coherence
# p-val=0.35, no coherence, but nestd subcompartments?

###################################################
######## Metacommunity Analyses: 2012 #############
###################################################

# Remove unneccessary columns, convert to matrix
PSRE_12_mat <- as.matrix(PSRE_12[, 9:27])
head(PSRE_12_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_12_mat)){
  if(sum(PSRE_12_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_12_mat <- PSRE_12_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_12_mat)){
  if(sum(PSRE_12_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_12_mat <- PSRE_12_mat[, -zeros]

nrow(PSRE_12_mat) # 231 indiv.
ncol(PSRE_12_mat) # 7 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_12_metacom <- metacommunity(PSRE_12_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_12_Comm <- PSRE_12_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_12_Comm, xlab="Parasite Species", ylab="Host Indiv."))
#looks sort of nested

# View coherence:
PSRE_12_metacom$Coherence
# pval = 0.069, no coherence

PSRE_12_metacom$Turnover
# pval = 0.59, fewer replacements than expected (quasi-nested)

PSRE_12_metacom$Boundary
# index=2.38, p<0.0001