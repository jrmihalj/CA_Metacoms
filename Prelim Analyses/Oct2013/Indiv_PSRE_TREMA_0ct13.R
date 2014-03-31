# California parasite metacommunities
# Analysis of metacommunity structure at the inidividual (PSRE)-level 
# Including only TREMATODE parasite species acquired within the pond
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
# ***Nyct - Protozoa ***Removed***
# ***Opal - Protist ***Removed***
# ***Pinw - Nematode ***Removed***
# Rib - Trematode
# Thic - Trematode


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

###################################################
######## Metacommunity Analyses: 2009 #############
###################################################

head(PSRE_9)

# Remove the three parasites above:
PSRE_9_trema <- PSRE_9[, -c(20,21,23)]
head(PSRE_9_trema)

# Remove unneccessary columns, convert to matrix
PSRE_9_trema_mat <- as.matrix(PSRE_9_trema[, 9:24])
head(PSRE_9_trema_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_9_trema_mat)){
  if(sum(PSRE_9_trema_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_9_trema_mat <- PSRE_9_trema_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_9_trema_mat)){
  if(sum(PSRE_9_trema_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_9_trema_mat <- PSRE_9_trema_mat[, -zeros]

nrow(PSRE_9_trema_mat) # 346 indiv.
ncol(PSRE_9_trema_mat) # 7 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_9_trema_metacom <- metacommunity(PSRE_9_trema_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_9_trema_Comm <- PSRE_9_trema_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_9_trema_Comm, xlab="Parasite Species", ylab="Host Indiv."))
# Maybe nested, but hard to say

# View coherence:
PSRE_9_trema_metacom$Coherence
# pval = 0.49, no coherence

###################################################
######## Metacommunity Analyses: 2010a ############
###################################################

# Remove the three parasites above:
PSRE_10a_trema <- PSRE_10[, -c(20,21,23)]
head(PSRE_10a_trema)

# Remove unneccessary columns, convert to matrix
PSRE_10a_trema_mat <- as.matrix(PSRE_10a_trema[, 9:24])
head(PSRE_10a_trema_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_10a_trema_mat)){
  if(sum(PSRE_10a_trema_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_10a_trema_mat <- PSRE_10a_trema_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_10a_trema_mat)){
  if(sum(PSRE_10a_trema_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_10a_trema_mat <- PSRE_10a_trema_mat[, -zeros]

nrow(PSRE_10a_trema_mat) # 276 indiv.
ncol(PSRE_10a_trema_mat) # 6 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_10a_trema_metacom <- metacommunity(PSRE_10a_trema_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_10a_trema_Comm <- PSRE_10a_trema_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_10a_trema_Comm, xlab="Parasite Species", ylab="Host Indiv."))

# View coherence:
PSRE_10a_trema_metacom$Coherence
# pval = 0.95, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

###################################################
######## Metacommunity Analyses: 2010b #############
###################################################

# Remove the three parasites above:
PSRE_10b_trema <- PSRE_10b[, -c(20,21,23)]
head(PSRE_10b_trema)

# Remove unneccessary columns, convert to matrix
PSRE_10b_trema_mat <- as.matrix(PSRE_10b_trema[, 9:24])
head(PSRE_10b_trema_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_10b_trema_mat)){
  if(sum(PSRE_10b_trema_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_10b_trema_mat <- PSRE_10b_trema_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_10b_trema_mat)){
  if(sum(PSRE_10b_trema_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_10b_trema_mat <- PSRE_10b_trema_mat[, -zeros]

nrow(PSRE_10b_trema_mat) # 421 indiv.
ncol(PSRE_10b_trema_mat) # 7 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_10b_trema_metacom <- metacommunity(PSRE_10b_trema_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_10b_trema_Comm <- PSRE_10b_trema_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_10b_trema_Comm, xlab="Parasite Species", ylab="Host Indiv."))

# View coherence:
PSRE_10b_trema_metacom$Coherence
# pval = 0.92, no coherence

###################################################
######## Metacommunity Analyses: 2011 #############
###################################################

# Remove the three parasites above:
PSRE_11_trema <- PSRE_11[, -c(20,21,23)]
head(PSRE_11_trema)

# Remove unneccessary columns, convert to matrix
PSRE_11_trema_mat <- as.matrix(PSRE_11_trema[, 9:24])
head(PSRE_11_trema_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_11_trema_mat)){
  if(sum(PSRE_11_trema_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11_trema_mat <- PSRE_11_trema_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_11_trema_mat)){
  if(sum(PSRE_11_trema_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11_trema_mat <- PSRE_11_trema_mat[, -zeros]

nrow(PSRE_11_trema_mat) # 470 indiv.
ncol(PSRE_11_trema_mat) # 8 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_11_trema_metacom <- metacommunity(PSRE_11_trema_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_11_trema_Comm <- PSRE_11_trema_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_11_trema_Comm, xlab="Parasite Species", ylab="Host Indiv."))

# View coherence:
PSRE_11_trema_metacom$Coherence
# pval = 0.37, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


###################################################
######## Metacommunity Analyses: 2011b #############
###################################################

# Remove the three parasites above:
PSRE_11b_trema <- PSRE_11b[, -c(20,21,23)]
head(PSRE_11b_trema)

# Remove unneccessary columns, convert to matrix
PSRE_11b_trema_mat <- as.matrix(PSRE_11b_trema[, 9:24])
head(PSRE_11b_trema_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_11b_trema_mat)){
  if(sum(PSRE_11b_trema_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11b_trema_mat <- PSRE_11b_trema_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_11b_trema_mat)){
  if(sum(PSRE_11b_trema_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11b_trema_mat <- PSRE_11b_trema_mat[, -zeros]

nrow(PSRE_11b_trema_mat) # 248 indiv.
ncol(PSRE_11b_trema_mat) # 7 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_11b_trema_metacom <- metacommunity(PSRE_11b_trema_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_11b_trema_Comm <- PSRE_11b_trema_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_11b_trema_Comm, xlab="Parasite Species", ylab="Host Indiv."))
# Looks nested

# View coherence:
PSRE_11b_trema_metacom$Coherence
# pval = 0.023, positive coherence

PSRE_11b_trema_metacom$Turnover
# p-val: 0.26 (many fewer replacements than mean) quasi-nested

PSRE_11b_trema_metacom$Boundary
# index: 2.15 p<0.0001 clumped losses

###################################################
######## Metacommunity Analyses: 2011_all #############
###################################################

# Different structures appeared in 2011a and 11b, so do an analysis with all 2011 samples

PSRE_11_all <- subset(PSRE, 
                      substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")

# Remove the three parasites above:
PSRE_11_all_trema <- PSRE_11_all[, -c(20,21,23)]
head(PSRE_11_all_trema)

# Remove unneccessary columns, convert to matrix
PSRE_11_all_trema_mat <- as.matrix(PSRE_11_all_trema[, 9:24])
head(PSRE_11_all_trema_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_11_all_trema_mat)){
  if(sum(PSRE_11_all_trema_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11_all_trema_mat <- PSRE_11_all_trema_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_11_all_trema_mat)){
  if(sum(PSRE_11_all_trema_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_11_all_trema_mat <- PSRE_11_all_trema_mat[, -zeros]

nrow(PSRE_11_all_trema_mat) # 569 indiv.
ncol(PSRE_11_all_trema_mat) # 8 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_11_all_trema_metacom <- metacommunity(PSRE_11_all_trema_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_11_all_trema_Comm <- PSRE_11_all_trema_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_11_all_trema_Comm, xlab="Parasite Species", ylab="Host Indiv."))
# Looks weird

# View coherence:
PSRE_11_all_trema_metacom$Coherence
# pval = 0.43, no coherence

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


###################################################
######## Metacommunity Analyses: 2012 #############
###################################################

# Remove the three parasites above:
PSRE_12_trema <- PSRE_12[, -c(20,21,23)]
head(PSRE_12_trema)

# Remove unneccessary columns, convert to matrix
PSRE_12_trema_mat <- as.matrix(PSRE_12_trema[, 9:24])
head(PSRE_12_trema_mat)

# Make sure no rows or columns sum to zero
zeros <- NULL
for (i in 1:nrow(PSRE_12_trema_mat)){
  if(sum(PSRE_12_trema_mat[i, ]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_12_trema_mat <- PSRE_12_trema_mat[-zeros, ]

zeros <- NULL
for (i in 1:ncol(PSRE_12_trema_mat)){
  if(sum(PSRE_12_trema_mat[, i]) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
PSRE_12_trema_mat <- PSRE_12_trema_mat[, -zeros]

nrow(PSRE_12_trema_mat) # 231 indiv.
ncol(PSRE_12_trema_mat) # 7 parasites

library(metacom)

# Run the metacommunity analyses
PSRE_12_trema_metacom <- metacommunity(PSRE_12_trema_mat, method="r1", sims=1000)

# Store and view the ordinated incidence matrix
PSRE_12_trema_Comm <- PSRE_12_trema_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(PSRE_12_trema_Comm, xlab="Parasite Species", ylab="Host Indiv."))
# Looks nested

# View coherence:
PSRE_12_trema_metacom$Coherence
# pval = 0.058, positive coherence

PSRE_12_trema_metacom$Turnover
# p-val: 0.59 (fewer replacements than mean) quasi-nested

PSRE_12_trema_metacom$Boundary
# index: 2.38 p<0.0001 clumped losses
