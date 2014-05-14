# California PSRE parasite metacommunities
# Using occupancy modeling in conjunction with the EMS framework

# Script Author: JR Mihaljevic
# Project co-authors: Bethany J. Hoye, Pieter T.J. Johnson

# Ordinating the data matrices (Y) and the estimated occurrence matrices (Z):

##############################################################################################################
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################

#########################################
#####      Y MATRIX STRUCTURE       #####
#########################################
# First use the metacom() package to ordinate and figure out structure of the 
# input matrices:

# 2009:

mat_2009 <- aperm(Y.obs[, 1:83], c(2, 1))
mat_2009 <- ifelse(mat_2009 > 0, 1, 0)

# Check column sums
# Remove species that occurred less than 5 times:
if(any(colSums(mat_2009) < 10) == TRUE){
  mat_2009 <- mat_2009[, -which(colSums(mat_2009) < 10)]
}

# Check row sums
if(any(rowSums(mat_2009) == 0) == TRUE){
  mat_2009 <- mat_2009[-which(rowSums(mat_2009) == 0), ]
}

meta_2009 <- Metacommunity(mat_2009, method="r1", sims=1000, allow.empty=TRUE)
meta_2009
quartz(height=6, width=3)
print(Matrix_Plot(meta_2009$Comm, xlab="Parasite Species", ylab="Sites"))

# Summary:
# Significant coherence but looks like two distinct compartments.
# Compartments seem to be driving by presence of MANO. 
##############################################################################################################
##############################################################################################################

# Separate compartment A and B:
compA_2009 <- meta_2009$Comm[1:18, 2:5] #18 sites (no alaria)
compB_2009 <- meta_2009$Comm[19:nrow(meta_2009$Comm), 1:4] #(excluded 2 sites w/ Mano)

meta_compA_2009 <- Metacommunity(compA_2009, method="r1", sims=1000, allow.empty=TRUE)
meta_compA_2009
quartz(height=6, width=3)
print(Matrix_Plot(meta_compA_2009$Comm, xlab="Parasite Species", ylab="Sites"))

# CompA: Significant coherence, non-sig turnover (quasi-nested)

meta_compB_2009 <- Metacommunity(compB_2009, method="r1", sims=1000, allow.empty=TRUE)
meta_compB_2009
quartz(height=6, width=3)
print(Matrix_Plot(meta_compB_2009$Comm, xlab="Parasite Species", ylab="Sites"))

# CompB: Significantly nested 

# Summary:
# In 2009, two compartments present, mainly distinguished by the presence of
# Mano/Alaria. In the occupancy model, Mano had a positive effect of snail richness.
# Indeed, compartment A has a signficantly higher average snail richness compared
# to compartment B, which excludes Mano. Furthermore, in the model, Alaria showed
# a positive response to canopy cover, while Mano showed a weakly negative resonse
# to canopy cover. Indeed compartment B (with Alaria) has higher average canopy cover.

# Need to find the best way to show this graphically. 
#------------------------------------------------------------------------------------------------------------#

# 2010:

mat_2010 <- aperm(Y.obs[, 84:188], c(2, 1))
mat_2010 <- ifelse(mat_2010 > 0, 1, 0)

# Check column sums
# Remove species that occurred less than 5 times:
if(any(colSums(mat_2010) < 10) == TRUE){
  mat_2010 <- mat_2010[, -which(colSums(mat_2010) < 10)]
}

# Check row sums
if(any(rowSums(mat_2010) == 0) == TRUE){
  mat_2010 <- mat_2010[-which(rowSums(mat_2010) == 0), ]
}

meta_2010 <- Metacommunity(mat_2010, method="r1", sims=1000, allow.empty=TRUE)
meta_2010
quartz(height=6, width=3)
print(Matrix_Plot(meta_2010$Comm, xlab="Parasite Species", ylab="Sites"))

# Summary:
# Significant coherence but looks like two distinct compartments.
# Compartments seem to be driving by presence of MANO. 
##############################################################################################################
##############################################################################################################

# Separate compartment A and B:
compA_2010 <- meta_2010$Comm[which(meta_2010$Comm[,5]==1), c(1,3:5)] #19 sites
compB_2010 <- meta_2010$Comm[which(meta_2010$Comm[,5]==0), 1:4] #75 sites

meta_compA_2010 <- Metacommunity(compA_2010, method="r1", sims=1000, allow.empty=TRUE)
meta_compA_2010
quartz(height=6, width=3)
print(Matrix_Plot(meta_compA_2010$Comm, xlab="Parasite Species", ylab="Sites"))

# CompA: Too sparse to work...

meta_compB_2010 <- Metacommunity(compB_2010, method="r1", sims=1000, allow.empty=TRUE)
meta_compB_2010
quartz(height=6, width=3)
print(Matrix_Plot(meta_compB_2010$Comm, xlab="Parasite Species", ylab="Sites"))

# CompB: Significantly nested 

# Summary:
# In 2010, two compartments present, mainly distinguished by the presence of
# Mano/abs of Glob. In the occupancy model, Mano had a positive effect of snail richness.
# Indeed, compartment A has a signficantly higher average snail richness compared
# to compartment B, which excludes Mano. Furthermore, in the model, Glob had a negative
# effect of Total_N, while Mano had a weak positive effect. However, no sig. differnce in
# Total_N in sites w/ and without Glob in 2010. 

# Need to find the best way to show this graphically. 
#------------------------------------------------------------------------------------------------------------#

# 2011:

mat_2011 <- aperm(Y.obs[, 189:255], c(2, 1))
mat_2011 <- ifelse(mat_2011 > 0, 1, 0)

# Check column sums
# Remove species that occurred less than 5 times:
if(any(colSums(mat_2011) < 10) == TRUE){
  mat_2011 <- mat_2011[, -which(colSums(mat_2011) < 10)]
}

# Check row sums
if(any(rowSums(mat_2011) == 0) == TRUE){
  mat_2011 <- mat_2011[-which(rowSums(mat_2011) == 0), ]
}

meta_2011 <- Metacommunity(mat_2011, method="r1", sims=1000, allow.empty=TRUE)
meta_2011
quartz(height=6, width=3)
print(Matrix_Plot(meta_2011$Comm, xlab="Parasite Species", ylab="Sites"))

# Summary:
# No coherence. Somewhat separated based on pres/abs of Mano and Fib, but neither
# of these species share a common covariate response. Mano has strong pos effect of elev
# and strong neg effect of hydroperiod, but Fib has strong neg effect of pH.

# Lack of coherence could be result of no common covariate effecting multiple
# species. 
#------------------------------------------------------------------------------------------------------------#

# 2012:

mat_2012 <- aperm(Y.obs[, 256:289], c(2, 1))
mat_2012 <- ifelse(mat_2012 > 0, 1, 0)

# Check column sums
# Remove species that occurred less than 5 times:
if(any(colSums(mat_2012) < 5) == TRUE){
  mat_2012 <- mat_2012[, -which(colSums(mat_2012) < 5)]
}

# Check row sums
if(any(rowSums(mat_2012) == 0) == TRUE){
  mat_2012 <- mat_2012[-which(rowSums(mat_2012) == 0), ]
}

meta_2012 <- Metacommunity(mat_2012, method="r1", sims=1000, allow.empty=TRUE)
meta_2012
quartz(height=6, width=3)
print(Matrix_Plot(meta_2012$Comm, xlab="Parasite Species", ylab="Sites"))

# Summary:
# No coherence. Might have two compartments based on Glob pres?
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################
##############################################################################################################
# Separate compartment A and B:
compA_2012 <- meta_2012$Comm[which(meta_2012$Comm[,1]==1), ] #18 sites
compB_2012 <- meta_2012$Comm[which(meta_2012$Comm[,1]==0), 2:5] #13 sites

meta_compA_2012 <- Metacommunity(compA_2012, method="r1", sims=1000, allow.empty=TRUE)
meta_compA_2012
quartz(height=6, width=3)
print(Matrix_Plot(meta_compA_2012$Comm, xlab="Parasite Species", ylab="Sites"))

# CompA: Quasi-nested

meta_compB_2012 <- Metacommunity(compB_2012, method="r1", sims=1000, allow.empty=TRUE)
meta_compB_2012
quartz(height=6, width=3)
print(Matrix_Plot(OrderMatrix(compB_2012), xlab="Parasite Species", ylab="Sites"))

# CompB: Too sparse to work... but looks very nested. 

# Summary: Glob negatively affected by pH, but so are most other parasites
# Not confident in that answer. Probably sample size too small. 

##############################################################################################################
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################

#########################################
#####      Z MATRIX STRUCTURE       #####
#########################################

# Now I'll represent Z-matrix structure (structure of Bayesian estimates for true occurrence)

# z matrices stored in data.frame "post.z"
Zposterior <- mat.or.vec(nr=nrow(Y.obs), nc=ncol(Y.obs))

for(i in 1:1000){ # 1000 draws from the posterior
  
  # Randomly choose from which chain the sample will originate:
  chain <- NULL
  chain <- sample(c(1:3), 1)
  
  # Create a subset vector of the values for all z[n, k]:
  subset <- NULL
  subset <- subset(post.z, Chain==chain & Iteration==paste(i))$value
  # Store this vector as a matrix
  Z_iter <- NULL
  Z_iter <- matrix(subset, nrow=nrow(Y.obs), ncol=ncol(Y.obs), byrow=F)
  
  Zposterior <- Zposterior + Z_iter
}  

Zposterior <- Zposterior/695 # This will make an average prob. occurrence
rownames(Zposterior) <- rownames(Y.obs)
colnames(Zposterior) <- colnames(Y.obs)
#------------------------------------------------------------------------------------------------------------#
# Look at 2009:
Z_2009 <- aperm(Zposterior, c(2,1))[1:83, ]

# Remove Thic: too problematic
Z_2009 <- Z_2009[, -8]
# Check column sums
# Remove species that occurred less than 10 times:
if(any(colSums(Z_2009) < 10) == TRUE){
  Z_2009 <- Z_2009[, -which(colSums(Z_2009) < 10)]
}

# Check row sums
if(any(rowSums(Z_2009) < 1) == TRUE){
  Z_2009 <- Z_2009[-which(rowSums(Z_2009) < 1), ]
}

# Conduct a CCA and order the sites based on the first axes scores:
library(vegan)

Z_2009.CCA <- decorana(Z_2009, ira=0) # CHECK CONVERGENCE, ETC.
Z_2009.CCA.sites <- Z_2009.CCA$rproj[, 1]
Z_2009.CCA.spp <- Z_2009.CCA$cproj[, 1]

Z_2009_ord <- Z_2009[order(Z_2009.CCA.sites, decreasing = FALSE), 
                     order(Z_2009.CCA.spp, decreasing = FALSE)]

quartz(height=6, width=3)
print(Matrix_HeatMap_NMDS(Z_2009_ord, xlab="Parasite Species", ylab="Sites"))
#------------------------------------------------------------------------------------------------------------#
# Look at 2010:
Z_2010 <- aperm(Zposterior, c(2,1))[84:188, ]

# Check effect of removing Fib:
Z_2010 <- Z_2010[, -3]

# Check column sums
# Remove species that occurred less than 10 times:
if(any(colSums(Z_2010) < 10) == TRUE){
  Z_2010 <- Z_2010[, -which(colSums(Z_2010) < 10)]
}

# Check row sums
if(any(rowSums(Z_2010) < 1) == TRUE){
  Z_2010 <- Z_2010[-which(rowSums(Z_2010) < 1), ]
}

# Conduct a CCA and order the sites based on the first axes scores:
library(vegan)

Z_2010.CCA <- decorana(Z_2010, ira=0) # CHECK CONVERGENCE, ETC.
Z_2010.CCA.sites <- Z_2010.CCA$rproj[, 1]
Z_2010.CCA.spp <- Z_2010.CCA$cproj[, 1]

Z_2010_ord <- Z_2010[order(Z_2010.CCA.sites, decreasing = FALSE), 
                     order(Z_2010.CCA.spp, decreasing = FALSE)]

quartz(height=6, width=3)
print(Matrix_HeatMap_NMDS(Z_2010_ord, xlab="Parasite Species", ylab="Sites"))

#------------------------------------------------------------------------------------------------------------#
# Look at 2011:
Z_2011 <- aperm(Zposterior, c(2,1))[189:255, ]

# Check column sums
# Remove species that occurred less than 10 times:
if(any(colSums(Z_2011) < 10) == TRUE){
  Z_2011 <- Z_2011[, -which(colSums(Z_2011) < 10)]
}

# Check row sums
if(any(rowSums(Z_2011) < 1) == TRUE){
  Z_2011 <- Z_2011[-which(rowSums(Z_2011) < 1), ]
}

# Conduct a CCA and order the sites based on the first axes scores:
library(vegan)

Z_2011.CCA <- decorana(Z_2011, ira=0) # CHECK CONVERGENCE, ETC.
Z_2011.CCA.sites <- Z_2011.CCA$rproj[, 1]
Z_2011.CCA.spp <- Z_2011.CCA$cproj[, 1]

Z_2011_ord <- Z_2011[order(Z_2011.CCA.sites, decreasing = FALSE), 
                     order(Z_2011.CCA.spp, decreasing = FALSE)]

quartz(height=6, width=3)
print(Matrix_HeatMap_NMDS(Z_2011_ord, xlab="Parasite Species", ylab="Sites"))

#------------------------------------------------------------------------------------------------------------#
# Look at 2012:
Z_2012 <- aperm(Zposterior, c(2,1))[256:289, ]

#Thic is problematic:
Z_2012 <- Z_2012[, -8]
# Check column sums
# Remove species that occurred less than 10 times:
if(any(colSums(Z_2012) < 5) == TRUE){
  Z_2012 <- Z_2012[, -which(colSums(Z_2012) < 5)]
}

# Check row sums
if(any(rowSums(Z_2012) < 1) == TRUE){
  Z_2012 <- Z_2012[-which(rowSums(Z_2012) < 1), ]
}

# Conduct a CCA and order the sites based on the first axes scores:
library(vegan)

Z_2012.CCA <- decorana(Z_2012, ira=0) # CHECK CONVERGENCE, ETC.
Z_2012.CCA.sites <- Z_2012.CCA$rproj[, 1]
Z_2012.CCA.spp <- Z_2012.CCA$cproj[, 1]

Z_2012_ord <- Z_2012[order(Z_2012.CCA.sites, decreasing = FALSE), 
                     order(Z_2012.CCA.spp, decreasing = FALSE)]

quartz(height=6, width=3)
print(Matrix_HeatMap_NMDS(Z_2012_ord, xlab="Parasite Species", ylab="Sites"))