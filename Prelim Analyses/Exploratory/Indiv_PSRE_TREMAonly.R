# California parasite metacommunities
# Analysis of metacommunity structure at the individual host level 
# Including only trematode species, and only PSRE
# Includes only sites with >10 PSRE dissected. 
# August 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Remove unnecessary data, and empty columns
indiv_trema <- indiv_PSRE[, -c(2:8)]
summary(indiv_trema)
zeros <- NULL
for (j in 2:ncol(indiv_trema)){
  if(sum(indiv_trema[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
indiv_trema <- indiv_trema[, -zeros]

# Get rid of non-trematode species
indiv_trema <- indiv_trema[, -c(12:14, 17)]
head(indiv_trema)

#Remove first row
indiv_trema_mat <- indiv_trema[, -1]
indiv_trema_mat <- as.matrix(indiv_trema_mat)

# Remove rows with sum = 0
zeros <- NULL
for (j in 1:nrow(indiv_trema_mat)){
  if(sum(indiv_trema_mat[j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros

indiv_trema_mat <- indiv_trema_mat[-zeros ,]

colnames(indiv_trema_mat) <- NULL
head(indiv_trema_mat)

ncol(indiv_trema_mat) #12 parasite species
nrow(indiv_trema_mat) #2130 hosts

##### METACOMMUNITY ANALYSES ####
# for all PSRE, only trematodes

library(metacom)
ordinated_indiv_trema <- OrderMatrix(indiv_trema_mat)
image(apply(ordinated_indiv_trema,2,rev), col=c(0,'black'), axes=FALSE, xlab='species', ylab='sites', cex.lab=1.4)
box(lwd=2,col=1)

#Determine coherence. Must be significant to continue with analyses. Otherwise
#the metacommunity is assembled randomly
indiv_trema_coherence <- Coherence(indiv_trema_mat, 
                                  sims=1000, method="r1", allow.empty=F)
indiv_trema_coherence
# p-val = 0.96...definitely no coherence there

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for PSRE hosts:
# All PSRE dissected in 2009, only trematodes, sites with >=10 dissected PSRE

# Need to convert 'site.yr' as character string, then partition dataset
indiv_trema$site.yr <- as.character(indiv_trema$site.yr)
indiv_trema_2009 <- subset(indiv_trema, 
                          substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2009")
head(indiv_trema_2009)

# Need to remove headings and first column
indiv_trema_2009 <- indiv_trema_2009[, -1]
colnames(indiv_trema_2009) <- NULL
head(indiv_trema_2009)

# Turn into a matrix
indiv_trema_2009 <- as.matrix(indiv_trema_2009)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_trema_2009)){
  if(sum(indiv_trema_2009[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_trema_2009 <- indiv_trema_2009[, -zeros]

# Check row sums as well
zeros <- NULL
for (j in 1:nrow(indiv_trema_2009)){
  if(sum(indiv_trema_2009[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
indiv_trema_2009 <- indiv_trema_2009[-zeros ,]

ncol(indiv_trema_2009) #9 parasite species
nrow(indiv_trema_2009) #742 hosts


#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
indiv_trema_2009_metacom <- metacommunity(indiv_trema_2009, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_indiv_trema_2009 <- indiv_trema_2009_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(ordinated_indiv_trema_2009))
# Looks sorta nested

indiv_trema_2009_metacom$Coherence
# pval = 0.46...no coherence   

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for PSRE hosts:
# All PSRE dissected in 2010, only trematodes, sites with >=10 dissected PSRE

# Need to partition dataset

indiv_trema_2010 <- subset(indiv_trema, 
                           substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2010")
head(indiv_trema_2010)

# Need to remove headings and first column
indiv_trema_2010 <- indiv_trema_2010[, -1]
head(indiv_trema_2010)

# Turn into a matrix
indiv_trema_2010 <- as.matrix(indiv_trema_2010)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_trema_2010)){
  if(sum(indiv_trema_2010[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_trema_2010 <- indiv_trema_2010[, -zeros]

# Check row sums as well
zeros <- NULL
for (j in 1:nrow(indiv_trema_2010)){
  if(sum(indiv_trema_2010[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
indiv_trema_2010 <- indiv_trema_2010[-zeros ,]

head(indiv_trema_2010)
ncol(indiv_trema_2010) #10 parasite species
nrow(indiv_trema_2010) #654 hosts


#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
indiv_trema_2010_metacom <- metacommunity(indiv_trema_2010, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_indiv_trema_2010 <- indiv_trema_2010_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(ordinated_indiv_trema_2010))
# Looks nested

indiv_trema_2010_metacom$Coherence
# pval = 0.061...no coherence, but borderline, definitely looks nested.

indiv_trema_2010_metacom$Turnover
#p-val = 0.18...no turnover (quasi-nested)

indiv_trema_2010_metacom$Boundary
# index=3.67, p<0.0001, clumped losses

# Save the site scores of the ordinated matrix 
scores_PSRE_2010_trema <- OrderMatrix(ordinated_indiv_trema_2010, output.scores=T)
sitescores_PSRE_2010_trema <- scores_PSRE_2010_trema$sitescores
sitescores_PSRE_2010_trema <- melt(sitescores_PSRE_2010_trema)
HostID <- rownames(sitescores_PSRE_2010_trema)
sitescores_PSRE_2010_trema <- data.frame(HostID, sitescores_PSRE_2010_trema)
head(sitescores_PSRE_2010_trema)

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for PSRE hosts:
# All PSRE dissected in 2011, only trematodes, sites with >=10 dissected PSRE

# Need to partition dataset

indiv_trema_2011 <- subset(indiv_trema, 
                           substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")
head(indiv_trema_2011)

# Need to remove headings and first column
indiv_trema_2011 <- indiv_trema_2011[, -1]
colnames(indiv_trema_2011) <- NULL
head(indiv_trema_2011)

# Turn into a matrix
indiv_trema_2011 <- as.matrix(indiv_trema_2011)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_trema_2011)){
  if(sum(indiv_trema_2011[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_trema_2011 <- indiv_trema_2011[, -zeros]

# Check row sums as well
zeros <- NULL
for (j in 1:nrow(indiv_trema_2011)){
  if(sum(indiv_trema_2011[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
indiv_trema_2011 <- indiv_trema_2011[-zeros ,]

head(indiv_trema_2011)
ncol(indiv_trema_2011) #7 parasite species
nrow(indiv_trema_2011) #489 hosts


#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
indiv_trema_2011_metacom <- metacommunity(indiv_trema_2011, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_indiv_trema_2011 <- indiv_trema_2011_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(ordinated_indiv_trema_2011))
# Looks nested, with two compartments?

indiv_trema_2011_metacom$Coherence
# pval = 0.0024, positive coherence

indiv_trema_2011_metacom$Turnover
# pval = 0.016, positive turnover (more than expected)

indiv_trema_2011_metacom$Boundary
# index=2.13, p<0.0001, so clementsian gradient. 

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for PSRE hosts:
# All PSRE dissected in 2012, only trematodes, sites with >=10 dissected PSRE

# Need to partition dataset

indiv_trema_2012 <- subset(indiv_trema, 
                           substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2012")
head(indiv_trema_2012)

# Need to remove headings and first column
indiv_trema_2012 <- indiv_trema_2012[, -1]
head(indiv_trema_2012)

# Turn into a matrix
indiv_trema_2012 <- as.matrix(indiv_trema_2012)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_trema_2012)){
  if(sum(indiv_trema_2012[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_trema_2012 <- indiv_trema_2012[, -zeros]

# Check row sums as well
zeros <- NULL
for (j in 1:nrow(indiv_trema_2012)){
  if(sum(indiv_trema_2012[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
indiv_trema_2012 <- indiv_trema_2012[-zeros ,]

ncol(indiv_trema_2012) #9 parasite species
nrow(indiv_trema_2012) #245 hosts


#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
indiv_trema_2012_metacom <- metacommunity(indiv_trema_2012, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_indiv_trema_2012 <- indiv_trema_2012_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(ordinated_indiv_trema_2012))
# Looks slightly clementsian? 

indiv_trema_2012_metacom$Coherence
# pval = 0.00038, positive coherence

indiv_trema_2012_metacom$Turnover
# pval = 0.504, not significant, but turnover tends towards less than expected

indiv_trema_2012_metacom$Boundary
# index=2.69, p<0.0001, clumped species losses

# So, we have quasi-nested with clumped species losses. 
# However, looks like there are two compartments, both nested?

######################################
######################################
# First compartment, rows 1-112

indiv_trema_compA2012 <- ordinated_indiv_trema_2012[1:112 ,]
# Check for column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_trema_compA2012)){
  if(sum(indiv_trema_compA2012[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_trema_compA2012 <- indiv_trema_compA2012[, -zeros]

nrow(indiv_trema_compA2012) # 112 hosts
ncol(indiv_trema_compA2012) # 4 parasites

### Metacommunity Analyses ###
indiv_trema_compA2012_metacom <- metacommunity(indiv_trema_compA2012, 
                                               order=F, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ord_indiv_trema_compA2012 <- indiv_trema_compA2012_metacom$Comm
quartz(height=3, width=3)
print(Matrix_Plot(ord_indiv_trema_compA2012))
# Looks  

indiv_trema_compA2012_metacom$Coherence
# pval = 0.063, borderline positive coherence

indiv_trema_compA2012_metacom$Turnover
# pval = 0.13, not significant, quasi-nested

indiv_trema_compA2012_metacom$Boundary
# index=1.62, p<0.0001, clumped species losses

######################################
######################################
# First compartment, rows 1-112

indiv_trema_compB2012 <- ordinated_indiv_trema_2012[113:245 ,]
# Check for column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_trema_compB2012)){
  if(sum(indiv_trema_compB2012[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_trema_compB2012 <- indiv_trema_compB2012[, -zeros]

nrow(indiv_trema_compB2012) #133 hosts
ncol(indiv_trema_compB2012) #7 parasites

### Metacommunity Analyses ###
indiv_trema_compB2012_metacom <- metacommunity(indiv_trema_compB2012, 
                                               order=F, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ord_indiv_trema_compB2012 <- indiv_trema_compB2012_metacom$Comm
quartz(height=3, width=3)
print(Matrix_Plot(ord_indiv_trema_compB2012))
# Looks  

indiv_trema_compB2012_metacom$Coherence
# pval = 0.010, positive coherence

indiv_trema_compB2012_metacom$Turnover
# pval = 0.55, not significant, but turnover tends towards more than expected

indiv_trema_compB2012_metacom$Boundary
# index=1.81, p<0.0001, clumped species losses
# quasi-Clementsian structure?

######################################
# Save the site scores of the ordinated matrix 
scores_PSRE_2012_trema <- OrderMatrix(ordinated_indiv_trema_2012, output.scores=T)
sitescores_PSRE_2012_trema <- scores_PSRE_2012_trema$sitescores
sitescores_PSRE_2012_trema <- melt(sitescores_PSRE_2012_trema)
HostID <- rownames(sitescores_PSRE_2012_trema)
sitescores_PSRE_2012_trema <- data.frame(HostID, sitescores_PSRE_2012_trema)
head(sitescores_PSRE_2012_trema)

