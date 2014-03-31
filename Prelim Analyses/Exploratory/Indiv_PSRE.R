# California parasite metacommunities
# Analysis of metacommunity structure at the individual host level 
# Including all parasite species, but only PSRE
# Includes only sites with >10 PSRE dissected. 
# August 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# First import a dataset with all aspects of individual host infection/dissection

# File name: 'Indiv_ALL_Hosts.csv'
indiv_allhosts <- read.csv(file.choose(), header=T, na.strings="")
head(indiv_allhosts)

# Now, subset so that we only have PSRE data

indiv_PSRE <- subset(indiv_allhosts, Spp=="PSRE")
summary(indiv_PSRE$Spp)

# Now, figure out if the PSRE individuals come from sites that were used in the 
# site-level analyses (i.e. sites with >= 10 dissected PSRE)
zeros <- NULL
for ( i in 1:nrow(indiv_PSRE)){
  out <- NULL
  out <- as.numeric(all_sites$site.yr %in% indiv_PSRE[i, 1])
  if (sum(out) == 0){
    zeros <- c(zeros, i)
  }
}
zeros
# Get rid of these rows with zeros

indiv_PSRE <- indiv_PSRE[-zeros, ]
nrow(indiv_PSRE)
head(indiv_PSRE)

# Make a matrix with just the parasite incidence in each individual

indiv_PSRE_mat <- indiv_PSRE[, -c(1:8)]
head(indiv_PSRE_mat)
indiv_PSRE_mat <- as.matrix(indiv_PSRE_mat)

# Get rid of rows that sum to zero (i.e. uninfected individuals)

zeros <- NULL
for (j in 1:nrow(indiv_PSRE_mat)){
  if(sum(indiv_PSRE_mat[j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of uninfected rows
indiv_PSRE_mat <- indiv_PSRE_mat[-zeros ,]

# Get rid of columns that sum to zero (i.e. parasites that are not present anywhere)

zeros <- NULL
for (j in 1:ncol(indiv_PSRE_mat)){
  if(sum(indiv_PSRE_mat[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of uninfected rows
indiv_PSRE_mat <- indiv_PSRE_mat[, -zeros]

ncol(indiv_PSRE_mat) # 16 parasite species present in full dataset
nrow(indiv_PSRE_mat) # 2328 hosts (sites) present in full dataset

##### METACOMMUNITY ANALYSES ####
# for all PSRE hosts, all Parasites

library(metacom)
ordinated_allPSRE <- OrderMatrix(indiv_PSRE_mat)
ordinated_allPSRE_plot <- Matrix_Plot(ordinated_allPSRE)
quartz(height=6, width=3)
print(ordinated_allPSRE_plot) 

#Determine coherence. Must be significant to continue with analyses. Otherwise
#the metacommunity is assembled randomly
indiv_all_coher <- Coherence(indiv_PSRE_mat, 
                                  sims=1000, method="r1", allow.empty=F)
indiv_all_coher
# p-val = 0.53...no coherence. 


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for PSRE hosts:
# All PSRE dissected in 2009, all parasites, sites with >=10 dissected PSRE

# Need to convert 'site.yr' as character string, then partition dataset
indiv_PSRE$site.yr <- as.character(indiv_PSRE$site.yr)
indiv_PSRE_2009 <- subset(indiv_PSRE, 
                     substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2009")
head(indiv_PSRE_2009)

# Need to remove headings and first column
indiv_PSRE_2009 <- indiv_PSRE_2009[, -c(1:8)]
colnames(indiv_PSRE_2009) <- NULL
head(indiv_PSRE_2009)

# Turn into a matrix
indiv_PSRE_2009 <- as.matrix(indiv_PSRE_2009)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_PSRE_2009)){
  if(sum(indiv_PSRE_2009[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_PSRE_2009 <- indiv_PSRE_2009[, -zeros]

# Check row sums as well
zeros <- NULL
for (j in 1:nrow(indiv_PSRE_2009)){
  if(sum(indiv_PSRE_2009[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
indiv_PSRE_2009 <- indiv_PSRE_2009[-zeros ,]

ncol(indiv_PSRE_2009) #12 parasite species
nrow(indiv_PSRE_2009) #754 hosts


#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
indiv_PSRE_2009_metacom <- metacommunity(indiv_PSRE_2009, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_indiv_PSRE_2009 <- indiv_PSRE_2009_metacom$Comm
image(apply(ordinated_indiv_PSRE_2009,2,rev), col=c(0,'black'), axes=FALSE, xlab='species', ylab='sites', cex.lab=1.4)
box(lwd=2,col=1)
# Maybe a pattern? 

indiv_PSRE_2009_metacom$Coherence
# pval = 0.51...no coherence. Shoot. 

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for PSRE hosts:
# All PSRE dissected in 2010, all parasites, sites with >=10 dissected PSRE

# Parition dataset:
indiv_PSRE_2010 <- subset(indiv_PSRE, 
                          substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2010")
head(indiv_PSRE_2010)

# Need to remove headings and first column
indiv_PSRE_2010 <- indiv_PSRE_2010[, -c(1:8)]
head(indiv_PSRE_2010)

# Turn into a matrix
indiv_PSRE_2010 <- as.matrix(indiv_PSRE_2010)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_PSRE_2010)){
  if(sum(indiv_PSRE_2010[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_PSRE_2010 <- indiv_PSRE_2010[, -zeros]
head(indiv_PSRE_2010)

# Check row sums as well
zeros <- NULL
for (j in 1:nrow(indiv_PSRE_2010)){
  if(sum(indiv_PSRE_2010[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
indiv_PSRE_2010 <- indiv_PSRE_2010[-zeros ,]

ncol(indiv_PSRE_2010) #13 parasite species
nrow(indiv_PSRE_2010) #795 individuals


#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
indiv_PSRE_2010_metacom <- metacommunity(indiv_PSRE_2010, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_indiv_PSRE_2010 <- indiv_PSRE_2010_metacom$Comm
ord_indiv_2010_plot <- Matrix_Plot(ordinated_indiv_PSRE_2010)
quartz(height=6, width=3)
print(ord_indiv_2010_plot)

indiv_PSRE_2010_metacom$Coherence
# p-val=0.078, very close to positive coherence...need to split into two compartments

indiv_PSRE_2010_metacom$Turnover
# p-val=0.0047, high turnover (anti-nested)..but looks like nested compartments

indiv_PSRE_2010_metacom$Boundary
# index=2.811, p<0.0001, highly clumped species losses

############
# Overall, shows clemensian structure: positive coherence, high 
# (positive) turnover, and clumped boundary losses.
# However, two compartments clearly visible and need to be analyzed
# separately: appear to be nested (positive coherence, negative turnover
# and evenly or clumped species losses). 

###############################################
# First compartment: rows 1-444
comp1_2010_metacom <- metacommunity(ordinated_indiv_PSRE_2010[1:444 ,], method="r1", sims=1000,
                                    order=F) #DON'T ORDER, BECAUSE IT ALREADY IS

# Store and view the ordinated occupancy matrix
ordinated_comp1_2010_Comm <- comp1_2010_metacom$Comm
comp1_plot <- Matrix_Plot(ordinated_comp1_2010_Comm)
quartz(height=3, width=3)
print(comp1_plot)

# Coherence
comp1_2010_metacom$Coherence
# p-val=0.061, very close to positive coherence

comp1_2010_metacom$Turnover
# p-val=0.10, no turnover, but more negative (way lower than the mean, but high variance)

comp1_2010_metacom$Boundary
# index=4.29, very clumped
# So, quasi-nested with clumped losses, but coherence not technically significant

###############################################
# Second compartment: rows 445-795, col 1-6
comp2_2010_metacom <- metacommunity(ordinated_indiv_PSRE_2010[445:795 , 1:6], method="r1", 
                                    sims=1000,
                                    order=F) #DON'T ORDER, BECAUSE IT ALREADY IS

# Store and view the ordinated occupancy matrix
ordinated_comp2_2010_Comm <- comp2_2010_metacom$Comm
comp2_plot <- Matrix_Plot(ordinated_comp2_2010_Comm)
quartz(height=3, width=3)
print(comp2_plot)

# Coherence
comp2_2010_metacom$Coherence
# p-val=0.014, positive coherence

comp2_2010_metacom$Turnover
# p-val=0.14, no turnover, but more negative (way lower than the mean, but high variance)

comp2_2010_metacom$Boundary
# index=1.18, p<0.0001, very clumped
# So, quasi-nested with clumped losses


# Save the site scores of the ordinated matrix 
scores_PSRE_2010 <- OrderMatrix(ordinated_indiv_PSRE_2010, output.scores=T)
sitescores_PSRE_2010 <- scores_PSRE_2010$sitescores
sitescores_PSRE_2010 <- melt(sitescores_PSRE_2010)
HostID <- rownames(sitescores_PSRE_2010)
sitescores_PSRE_2010 <- data.frame(HostID, sitescores_PSRE_2010)
head(sitescores_PSRE_2010)

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for PSRE hosts:
# All PSRE dissected in 2011, all parasites, sites with >=10 dissected PSRE

# Parition dataset:
indiv_PSRE_2011 <- subset(indiv_PSRE, 
                          substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2011")
head(indiv_PSRE_2011)

# Need to remove headings and first column
indiv_PSRE_2011 <- indiv_PSRE_2011[, -c(1:8)]
colnames(indiv_PSRE_2011) <- NULL
head(indiv_PSRE_2011)

# Turn into a matrix
indiv_PSRE_2011 <- as.matrix(indiv_PSRE_2011)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_PSRE_2011)){
  if(sum(indiv_PSRE_2011[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_PSRE_2011 <- indiv_PSRE_2011[, -zeros]

# Check row sums as well
zeros <- NULL
for (j in 1:nrow(indiv_PSRE_2011)){
  if(sum(indiv_PSRE_2011[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
indiv_PSRE_2011 <- indiv_PSRE_2011[-zeros ,]

ncol(indiv_PSRE_2011) #11 parasite species
nrow(indiv_PSRE_2011) #534 hosts


#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
indiv_PSRE_2011_metacom <- metacommunity(indiv_PSRE_2011, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_indiv_PSRE_2011 <- indiv_PSRE_2011_metacom$Comm
image(apply(ordinated_indiv_PSRE_2011,2,rev), col=c(0,'black'), axes=FALSE, xlab='species', ylab='sites', cex.lab=1.4)
box(lwd=2,col=1)
# Matrix looks nested

indiv_PSRE_2011_metacom$Coherence
#p-val=0.67, even though looks nested. Perhaps would work if only trematodes?

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now check the coherence of within-year data sets for PSRE hosts:
# All PSRE dissected in 2012, all parasites, sites with >=10 dissected PSRE

# Parition dataset:
indiv_PSRE_2012 <- subset(indiv_PSRE, 
                          substr(site.yr, nchar(site.yr)-3, nchar(site.yr)) %in% "2012")
head(indiv_PSRE_2012)

# Need to remove headings and first column
indiv_PSRE_2012 <- indiv_PSRE_2012[, -c(1:8)]
colnames(indiv_PSRE_2012) <- NULL
head(indiv_PSRE_2012)

# Turn into a matrix
indiv_PSRE_2012 <- as.matrix(indiv_PSRE_2012)

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_PSRE_2012)){
  if(sum(indiv_PSRE_2012[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_PSRE_2012 <- indiv_PSRE_2012[, -zeros]

# Check row sums as well
zeros <- NULL
for (j in 1:nrow(indiv_PSRE_2012)){
  if(sum(indiv_PSRE_2012[ j ,]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
indiv_PSRE_2012 <- indiv_PSRE_2012[-zeros ,]

ncol(indiv_PSRE_2012) #10 parasite species
nrow(indiv_PSRE_2012) #245 hosts


#### METACOMMUNITY ANALYSES #####

# Run the gamet of metacommunity analyses to determine structure
indiv_PSRE_2012_metacom <- metacommunity(indiv_PSRE_2012, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_indiv_PSRE_2012 <- indiv_PSRE_2012_metacom$Comm
quartz(height=6, width=3)
print(Matrix_Plot(ordinated_indiv_PSRE_2012))
# Matrix looks nested w/ two compartments?

indiv_PSRE_2012_metacom$Coherence
#p-val = .0061, positive coherence

indiv_PSRE_2012_metacom$Turnover
#p-val = 0.42, no significant turnover (quasi-structure)
#Turnover slighly lower than the mean, so quasi-nested

indiv_PSRE_2012_metacom$Boundary
#index=3.07, p<0.0001
#Clumped species loss

########################################
# First compartment: rows 1-112
indiv_PSRE_compA2012 <- ordinated_indiv_PSRE_2012[1:112 ,]

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_PSRE_compA2012)){
  if(sum(indiv_PSRE_compA2012[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_PSRE_compA2012 <- indiv_PSRE_compA2012[, -zeros]


# Run the gamet of metacommunity analyses to determine structure
indiv_PSRE_compA2012_metacom <- metacommunity(indiv_PSRE_compA2012, order=F, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_indiv_PSRE_compA2012_metacom <- indiv_PSRE_compA2012_metacom$Comm
quartz(height=3, width=3)
print(Matrix_Plot(ordinated_indiv_PSRE_compA2012_metacom))
# Matrix looks nested w/ two compartments?

indiv_PSRE_compA2012_metacom$Coherence
#p-val = .059, positive coherence

indiv_PSRE_compA2012_metacom$Turnover
#p-val = 0.13, no significant turnover (quasi-structure)
#Turnover lower than the mean, so quasi-nested

indiv_PSRE_compA2012_metacom$Boundary
#index=1.62, p<0.0001
#Clumped species loss

########################################
# Second compartment: rows 113-245
indiv_PSRE_compB2012 <- ordinated_indiv_PSRE_2012[113:245 ,]

# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(indiv_PSRE_compB2012)){
  if(sum(indiv_PSRE_compB2012[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
indiv_PSRE_compB2012 <- indiv_PSRE_compB2012[, -zeros]


# Run the gamet of metacommunity analyses to determine structure
indiv_PSRE_compB2012_metacom <- metacommunity(indiv_PSRE_compB2012, order=F, method="r1", sims=1000)

# Store and view the ordinated occupancy matrix
ordinated_indiv_PSRE_compB2012_metacom <- indiv_PSRE_compB2012_metacom$Comm
quartz(height=3, width=3)
print(Matrix_Plot(ordinated_indiv_PSRE_compB2012_metacom))
# Matrix looks nested w/ two compBrtments?

indiv_PSRE_compB2012_metacom$Coherence
#p-val = .024, positive coherence

indiv_PSRE_compB2012_metacom$Turnover
#p-val = 0.90, no significant turnover (quasi-structure)
#Turnover lower than the mean, so quasi-nested

indiv_PSRE_compB2012_metacom$Boundary
#index=1.64, p<0.0001
#Clumped species loss

########################################
# Save the site scores of the ordinated matrix 
scores_PSRE_2012 <- OrderMatrix(ordinated_indiv_PSRE_2012, output.scores=T)
sitescores_PSRE_2012 <- scores_PSRE_2012$sitescores
sitescores_PSRE_2012 <- melt(sitescores_PSRE_2012)
HostID <- rownames(sitescores_PSRE_2012)
sitescores_PSRE_2012 <- data.frame(HostID, sitescores_PSRE_2012)
head(sitescores_PSRE_2012)
