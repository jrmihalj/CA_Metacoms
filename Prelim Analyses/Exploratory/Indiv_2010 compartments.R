# This code will analyze the two separate metacommunity compartments within the 
# individual PSRE 2010 dataset. 
# August 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# First, I plucked out which HostID's are associated with each compartment
# First (upper) compartment:
# file = '2010_First_comp.csv'
comp1_2010 <- read.csv(file.choose(), header=T)
head(comp1_2010)

# Save as matrix
comp1_2010_Comm <- as.matrix(comp1_2010[, -1])

ncol(comp1_2010_Comm) #13 parasite species
nrow(comp1_2010_Comm) #444 individuals

#### METACOMMUNITY ANALYSES #####
library(metacom)
# Run the gamet of metacommunity analyses to determine structure
comp1_2010_metacom <- metacommunity(comp1_2010_Comm, method="r1", sims=1000,
                                    order=F) #DON'T ORDER, BECAUSE IT ALREADY IS

# Store and view the ordinated occupancy matrix
ordinated_comp1_2010_Comm <- comp1_2010_metacom$Comm
comp1_plot <- Matrix_Plot(ordinated_comp1_2010_Comm)
quartz(height=3, width=3)
print(comp1_plot)

# Coherence
comp1_2010_metacom$Coherence
# p-val=0.059, very close to positive coherence

comp1_2010_metacom$Turnover
# p-val=0.10, no turnover, but more negative (way lower than the mean, but high variance)

comp1_2010_metacom$Boundary
# index=4.29, very clumped
# So, quasi-nested with clumped losses, but coherence not technically significant


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# Second (lower) compartment:
# file = '2010_Second_comp.csv'
comp2_2010 <- read.csv(file.choose(), header=T)
head(comp2_2010)

# Save as matrix
comp2_2010_Comm <- as.matrix(comp2_2010[, -1])

# Get rid of columns (parasite species) with no presence
# Check if any column sums = 0 and remove
zeros <- NULL
for (j in 1:ncol(comp2_2010_Comm)){
  if(sum(comp2_2010_Comm[, j]) == 0){
    zeros <- c(zeros, j)
  }
}
zeros
# Get rid of zeros
comp2_2010_Comm <- comp2_2010_Comm[, -zeros]

ncol(comp2_2010_Comm) #6 parasite species
nrow(comp2_2010_Comm) #351 individuals

#### METACOMMUNITY ANALYSES #####
library(metacom)
# Run the gamet of metacommunity analyses to determine structure
comp2_2010_metacom <- metacommunity(comp2_2010_Comm, method="r1", sims=1000,
                                    order=F) # Already ordered

# Store and view the ordinated occupancy matrix
ordinated_comp2_2010_Comm <- comp2_2010_metacom$Comm
comp2_plot <- Matrix_Plot(ordinated_comp2_2010_Comm)
quartz(height=3, width=3)
print(comp2_plot)

# Coherence
comp2_2010_metacom$Coherence
# p-val=0.022, positive coherence

comp2_2010_metacom$Turnover
# p-val=0.13, no turnover (but lower than mean, high variance)

comp2_2010_metacom$Boundary
# index=1.18, p<0.0001, clumped
# So quasi-nested, with clumped losses

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Prelim analyses:

# subset out the sites that have really large influence?
wo_large <- subset(alldata_indiv_2010, value < 5)

# Pairs plot
pairs(~value+no_malf+SVL+Lat+Long+Elev+Slope+Aspect+hydro+AGG+DEV+FOR+SSG+WET, data=wo_large)
# no_malf looks pretty interesting
pairs(~value+area+depth+pH+tempW+Troph+veg_s+fish, data=wo_large)
pairs(~value+Cond+TDS+Salin+DO+DOmg+DOC+NH4+DON+PO4+DOP+NO3, data=wo_large)
pairs(~value+TATO+RADR+BUBO+AMCA+RACA+PSRE+PHySA+LyMN+HELI+GYRA+RADIX, data=wo_large)
# TATO and PSRE look interesting

###### need to look at total density/richness of hosts at a site! ######

wo_large$no_malf <- as.numeric(as.character(wo_large$no_malf))
plot(value~no_malf, type="p", data=wo_large)
#Definitely looks like there's a pattern with number of malformations. 

wo_large$per_malf <- as.numeric(as.character(wo_large$per_malf))
#Slight pattern with percent malformation
plot(value~per_malf, data=wo_large)

boxplot(value~Malf, data=wo_large)
#Definitely seems like a pattern here

plot(value~TATO, data=wo_large)
#Maybe a pattern here, maybe driven by an outlier?

plot(value~PSRE, data=wo_large)
#Potential pattern here as well. 


# Look at the separate compartments
# Create a new column for Compartment
alldata_indiv_2010$Comp <- rep(0, nrow(alldata_indiv_2010))

# Figure out which hosts belong in which compartment
for (i in 1:nrow(alldata_indiv_2010)){
  out <- NULL
  out <- comp1_2010$ID %in% paste(alldata_indiv_2010$HostID[i]) 
  out <- as.numeric(out)
  
  if(sum(out)==1){
    alldata_indiv_2010$Comp[i] <- 1
  } else {
    alldata_indiv_2010$Comp[i] <- 2
  }
}
# subset out the sites that have really large influence?
wo_large <- subset(alldata_indiv_2010, value < 5)


alldata_indiv_2010$Comp <- as.numeric(as.character(alldata_indiv_2010$Comp))
plot(no_malf~factor(Comp), data=wo_large)
quartz(height=3.5, width=5.5)
plot(value~no_malf, data=subset(wo_large, Comp==1))
plot(value~no_malf, data=subset(wo_large, Comp==2))

# Values are lowest for Comp2, because it is on the bottom of the ordinated matrix

#Replace the NAs in no_malf with zeros
wo_large$no_malf[is.na(wo_large$no_malf)] <- 0
plot(as.numeric(as.character(no_malf))~factor(Comp), data=wo_large)
ggplot(wo_large, aes(x=no_malf, fill=factor(Comp)))+
  geom_histogram(position="dodge", binwidth=0.5)

# Pairs plot
pairs(~Comp+no_malf+SVL+Lat+Long+Elev+Slope+Aspect+hydro+AGG+DEV+FOR+SSG+WET, data=wo_large)
# no_malf looks pretty interesting
pairs(~Comp+area+depth+pH+tempW+Troph+veg_s+fish, data=wo_large)
pairs(~Comp+Cond+TDS+Salin+DO+DOmg+DOC+NH4+DON+PO4+DOP+NO3, data=wo_large)
pairs(~Comp+TATO+RADR+BUBO+AMCA+RACA+PSRE+PHySA+LyMN+HELI+GYRA+RADIX, data=wo_large)











##########################
##########################
# UNUSED CODE:

# comp1_2010_Comm <- subset(indiv_allhosts, HostID %in% comp1_2010$First_comp)
# 
# nrow(comp1_2010_Comm) #444 hosts
# 
# # Clean up this data.frame to only include parasite information
# comp1_2010_Comm <- comp1_2010_Comm[, -c(1:8, 28)]
# head(comp1_2010_Comm)
# 
# # Get rid of columns (parasite species) with no presence
# # Check if any column sums = 0 and remove
# zeros <- NULL
# for (j in 1:ncol(comp1_2010_Comm)){
#   if(sum(comp1_2010_Comm[, j]) == 0){
#     zeros <- c(zeros, j)
#   }
# }
# zeros
# # Get rid of zeros
# comp1_2010_Comm <- comp1_2010_Comm[, -zeros]
# 
# # Subset full dataset to include only these HostID's:
# indiv_allhosts$HostID <- rownames(indiv_allhosts)
# comp2_2010_Comm <- subset(indiv_allhosts, HostID %in% comp2_2010$Second_comp)
# 
# nrow(comp2_2010_Comm) #351 hosts
# 
# # Clean up this data.frame to only include parasite information
# comp2_2010_Comm <- comp2_2010_Comm[, -c(1:8, 28)]
# head(comp2_2010_Comm)
# 
# # Get rid of columns (parasite species) with no presence
# # Check if any column sums = 0 and remove
# zeros <- NULL
# for (j in 1:ncol(comp2_2010_Comm)){
#   if(sum(comp2_2010_Comm[, j]) == 0){
#     zeros <- c(zeros, j)
#   }
# }
# zeros
# # Get rid of zeros
# comp2_2010_Comm <- comp2_2010_Comm[, -zeros]
# 
# 
