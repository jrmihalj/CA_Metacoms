# California PSRE parasite metacommunities
# Using occupancy modeling in conjunction with the EMS framework

# Script Author: JR Mihaljevic
# Project co-authors: Bethany J. Hoye, Pieter T.J. Johnson

# Script 1: 
# Formatting data set

# Goals:
# 1. Create a parasite species-by-site-by-year array:
#    - Years of interest: 2009 - 2012
#    - Need to use sites that were sampled at least two years in order to estimate colonization
#      and persistence parameters
#    - Create matrix of J (resampling surveys) per site (i.e. number of frogs sampled per year)
#    - Create one array for all parasite species, and one with trematode-only parasites:
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
#
# Gorg, Glyp, Hae, Hali, Mega, Ozwa, Rhab and Spir are all thought to be 
# acquired post-metamorphosis.
#
# 2. Create a site-by-covariate-by-year array:
#    - One for each array above
#    - Normalize covariate values: (Value - mean)/st.dev.
#
##############################################################################################################
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################

# 1. Create a parasite species-by-site-by-year array:

# Step 1:
# Figure out which sites were sampled in at least 2 years:
# file name: SitesAndYears.csv
SitesYears <- read.csv(file.choose(), header=T) #SitesAndYears.csv
levels(SitesYears$site)

# Store the sites that were sampled more than one year
include <- NULL 
for(i in 1:length(levels(SitesYears$site))){
  sub <- NULL
  sub <- subset(SitesYears, site==levels(SitesYears$site)[i])
  
  if(nrow(sub) > 1){
    include <- c(include, levels(SitesYears$site)[i])
  }
}
#------------------------------------------------------------------------------------------------------------#

# Step 2:
# Subset out all the sites sampled more than one year:

# Import the dataset that has all individual PSRE dissected
# file name: 'Indiv_ALL_Host.csv'
PSRE <- read.csv(file.choose(), header=T)
head(PSRE)

# Subset only PSRE:
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

# Remove the unused columns
PSRE <- PSRE[, -c(3:8, rem_cols)]

# Remove non-trematode species:
nontrem <- c(9:11)
PSRE <- PSRE[, -nontrem]


# Make columns for site and year labels:
PSRE$site.yr <- as.character(PSRE$site.yr)
PSRE$site <- substr(PSRE$site.yr, 1, nchar(PSRE$site.yr)-5)
PSRE$year <- substr(PSRE$site.yr, nchar(PSRE$site.yr)-3, nchar(PSRE$site.yr))

# Subset out the sites that were sampled in more than one year:
PSRE.multiyear <- data.frame()
for(i in 1:length(include)){
  sub <- NULL
  sub <- subset(PSRE, site==include[i])
  PSRE.multiyear <- rbind(PSRE.multiyear, sub)
}

nrow(PSRE.multiyear) #2055 records

# Remove all but relevant columns from PSRE.multiyear
PSRE.multiyear <- PSRE.multiyear[, 3:ncol(PSRE.multiyear)]

# Remove any individuals that have zero parasites
zero.sums <- NULL
for(i in 1:nrow(PSRE.multiyear)){
  if(sum(PSRE.multiyear[i, 1:8])==0){
    zero.sums <- c(zero.sums, i)
  }
}
length(zero.sums) #392 individuals

PSRE.multiyear <- PSRE.multiyear[-zero.sums, ]
nrow(PSRE.multiyear) #1663
#------------------------------------------------------------------------------------------------------------#

# Step 3:
# Assemble the occurrences into a matrix for each year

Y.obs <- array(0, dim=c(8, length(include), 4)) # species, sites, years
colnames(Y.obs) <- include
rownames(Y.obs) <- colnames(PSRE.multiyear)[1:8]

# Add up occurrences for each site:
years <- c("2009", "2010", "2011", "2012")

# Create a matrix of J:
J <- array(0, dim=c(length(include), 4)) # the number of individuals dissected at each site

for(t in 1:4){ # number of years:
  for (i in 1:length(include)){ # number of sites:
  sub <- NULL
  sub <- subset(PSRE.multiyear, year==years[t] & site==include[i])
  J[i, t] <- nrow(sub)
  
    for(j in 1:8){ # number of species:
      Y.obs[j, i, t] <- sum(sub[, j]) 
    } 
  }
}

##############################################################################################################
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################

# 2. Create a site-by-covariate-by-year array:

# Import covariate file:
# Covariates.csv

Cov <- read.csv(file.choose(), header=T)

# I've used the following covariates:
colnames(Cov)
# [1] "Site"       "Year"       "Elev"       "WET"       
# [5] "pH"         "veg_s"      "Canopy"     "Cond"      
# [9] "Salin"      "DOmg"       "TotalN"     "TotalP"    
# [13] "Amph_Rich"  "Snail_Rich" "AMCA"       "BUBO"      
# [17] "PSRE"       "RACA"       "RADR"       "TATO"      
# [21] "PHySA"      "LyMN"       "HELI"       "GyRA"      
# [25] "RADIX"      "hydro" 

# Columns 15 - 25 are pres/abs of host sp

# Column 26 is hydro with levels: Permanent, Semi-permanent, Temporary
# Make hydro into dummy variables:
levels(Cov$hydro)[1] <- 1 # Permanent
levels(Cov$hydro)[2] <- 2 # Semi-permanent
levels(Cov$hydro)[3] <- 3 # Temporary

# Subset all the relevant sites:
Cov.multiyear <- data.frame()
for(i in 1:length(include)){
  sub <- NULL
  sub <- subset(Cov, Site==include[i])
  Cov.multiyear <- rbind(Cov.multiyear, sub)
}

# Store values into a site x covariate x year array to match Y.obs:
X <- NULL
X <- array(0, dim=c(length(include), ncol(Cov.multiyear)-2, length(years)))
colnames(X) <- colnames(Cov.multiyear)[3:ncol(Cov.multiyear)]
rownames(X) <- include

for(t in 1:length(years)){
  for(i in 1:length(include)){
    sub <- NULL
    sub <- subset(Cov.multiyear, Site==include[i] & Year==years[t], select=Elev:hydro)
    for(j in 1:ncol(X)){
    X[i, j, t] <- sub[1, j]
    }
  }
}

# Standardize each covariate: (value - mean)
Xmeans <- array(0, dim=c(12, length(years)))
for(t in 1:length(years)){
  for(j in 1:12){ # all the non-factor level covariates
      Xmeans[j, t] <- mean(X[, j, t], na.rm=T)
    for(i in 1:nrow(X)){
      X[i, j, t] <- X[i, j, t] - Xmeans[j, t]
    }
  }
}
X[which(X=="NaN")] <- NA

# Look for collinearity:
cor(X[, c(1:12), 1], use="complete.obs")
quartz(height=10, width=10)
pairs(X[,c(1:12),2])

# Conductivity and Salinity very correlated: Remove Conductivity
# Snail_Rich and TotalN very correlated: Remove TotalN
X <- X[, -c(6,9), ]

# Finalized Covariates:
colnames(X)


##############################################################################################################
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################

##################################
#####     RUN JAGS MODEL     #####
##################################
library(rjags)
jags_d <- NULL
jags_d <- list(x=X,
               y=Y.obs,
               K=ncol(Y.obs),
               n=nrow(Y.obs),
               J=J,
               T=length(years),
               ncovs=ncol(X))

# Set initial parameters:
# Z values (unobserved), mu.x
zinit <- NULL
zinit <- ifelse(Y.obs > 0, 1, 0)

# Start the model
params <- c("b0", "b", "c", "d", "p")

mod <- NULL
mod <- jags.model(file = "OccMod_MultiYear_SSVS.txt", 
                  data = jags_d, n.chains = 3, n.adapt=1000,
                  inits = list(z=zinit))
update(mod, n.iter=5000) # 5000 burn-in

out <- NULL
out <- coda.samples(mod, n.iter = 5000, variable.names = params, thin=5)

#Plots:
library(mcmcplots)

lablims <- array(0, dim=c(7, 2))
start <- 1
for(i in 1:7){
  lablims[i,] <- c(start, start+42)
  start <- start+43
}

quartz(height=4, width=11)
caterplot(out, parms="b0", horizontal=F)

for(i in 1:nrow(lablims)){
  quartz(height=4, width=11)
  caterplot(out, parms="b", 
            lab.lim=lablims[i, ], 
            horizontal=F)
}

for(i in 1:nrow(lablims)){
  quartz(height=4, width=11)
  caterplot(out, parms="c", 
            lab.lim=lablims[i, ], 
            horizontal=F)
}

for(i in 1:nrow(lablims)){
  quartz(height=4, width=11)
  caterplot(out, parms="d", 
            lab.lim=lablims[i, ], 
            horizontal=F)
}

