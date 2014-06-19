# California PSRE parasite metacommunities
# Using occupancy modeling in conjunction with the EMS framework

# Script Author: JR Mihaljevic
# Project co-authors: Bethany J. Hoye, Pieter T.J. Johnson

# Script 1: 
# Formatting data set for a single time point analyses.
# This will use all site-years, in order to train the model to better estimate covariate effects
# on occurrence probability. 

# Goals:
# 1. Create a parasite species-by-site array:
#    - Years of interest: 2009 - 2012
#    - I will use all sites that were sampled to estimate covariate effects
#    - Include sites that had at least 5 PSRE dissected each year.
#    - Create vector of J (resampling surveys) per site (i.e. number of frogs sampled)
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
# 2. Create a site-by-covariate array:
#    - Normalize covariate values: (Value - mean)/(2*st.dev.)
#
##############################################################################################################
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################

# 1. Create a parasite species-by-site-by-year array:

# Step 1:
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
PSRE <- PSRE[-unique(remove), ]
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

# Subset out the sites that had at least 8 PSRE dissected each time

years <- c("2009", "2010", "2011", "2012")

PSRE.multiyear <- data.frame()
for(i in 1:length(unique(PSRE$site))){
  for(t in 1:length(years)){
    sub <- NULL
    sub <- subset(PSRE, site==unique(PSRE$site)[i] & year==years[t])
    if(nrow(sub) >= 8){
      PSRE.multiyear <- rbind(PSRE.multiyear, sub)
    }
  }
}

nrow(PSRE.multiyear) # 2939 records

# Sort by year:
PSRE.multiyear <- PSRE.multiyear[order(PSRE.multiyear$year), ]


# Remove all but relevant columns from PSRE.multiyear
PSRE.included <- PSRE.multiyear[, c(1, 3:(ncol(PSRE)-2))]
unique(PSRE.included$site.yr) # 271 sites.yr's

#------------------------------------------------------------------------------------------------------------#

# Step 3:
# Assemble the occurrences into a matrix for each year

Y.obs <- array(0, dim=c(ncol(PSRE.included)-1, length(unique(PSRE.included$site.yr)))) # species, sites, years
colnames(Y.obs) <- unique(PSRE.included$site.yr)
rownames(Y.obs) <- colnames(PSRE.included)[2:ncol(PSRE.included)]

# Create a vector of J:
J <- NULL # the number of individuals dissected at each site

for (i in 1:length(unique(PSRE.included$site.yr))){ # number of sites:
  sub <- NULL
  sub <- subset(PSRE.included, site.yr==unique(PSRE.included$site.yr)[i])
  J[i] <- nrow(sub)
  
  if(nrow(sub)==0){
    Y.obs[, i] <- NA
  }else{
    for(j in 1:(ncol(PSRE.included)-1)){ # number of species:
      Y.obs[j, i] <- sum(sub[, j+1]) 
    } 
  }
}

# Need to remove some sites that have too many missing covariates:
# BigPR_2009 (4), CA-CORTE_2009 (19),
# GDPND013_2010 (117), GDPND015_2010 (119)
# NTalkGCP_2012 (254)
Y.obs <- Y.obs[, -c(4,19,117,119,254)]
J <- J[-c(4,19,117,119,254)]
##############################################################################################################
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################

# 2. Create a site-by-covariate-by-year array:

# Import covariate file:
# Covariates2.csv

Cov <- read.csv(file.choose(), header=T)

# I've used the following covariates:
colnames(Cov)
# [1] "Site"       "Year"       "Lat"        "Long"       "Elev"      
# [6] "Slope"      "Aspect"     "FOR"        "SSG"        "area"      
# [11] "veg_s"      "OpenW"      "fish"       "Cond"       "TDS"       
# [16] "DOmg"       "Amph_Rich"  "Snail_Rich" "AMCA"       "BUBO"      
# [21] "PSRE"       "RACA"       "RADR"       "TATO"       "PHySA"     
# [26] "LyMN"       "HELI"       "GyRA"       "RADIX"      "hydro"

# Make hydro into dummy variables:
levels(Cov$hydro)[1] <- 1 # Permanent
levels(Cov$hydro)[2] <- 2 # Semi-permanent
levels(Cov$hydro)[3] <- 3 # Temporary

# Make a site.yr column in Cov:
Cov$site.yr <- paste(Cov$Site, "_", Cov$Year, sep="")

# Subset all the relevant sites:
Cov.multiyear <- data.frame()
for(i in 1:ncol(Y.obs)){
  sub <- NULL
  sub <- subset(Cov, site.yr==colnames(Y.obs)[i])
  Cov.multiyear <- rbind(Cov.multiyear, sub)
}

#Remove "Site" and "Year" columns
Cov.multiyear <- Cov.multiyear[, -c(1:2)]

# Store values into a site x covariate x year array to match Y.obs:
X <- NULL
X <- array(0, dim=c(ncol(Y.obs), ncol(Cov.multiyear)-1))
colnames(X) <- colnames(Cov.multiyear)[1:ncol(Cov.multiyear)-1]
rownames(X) <- colnames(Y.obs)

for(i in 1:ncol(Y.obs)){
  sub <- NULL
  sub <- subset(Cov.multiyear, site.yr==colnames(Y.obs)[i], select=Lat:hydro)
  
  for(j in 1:ncol(X)){
    X[i, j] <- sub[1, j]
  }
}

# DATA CLEAN UP:

# Need to log transform:
# Area (8), OpenW (10), Cond (12), TDS (13), DOmg (14)
# Square-root trans: Elev (3), Slope(4), Amph_Rich (15),

need.log <- c(8,10,12,13,14)
need.sqrt <- c(3,4,15)

for(i in need.log){X[,i] <- log(X[,i]+1)}
for(i in need.sqrt){X[,i] <- sqrt(X[,i])}


# Insert the NMDS results:
X <- data.frame(X, Amph_NMDS_fixed[,2:3], Snails_NMDS_fixed[,2:3])

# Insert the reciprocal averaging results:
X <- data.frame(X, Amph_RA, Snails_RA_fixed[,2:3])

# Remove the pres/abs of hosts:
X <- X[, -c(17:27)]

# Rearrange a bit:
X <- X[, c(1:10, 12:16, 18:25, 11, 17)]

# Finalized Covariates:
colnames(X)
# [1] "Lat"         "Long"        "Elev"        "Slope"       "Aspect"     
# [6] "FOR"         "SSG"         "area"        "veg_s"       "OpenW"      
# [11] "Cond"        "TDS"         "DOmg"        "Amph_Rich"   "Snail_Rich" 
# [16] "Amph_MDS1"   "Amph_MDS2"   "Snails_MDS1" "Snails_MDS2" "Amph_RA1"   
# [21] "Amph_RA2"    "Snails_RA1"  "Snails_RA2"  "fish"        "hydro"

##############################################################################################################
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################

# Format for Bayesian model:
# First Separate by year:

Xcov_2009 <- as.matrix(X[1:77, ])
Xcov_2010 <- as.matrix(X[78:175, ])
Xcov_2011 <- as.matrix(X[176:236, ])
Xcov_2012 <- as.matrix(X[237:266, ])
Xcov_all <- as.matrix(X)

Yobs_2009 <- Y.obs[, 1:77] # Remove Clin, Fib, Thic
Yobs_2010 <- Y.obs[, 78:175] # ^ Same
Yobs_2011 <- Y.obs[, 176:236] # Remove Clin, Thic
Yobs_2012 <- Y.obs[, 237:266] # Remove Clin, Fib, Thic
Yobs_all <- Y.obs #Remove Thic

Yobs_2009 <- Yobs_2009[-c(2,3,8), ]
Yobs_2010 <- Yobs_2010[-c(2,3,8), ]
Yobs_2011 <- Yobs_2011[-c(2,8), ]
Yobs_2012 <- Yobs_2012[-c(2,3,8), ]
Yobs_all <- Yobs_all[-8, ]

Nsite_2009 <- nrow(Xcov_2009)
Nsite_2010 <- nrow(Xcov_2010)
Nsite_2011 <- nrow(Xcov_2011)
Nsite_2012 <- nrow(Xcov_2012)
Nsite_all <- nrow(Xcov_all)

Nspecies_2009 <- nrow(Yobs_2009)
Nspecies_2010 <- nrow(Yobs_2010)
Nspecies_2011 <- nrow(Yobs_2011)
Nspecies_2012 <- nrow(Yobs_2012)
Nspecies_all <- nrow(Yobs_all)

# All sites surveyed same # times:
J_2009 <- rep(J[1:77], times=Nspecies_2009*Nsite_2009)
J_2010 <- rep(J[78:175], times=Nspecies_2010*Nsite_2010)
J_2011 <- rep(J[176:236], times=Nspecies_2011*Nsite_2011)
J_2012 <- rep(J[237:266], times=Nspecies_2012*Nsite_2012)
J_all <- rep(J, times=Nspecies_all*Nsite_all)

Ncov <- ncol(X)

# X needs to have repeated covariates for each species, long form

X_2009 <- array(0, dim=c(Nsite_2009*Nspecies_2009, Ncov))
X_2010 <- array(0, dim=c(Nsite_2010*Nspecies_2010, Ncov))
X_2011 <- array(0, dim=c(Nsite_2011*Nspecies_2011, Ncov))
X_2012 <- array(0, dim=c(Nsite_2012*Nspecies_2012, Ncov))
X_all <- array(0, dim=c(Nsite_all*Nspecies_all, Ncov))

# Do the following for each year:
t <- 1; i <- 1
TT <- Nsite_2009
while(i <= Nspecies_2009){
  X_2009[t:TT, ] <- Xcov_2009
  t <- t + Nsite_2009
  TT <- TT + Nsite_2009
  i <- i + 1
}

# Species
Species_2009 <- rep(c(1:Nspecies_2009), each=Nsite_2009)
Species_2010 <- rep(c(1:Nspecies_2010), each=Nsite_2010)
Species_2011 <- rep(c(1:Nspecies_2011), each=Nsite_2011)
Species_2012 <- rep(c(1:Nspecies_2012), each=Nsite_2012)
Species_all <- rep(c(1:Nspecies_all), each=Nsite_all)

# Observations/data:
Y_2009 <- NULL
Y_2010 <- NULL
Y_2011 <- NULL
Y_2012 <- NULL
Y_all <- NULL

# Do this for each year:
for(i in 1:Nspecies_all){
  Y_all <- c(Y_all, Yobs_all[i, ])
}

# Number of total observations
Nobs_2009 <- Nspecies_2009*Nsite_2009
Nobs_2010 <- Nspecies_2010*Nsite_2010
Nobs_2011 <- Nspecies_2011*Nsite_2011
Nobs_2012 <- Nspecies_2012*Nsite_2012
Nobs_all <- Nspecies_all*Nsite_all

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
