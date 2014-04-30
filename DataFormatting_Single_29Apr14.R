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

# Subset out the sites that were sampled in more than one year
# AND had at least 5 PSRE dissected each time

years <- c("2009", "2010", "2011", "2012")

PSRE.multiyear <- data.frame()
for(i in 1:length(unique(PSRE$site))){
  for(t in 1:length(years)){
    sub <- NULL
    sub <- subset(PSRE, site==unique(PSRE$site)[i] & year==years[t])
    if(nrow(sub) >= 5){
      PSRE.multiyear <- rbind(PSRE.multiyear, sub)
    }
  }
}

nrow(PSRE.multiyear) #3048 records

# Figure out which sites were sampled in at least 2 years:
# file name: SitesAndYears.csv

# # Store the sites that were sampled more than one year
# PSRE.included <- data.frame()
# for(i in 1:length(unique(PSRE.multiyear$site))){
#   sub <- NULL
#   sub <- subset(PSRE.multiyear, site==unique(PSRE.multiyear$site)[i])
#   if(length(unique(sub$year)) > 1){ # sampled in more than one year
#     PSRE.included <- rbind(PSRE.included, sub)
#   }
# }

# unique(PSRE.included$site) #69 sites

# Remove all but relevant columns from PSRE.multiyear
PSRE.included <- PSRE.multiyear[, c(1, 3:13)]
unique(PSRE.included$site.yr) # 289 sites.yr's

# # Remove any individuals that have zero parasites
# zero.sums <- NULL
# for(i in 1:nrow(PSRE.included)){
#   if(sum(PSRE.included[i, 1:8])==0){
#     zero.sums <- c(zero.sums, i)
#   }
# }
# length(zero.sums) #379 individuals
# 
# PSRE.included <- PSRE.included[-zero.sums, ]
# nrow(PSRE.included) #1635
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

# Make a site.yr column in Cov:
Cov$site.yr <- paste(Cov$Site, "_", Cov$Year, sep="")

# Subset all the relevant sites:
Cov.multiyear <- data.frame()
for(i in 1:length(unique(PSRE.included$site.yr))){
  sub <- NULL
  sub <- subset(Cov, site.yr==unique(PSRE.included$site.yr)[i])
  Cov.multiyear <- rbind(Cov.multiyear, sub)
}

#Remove "Site" and "Year" columns
Cov.multiyear <- Cov.multiyear[, -c(1:2)]

# Store values into a site x covariate x year array to match Y.obs:
X <- NULL
X <- array(0, dim=c(length(unique(PSRE.included$site.yr)), ncol(Cov.multiyear)-1))
colnames(X) <- colnames(Cov.multiyear)[1:ncol(Cov.multiyear)-1]
rownames(X) <- unique(PSRE.included$site.yr)

for(i in 1:length(unique(PSRE.included$site.yr))){
  sub <- NULL
  sub <- subset(Cov.multiyear, site.yr==unique(PSRE.included$site.yr)[i], select=Elev:hydro)
  for(j in 1:ncol(X)){
    X[i, j] <- sub[1, j]
  }
}

# Standardize each covariate: (value - mean) / 2sd
Xmeans <- NULL
Xsds <- NULL
for(j in 1:12){ # all the non-factor level covariates
  Xmeans[j] <- mean(X[, j], na.rm=T)
  Xsds[j] <- 2*sqrt(var(X[, j], na.rm=T))
  for(i in 1:nrow(X)){
    X[i, j] <- (X[i, j] - Xmeans[j]) / Xsds[j]
  }
}
# X[which(X=="NaN")] <- NA

# Look for collinearity:
cor(X[, c(1:12)], use="complete.obs")
quartz(height=10, width=10)
pairs(X[,c(1:12)])

# Conductivity and Salinity very correlated: Remove Conductivity
X <- X[, -c(6)]

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
               #T=length(years),
               ncovs=ncol(X))

# Set initial parameters:
# Z values (unobserved), mu.x
zinit <- NULL
zinit <- ifelse(Y.obs > 0, 1, 0)

# Start the model
params <- c("b0", "b", "p", "psiMean", "pMean")

# nstore<-1000
# nadap<-1000
# nburn<-5000
# thin<-10
# 
# mod <- NULL
# mod <- jags.model(file = "OccMod_Single_SSVS.txt", 
#                   data = jags_d, n.chains = 3, n.adapt=nadap,
#                   inits = list(z=zinit))
# update(mod, n.iter=nburn)
# 
# out <- NULL
# out <- coda.samples(mod, n.iter = nstore*thin, variable.names = params, thin=thin)


# Run model in parallel:
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

jags.parsamps <- NULL
jags.parsamps <- foreach(i=1:getDoParWorkers()) %dopar% {
  library(rjags)
  setwd("~/GitHub/CA_Metacoms")
  nstore<-1000
  nadap<-30000
  nburn<-30000
  thin<-50
  mod <- jags.model(file = "OccMod_Single_SSVS.txt", 
                    data = jags_d, n.chains = 1, n.adapt=nadap,
                    inits = list(z=zinit))
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = nstore*thin, variable.names = params, thin=thin)
  return(out)
}

bundle <- NULL
bundle <- list(jags.parsamps[[1]][[1]],
               jags.parsamps[[2]][[1]],
               jags.parsamps[[3]][[1]])

class(bundle) <- "mcmc.list"

stopCluster(cl)

# plot(bundle[, "p[8,1]"])
#Plots:
library(mcmcplots)

#mcmcplot(bundle, parms="b") # creates HTML of diagnostic plots

lablims <- array(0, dim=c(7, 2))
start <- 1
for(i in 1:nrow(lablims)){
  lablims[i,] <- c(start, start+42)
  start <- start+43
}

x11(height=4, width=11)
caterplot(bundle, parms="psiMean", horizontal=F)

for(i in 1:nrow(lablims)){
  x11(height=4, width=11)
  caterplot(bundle, parms="b", 
            lab.lim=lablims[i, ], 
            horizontal=F)
}

# for(i in 1:nrow(lablims)){
#   x11(height=4, width=11)
#   caterplot(bundle, parms="c", 
#             lab.lim=lablims[i, ], 
#             horizontal=F)
# }
# 
# for(i in 1:nrow(lablims)){
#   x11(height=4, width=11)
#   caterplot(bundle, parms="d", 
#             lab.lim=lablims[i, ], 
#             horizontal=F)
# }

# Make some data frames...
library(ggmcmc)

post.z <- ggs(bundle, family="z")
post.c <- ggs(bundle, family="c")
post.d <- ggs(bundle, family="d")
post.b <- ggs(bundle, family="b")




