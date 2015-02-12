#BEAST R-code version
#Data: 2009

# First clear R's memory, so nodes don't have same seed:
rm(list=ls()) 
gc()

# Set working directory:
setwd("/home/jmihaljevic/PhD/CAmetacoms/")
# Load in the FormmattedData.RData:
load("FormattedData.RData")

# Some more data formatting:
colnames(X_2009) <- colnames(X)

# Standardize each covariate: (value - mean) / 2sd
for(j in 1:15){ # all the non-factor level covariates
  X_2009[, j] <- (X_2009[, j] - mean(X_2009[, j], na.rm=T)) / (2 * sd(X_2009[, j], na.rm=T))
}

#Collinearity:
#Lat-long 
#Long-Elev (quadratic)
#FOR-SSG
#Cond-TDS
#Amph_RA1 - SSG
#Amph_RA1 - Amph_MDS1
#Amph_Rich - Amph_RA2


# REMOVE: Long, Elev, SSG, FOR, Amph_MDS1, Amph_RA2, TDS

X_2009 <- X_2009[, -c(2,3,6,7,12,16,21)]

# Remove Fish:
X_2009 <- X_2009[, -17]

Ncov_2009 <- ncol(X_2009)
################################################
# JAGS model run:
################################################
library(rjags)
library(random)
load.module("lecuyer")
load.module("glm")

# Input the data:
Y=Y_2009
X=X_2009
Species=Species_2009
Nspecies=Nspecies_2009
Ncov=Ncov_2009
Nobs=Nobs_2009
J=J_2009

setwd("/home/jmihaljevic/PhD/CAmetacoms/2009/")

# create unique identifier for iteration
uniqueID <- paste(i, "-", runif(1, 0, 1), sep="")
pid <- Sys.getpid()
# make a data file with unique name
dump(c("X", "Y", "Species", "Nspecies", "Ncov", "Nobs", "J"), 
     file=paste(uniqueID, "data.R", sep="_"))

# Set initial parameters:
# Z values (unobserved)
zinit <- NULL
zinit <- ifelse(Y_2009 > 0, 1, 0)
z = zinit

# Write a file of the inits:
dump("z", file=paste(uniqueID, "inits.R", sep="_"))

# write a unique jags.cmd file that refers to the data and initial values
cat(paste('model in "MLMmodel2009_all.txt"
          data in "', uniqueID, '_data.R"
          compile, nchains(1)
          parameters in "', uniqueID, '_inits.R"
          initialize
          adapt 100
          update 10
          monitor alpha, thin(10)
          monitor betas, thin(10)
          monitor p.detect, thin(10)
          monitor sd.beta.post, thin(10)
          monitor mean.beta.post, thin(10)
          monitor z, thin(10)
          update 10
          coda *, stem(PID', pid, ')
          exit', 
          sep=""), 
    file=paste(uniqueID, "jags.cmd", sep=""))

# Run the JAGS script file:
system(paste("jags ", uniqueID, "jags.cmd", sep=""))

# Delete all newly created files and directory:
system(paste("rm ", uniqueID, "*", sep=""))