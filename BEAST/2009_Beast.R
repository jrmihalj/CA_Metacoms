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
# data:
jags_d <- list(Y=Y_2009,
               X=X_2009,
               Species=Species_2009,
               Nspecies=Nspecies_2009,
               Ncov=Ncov_2009,
               Nobs=Nobs_2009,
               J=J_2009)

# parameters:
params <- c("alpha", "betas", "p.detect", 
            "p.include", "sd.beta.post", "mean.beta.post",
            "z")

jinits <- function() {
  list(
    z=ifelse(Y_2009 > 0, 1, 0),
    .RNG.name=c("base::Super-Duper"),
    .RNG.seed=as.numeric(randomNumbers(n = 1, min = 1, max = 1e+06, col=1))
  )
}

JAGS_func <- function(store, adapt, burn, thin){
  mod <- jags.model(file = "MLMmodel2009_all.txt",
                    data = jags_d, n.chains = 1, n.adapt=adapt,
                    inits = jinits)
  update(mod, n.iter=burn)
  out <- coda.samples(mod, n.iter = store*thin, variable.names = params, thin=thin)
  return(out)
}

#Change wd so that JAGS can find model spec
setwd("/home/jmihaljevic/PhD/CAmetacoms/2009")

jags.parsamps <- NULL
jags.parsamps <- JAGS_func(store=100, adapt=100, burn=5, thin=5)

pid <- Sys.getpid()
save(jags.parsamps, file=paste("2009Chain", pid, sep="_"))