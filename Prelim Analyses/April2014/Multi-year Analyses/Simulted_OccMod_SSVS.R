# Simulating the multi-year occupancy model 
# Includes stochastic search variable selection

# Author: JR Mihaljevic
# Date: March 2014

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

# Establish some useful functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

##############################################
#####             PARAMETERS             #####
##############################################

# Basic parameters:

N <- 12  # species
K <- 75  # sites
J <- 4   # resurveys
Time <- 3   # years

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

# Initial occurence probabilities:
psiMean <- 0.60
psiSD <- 0.25
mu_b <- Logit(psiMean)

b0 <- rnorm(N, mu_b, psiSD)

# Detection probabilities:
pMean <- 0.75
pSD <- 0.75
mu_p <- Logit(pMean)

p <- array(0, dim=c(N, Time))
for(t in 1:Time){
  p[, t] <- AntiLogit(rnorm(N, mu_p, pSD))
}

# Persistence probabilities:
phiMean <- 0.80
phiSD <- 0.25
mu_phi <- Logit(phiMean)

c0 <- array(0, dim=c(N, Time-1))
for(t in 1:Time-1){
  c0[, t] <- rnorm(N, mu_phi, phiSD)
}

# Colonization probabilities:
gamMean <- 0.50
gamSD <- 0.25
mu_gam <- Logit(gamMean)

d0 <- array(0, dim=c(N, Time-1))
for(t in 1:Time-1){
  d0[, t] <- rnorm(N, mu_gam, gamSD)
}

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

# Covariate effects, coded in the SSVS framework:
c <- 500
p_incl <- 0.50

bSD <- 0.005
cSD <- 0.015
dSD <- 0.005

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

# Generate covariates

ncov <- 20
X <- array(0, dim=c(K, ncov, Time))

for (t in 1:Time){
  for(j in 1:ncov){
    X[, j, t] <- rnorm(K, 0, 1)
  }
}

included <- array(0, dim=c(N, ncov, 3)) # 3 covariate effect types (b, c, d)

for (i in 1:N){
  for(j in 1:3){ 
    included[i, ,j] <- rbinom(ncov, 1, p_incl)
  }
}

# Generate covariates
b.spp <- array(0, dim=c(N, ncov))
c.spp <- array(0, dim=c(N, ncov))
d.spp <- array(0, dim=c(N, ncov))  

for(i in 1:N){
  for(j in 1:ncov){
    b.spp[i, j] <- rnorm(n=1, mean=0,
                        sd=ifelse(included[i,j,1]==1,
                                  bSD * c,
                                  bSD))
    c.spp[i, j] <- rnorm(n=1, mean=0,
                        sd=ifelse(included[i,j,2]==1,
                                  cSD * c,
                                  cSD))
    d.spp[i, j] <- rnorm(n=1, mean=0,
                        sd=ifelse(included[i,j,3]==1,
                                  dSD * c,
                                  dSD))
  }
}
  

                      
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

#############################################
#####     SIMULATE OCCUPANCY STATES     #####
#############################################

# Storage:
temp.b <- array(0, dim=c(N, K, ncov, Time))
temp.c <- array(0, dim=c(N, Time-1, K, ncov))
temp.d <- array(0, dim=c(N, Time-1, K, ncov))

Y <- array(0, dim=c(N, K, Time))
z <- array(0, dim=c(N, K, Time))


for (i in 1:N) { # n=number of species
  
  ## Species-specific covariate effects on initial occurrence:
  for (cov in 1:ncov) {
    for(k in 1:K){
      temp.b[i, k, cov, 1] <- b.spp[i, cov] * X[k, cov, 1]
    }
  }
  
  for (k in 1:K) {
    lpsi <- NULL
    psi <- NULL
    
    lpsi <- b0[i] + sum(temp.b[i, k, , 1]) # Year-specific covariate values
    psi <- AntiLogit(lpsi)
    z[i,k,1] <- rbinom(1, 1, psi)
    
    Y[i,k,1] <- rbinom(1, J, p[i,1]*z[i,k,1]) # J=number of hosts sampled from each site
  }
  
  ## Model of changes in occupancy state for t=2, ..., T
  
  for (cov in 1:ncov) {
    for(t in 1:(Time-1)){
      for(k in 1:K){
        temp.c[i, t, k, cov] <- c.spp[i, cov] * X[k, cov, t]
        temp.d[i, t, k, cov] <- d.spp[i, cov] * X[k, cov, t]
      }
    }
  }
  
  for (t in 1:(Time-1)) {
    for (k in 1:K) { # K=number of sites
      lgam <- NULL
      gam <- NULL
      lphi <- NULL
      phi <- NULL
      mu.z <- NULL
      
      lgam <- c0[i,t] + sum(temp.c[i, t, k, ])
      gam <-  AntiLogit(lgam)
      lphi <- d0[i,t] + sum(temp.d[i, t, k, ])
      phi <-  AntiLogit(lphi)
      
      mu.z <- (phi * z[i,k,t] + gam * (1-z[i,k,t]))
      z[i,k,t+1] <- rbinom(1, 1, mu.z)
      
      Y[i,k,t+1] <- rbinom(1, J, p[i,t+1]*z[i,k,t+1])
    } 
  }
}

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

##################################
#####     RUN JAGS MODEL     #####
##################################
library(rjags)
jags_d <- NULL
jags_d <- list(x=X,
               y=Y,
               K=K,
               n=N,
               J=J,
               T=Time,
               ncovs=ncov)

# Set initial parameters:
# Z values (unobserved)
zinit <- NULL
zinit <- ifelse(Y > 0, 1, 0)

# Start the model
params <- c("b", "c", "d", "p")

mod <- NULL
mod <- jags.model(file = "OccMod_MultiYear_SSVS.txt", 
                  data = jags_d, n.chains = 3, n.adapt=1000,
                  inits = list(z=zinit))
update(mod, n.iter=5000) # 5000 burn-in

out <- NULL
out <- coda.samples(mod, n.iter = 5000, variable.names = params, thin=5)


library(mcmcplots)
quartz(height=4, width=11)
caterplot(out, parms="p", horizontal=F)
caterpoints(as.vector(p), horizontal=F)

quartz(height=4, width=11)
caterplot(out, parms="b", horizontal=F) #use lab.lim option to plot fewer parameters
caterpoints(as.vector(b.spp), horizontal=F)

quartz(height=4, width=11)
caterplot(out, parms="c", horizontal=F)
caterpoints(as.vector(c.spp), horizontal=F)

quartz(height=4, width=11)
caterplot(out, parms="d", horizontal=F)
caterpoints(as.vector(d.spp), horizontal=F)

