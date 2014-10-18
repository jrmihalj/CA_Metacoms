colnames(X_2012) <- colnames(X)

# Standardize each covariate: (value - mean) / 2sd
for(j in 1:15){ # all the non-factor level covariates
  X_2012[, j] <- (X_2012[, j] - mean(X_2012[, j], na.rm=T)) / (2 * sd(X_2012[, j], na.rm=T))
}

# Look for collinearity:
cor(X_2012, use="complete.obs")
x11(height=10, width=15)
pairs(X_2012[, c(1:13,25)])
pairs(X_2012[,c(6:10, 14:23)])

# lots of Collinearity

# REMOVE: Lat, Long, Elev, area, Cond, TDS, veg_s, 
# slope, FOR, SSG, SnailRich, SnailMDS2, SnailRA1, AmphRA1, AmphRA2, AmphMDS2

X_2012 <- X_2012[, -c(1:4,6:9,11,12,15,17:22)]

# REMOVE: Fish

X_2012 <- X_2012[, -7]

Ncov_2012 <- ncol(X_2012)
################################################
# JAGS model run:
################################################
library(rjags)
library(random)
load.module("lecuyer")
load.module("glm")
# data:
jags_d <- list(Y=Y_2012,
               X=X_2012,
               Species=Species_2012,
               Nspecies=Nspecies_2012,
               Ncov=Ncov_2012,
               Nobs=Nobs_2012,
               J=J_2012)

# parameters:
params <- c("alpha", "betas", "p.detect", 
            "p.include", "sd.beta.post", "mean.beta.post",
            "z")

jinits <- function() {
  list(
    z=ifelse(Y_2012 > 0, 1, 0),
    .RNG.name=c("base::Super-Duper"),
    .RNG.seed=as.numeric(randomNumbers(n = 1, min = 1, max = 1e+06, col=1))
  )
}

# initialize model:
# This assumes you can run three parallel chains. Change accordingly.
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

jags.parsamps <- NULL
jags.parsamps <- foreach(i=1:3, .packages=c('rjags','random')) %dopar% {
  #setwd("~/GitHub/MLM_EcologyInSilico")
  store<-1500
  nadap<-80000
  nburn<-80000
  thin<-50
  mod <- jags.model(file = "MLM_model.txt", 
                    data = jags_d, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, variable.names = params, thin=thin)
  return(out)
}

bundle_full <- NULL
bundle_full <- list(jags.parsamps[[1]][[1]],
                    jags.parsamps[[2]][[1]],
                    jags.parsamps[[3]][[1]])

class(bundle_full) <- "mcmc.list"

stopCluster(cl)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Check Convergence:
################################################
library(ggmcmc)
library(mcmcplots)
alpha.df <- ggs(bundle_full, family="alpha")
beta.df <- ggs(bundle_full, family="betas")
sd.beta.df <- ggs(bundle_full, family="sd.beta.post")
mean.beta.df <- ggs(bundle_full, family="mean.beta.post")
p.detect.df <- ggs(bundle_full, family="p.detect")

ggs_Rhat(alpha.df)
ggs_Rhat(beta.df)
ggs_Rhat(sd.beta.df)
ggs_Rhat(p.detect.df)
ggs_Rhat(mean.beta.df)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(bundle_full, parms="betas", horizontal=F)#, random=50)

caterplot(bundle_full, parms="mean.beta.post", horizontal=F)
caterplot(bundle_full, parms="sd.beta.post", horizontal=F)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Calculate WAIC for full model:
# ################################################
source(file="calc_waic.R")

WAIC_full <- calc_waic(bundle_full, jags_d)
WAIC_full$WAIC
# 517.0422

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Check 'best' model:
# ################################################

# Check to see which covariate effects (slopes) != 0

source(file="HDI.R")
hdi.mean <- array(0, dim=c(Ncov_2012, 2))

for(i in 1:Ncov_2012){
  sub <- subset(mean.beta.df, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean[i, ] <- hdi
}
hdi.mean

# Which do not include zero?
# OpenW (2), Amph_Rich(4), hydro (7)

# Check to see which covariate effects are significantly random (st.dev > 0)

hdi.sd <- array(0, dim=c(Ncov_2012, 2))

for(i in 1:Ncov_2012){
  sub <- subset(sd.beta.df, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd[i, ] <- hdi
}
hdi.sd

# Which are greater than zero?
# Aspect (1)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now refit reduced model:
################################################

# X_reduced for reduced covariates:

X_2012_reduced <- X_2012[, c(1,2,4,7)]
Ncov_2012_reduced <- ncol(X_2012_reduced)

# data:
jags_d_reduced <- list(Y=Y_2012,
                       X=X_2012_reduced,
                       Species=Species_2012,
                       Nspecies=Nspecies_2012,
                       Ncov=Ncov_2012_reduced,
                       Nobs=Nobs_2012,
                       J=J_2012)

# parameters:
params_reduced <- c("alpha", "betas", "p.detect", 
                    "sd.beta.post", "mean.beta.post", "psi", 
                    "z")

# initialize model:

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

jags.parsamps <- NULL
jags.parsamps <- foreach(i=1:3, .packages=c('rjags','random')) %dopar% {
  #setwd("C:\Users\Joe\Documents\GitHub\CA_Metacoms")
  store<-1500
  nadap<-80000
  nburn<-80000
  thin<-50
  mod <- jags.model(file = "MLM_model.txt", 
                    data = jags_d_reduced, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, 
                      variable.names = params_reduced, thin=thin)
  return(out)
}

bundle_reduced <- NULL
bundle_reduced <- list(jags.parsamps[[1]][[1]],
                       jags.parsamps[[2]][[1]],
                       jags.parsamps[[3]][[1]])

class(bundle_reduced) <- "mcmc.list"

stopCluster(cl)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Calculate WAIC for reduced model:
# ################################################
source(file="calc_waic.R")

WAIC_reduced <- calc_waic(bundle_reduced, jags_d_reduced)
WAIC_reduced$WAIC
# 498.7703

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Try Dropping Aspect        (3cov): best
# Try Dropping hydro instead (3covB): 
# Try Dropping OpenW instead (3covC): 
# Try Dropping Amph_Rich     (3covD): (no diff from 3cov)
 
# Try Dropping Aspect + OpenW     (2cov): best
# Try Dropping Aspect + hydro     (2covB): 
# Try Dropping Aspect + Amph_Rich (2covC):

# Try only Amph_Rich (1cov): better, but:
# Try only hydro     (1covB): 
################################################

X_2012_2cov <-(X_2012_reduced[, -c(3)])
head(X_2012_2cov)
Ncov_2012_2cov <- ncol(X_2012_2cov)

# data:
jags_d_2cov <- list(Y=Y_2012,
                    X=X_2012_2cov,
                    Species=Species_2012,
                    Nspecies=Nspecies_2012,
                    Ncov=Ncov_2012_2cov,
                    Nobs=Nobs_2012,
                    J=J_2012)


# initialize model:

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

jags.parsamps <- NULL
jags.parsamps <- foreach(i=1:3, .packages=c('rjags','random')) %dopar% {
  #setwd("C:\Users\Joe\Documents\GitHub\CA_Metacoms")
  store<-1500
  nadap<-80000
  nburn<-80000
  thin<-50
  mod <- jags.model(file = "MLM_model.txt", 
                    data = jags_d_2cov, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, 
                      variable.names = params, thin=thin)
  return(out)
}

bundle_2cov <- NULL
bundle_2cov <- list(jags.parsamps[[1]][[1]],
                     jags.parsamps[[2]][[1]],
                     jags.parsamps[[3]][[1]])

class(bundle_2cov) <- "mcmc.list"

stopCluster(cl)

caterplot(bundle_2cov, parms="betas", horizontal=F)
caterplot(bundle_2cov, parms="mean.beta.post", horizontal=F)
caterplot(bundle_2cov, parms="sd.beta.post", horizontal=F)
caterplot(bundle_2cov, parms="alpha", horizontal=F)
caterplot(bundle_2cov, parms="p.detect", horizontal=F)


################################################
# Calculate WAIC for step-wise removal models:
################################################
source(file="calc_waic.R")

WAIC_full <- calc_waic(bundle_full, jags_d)
WAIC_full$WAIC
# 517.0422

WAIC_reduced <- calc_waic(bundle_reduced, jags_d_reduced)
WAIC_reduced$WAIC
# 498.7703

###########
WAIC_3cov <- calc_waic(bundle_3cov, jags_d_3cov)
WAIC_3cov$WAIC
# 494.7796

WAIC_3covB <- calc_waic(bundle_3covB, jags_d_3covB)
WAIC_3covB$WAIC
# 496.8541 worse

WAIC_3covC <- calc_waic(bundle_3covC, jags_d_3covC)
WAIC_3covC$WAIC
# 495.3057 worse

WAIC_3covD <- calc_waic(bundle_3covD, jags_d_3covD)
WAIC_3covD$WAIC
# 494.73 no different than 3cov

###########
WAIC_2cov <- calc_waic(bundle_2cov, jags_d_2cov)
WAIC_2cov$WAIC
# 492.8701** better

WAIC_2covB <- calc_waic(bundle_2covB, jags_d_2covB)
WAIC_2covB$WAIC
# 494.2248 nope...

WAIC_2covC <- calc_waic(bundle_2covC, jags_d_2covC)
WAIC_2covC$WAIC
# 494.046 nope...

###########
WAIC_1cov <- calc_waic(bundle_1cov, jags_d_1cov)
WAIC_1cov$WAIC
# 492.0211

WAIC_1covB <- calc_waic(bundle_1covB, jags_d_1covB)
WAIC_1covB$WAIC
# 491.7456 ***

###########
WAIC_null <- calc_waic(bundle_null, jags_d_null)
WAIC_null$WAIC
# 492.3237

# BEST MODEL:
# Hydro OR 
# Null OR 
# Amph_Rich OR
# Amph_Rich + Hydro

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now fit a null model (only random intercepts):
################################################

# data:
jags_d_null <- list(Y=Y_2012,
                    Species=Species_2012,
                    Nspecies=Nspecies_2012,
                    Nobs=Nobs_2012,
                    J=J_2012)

# parameters:
params_null <- c("alpha", "p.detect", "psi", "z")


# initialize model:

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

jags.parsamps <- NULL
jags.parsamps <- foreach(i=1:3, .packages=c('rjags','random')) %dopar% {
  #setwd("C:\Users\Joe\Documents\GitHub\CA_Metacoms")
  store<-1500
  nadap<-80000
  nburn<-80000
  thin<-50
  mod <- jags.model(file = "MLM_model_Null.txt", 
                    data = jags_d_null, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, 
                      variable.names = params_null, thin=thin)
  return(out)
}

bundle_null <- NULL
bundle_null <- list(jags.parsamps[[1]][[1]],
                    jags.parsamps[[2]][[1]],
                    jags.parsamps[[3]][[1]])

class(bundle_null) <- "mcmc.list"

stopCluster(cl)


# ################################################
# # Calculate WAIC for reduced model:
# ################################################
source(file="calc_waic.R")

WAIC_null <- calc_waic(bundle_null, jags_d_null)
WAIC_null$WAIC
# 492.3237

################################################
# Check random vs. fixed of BEST MODEL:
################################################
mean.beta.df.best <- ggs(bundle_1covB, family="mean.beta.post")
sd.beta.df.best <- ggs(bundle_1covB, family="sd.beta.post")

# Check slopes (Fixed effects)

hdi.mean.best <- HDI(mean.beta.df.best[,4])
hdi.mean.best

# Which do not include zero?
# NONE: No fixed effects

# Check st.dev (Random Effects)

hdi.sd.best <- HDI(sd.beta.df.best[,4])
hdi.sd.best

# No random effects either
