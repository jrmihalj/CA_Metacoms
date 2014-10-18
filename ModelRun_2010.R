colnames(X_2010) <- colnames(X)

# Standardize each covariate: (value - mean) / 2sd
for(j in 1:15){ # all the non-factor level covariates
  X_2010[, j] <- (X_2010[, j] - mean(X_2010[, j], na.rm=T)) / (2 * sd(X_2010[, j], na.rm=T))
}

# Look for collinearity:
cor(X_2010, use="complete.obs")
x11(height=10, width=15)
pairs(X_2010[,11:23])
pairs(X_2010[,c(6:10, 14:23)])

#Collinearity:
#Lat-long 
#Long-Elev (quadratic)
#FOR-SSG
#Cond-TDS
#Slope-Area
#Snails_RA1 - Snails_MDS1
#Snails_RA1 - Snail_Rich


# REMOVE: Long, Elev, SSG, Slope, Snails_MDS1, Snails_RA1, TDS

X_2010 <- X_2010[, -c(2,3,4,7,12,18,22)]

# Remove fish:
X_2010 <- X_2010[, -17]

Ncov_2010 <- ncol(X_2010)
################################################
# JAGS model run:
################################################
library(rjags)
library(random)
load.module("lecuyer")
load.module("glm")
# data:
jags_d <- list(Y=Y_2010,
               X=X_2010,
               Species=Species_2010,
               Nspecies=Nspecies_2010,
               Ncov=Ncov_2010,
               Nobs=Nobs_2010,
               J=J_2010)

# parameters:
params <- c("alpha", "betas", "p.detect", 
            "p.include", "sd.beta.post", "mean.beta.post",
            "z")

jinits <- function() {
  list(
    z=ifelse(Y_2010 > 0, 1, 0),
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
#ggs_Rhat(p.include.df)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(bundle_full, parms="betas", horizontal=F)#, random=50)

caterplot(bundle_full, parms="sd.beta.post", horizontal=F)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Calculate WAIC for full model:
# ################################################
source(file="calc_waic.R")

WAIC_full <- calc_waic(bundle_full, jags_d)
WAIC_full$WAIC
# 2224.971

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Check best model:
# ################################################

# Check to see which covariate effects (slopes) != 0

source(file="HDI.R")
hdi.mean <- array(0, dim=c(Ncov_2010, 2))

for(i in 1:Ncov_2010){
  sub <- subset(mean.beta.df, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean[i, ] <- hdi
}
hdi.mean

# Which do not include zero?
# Lat(1), Snail_Rich(10), Amph_RA2(15), Snails_RA2(16)

# Check to see which covariate effects are significantly random (st.dev > 0)

hdi.sd <- array(0, dim=c(Ncov_2010, 2))

for(i in 1:Ncov_2010){
  sub <- subset(sd.beta.df, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd[i, ] <- hdi
}
hdi.sd

# Which are greater than zero?
# Snails_MDS2(13)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now refit reduced model:
################################################

# X_reduced for reduced covariates:

X_2010_reduced <- X_2010[, c(1, 10, 13, 15, 16)]
Ncov_2010_reduced <- ncol(X_2010_reduced)

# data:
jags_d_reduced <- list(Y=Y_2010,
                    X=X_2010_reduced,
                    Species=Species_2010,
                    Nspecies=Nspecies_2010,
                    Ncov=Ncov_2010_reduced,
                    Nobs=Nobs_2010,
                    J=J_2010)

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

################################################
# Check Convergence:
################################################
library(ggmcmc)
library(mcmcplots)
alpha.df.reduced <- ggs(bundle_reduced, family="alpha")
beta.df.reduced <- ggs(bundle_reduced, family="betas")
sd.beta.df.reduced <- ggs(bundle_reduced, family="sd.beta.post")
mean.beta.df.reduced <- ggs(bundle_reduced, family="mean.beta.post")
p.detect.df.reduced <- ggs(bundle_reduced, family="p.detect")
psi.df.reduced <- ggs(bundle_reduced, family="psi")
z.df.reduced <- ggs(bundle_reduced, family="z")

ggs_Rhat(alpha.df.reduced)
ggs_Rhat(beta.df.reduced)
ggs_Rhat(sd.beta.df.reduced)
ggs_Rhat(p.detect.df.reduced)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(bundle_reduced, parms="betas", horizontal=F)
caterplot(bundle_reduced, parms="mean.beta.post", horizontal=F)
caterplot(bundle_reduced, parms="sd.beta.post", horizontal=F)

# ################################################
# # Calculate WAIC for reduced model:
# ################################################
source(file="calc_waic.R")

WAIC_reduced <- calc_waic(bundle_reduced, jags_d_reduced)
WAIC_reduced$WAIC
# 2205.511

################################################
# Check random vs. fixed:
################################################

# Check slopes (Fixed effects)

hdi.mean.reduced <- array(0, dim=c(Ncov_2010_reduced, 2))

for(i in 1:Ncov_2010_reduced){
  sub <- subset(mean.beta.df.reduced, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean.reduced[i, ] <- hdi
}
hdi.mean.reduced

# Which do not include zero?
# Lat(1), Snail_Rich(2), Amph_RA2(4), Snails_RA2(5)

# Check st.dev (Random Effects)

hdi.sd.reduced <- array(0, dim=c(Ncov_2010_reduced, 2))

for(i in 1:Ncov_2010_reduced){
  sub <- subset(sd.beta.df.reduced, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd.reduced[i, ] <- hdi
}
hdi.sd.reduced

# Snails_MDS2(3)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Try Dropping Snails_RA2  (4cov): better, but
# Try Dropping Amph_RA2    (4covB): still better, but
# Try Dropping Snail_Rich  (4covC): much better, but
# Try Dropping Lat         (4covD): BAD, so
# Try Dropping Snails_MDS2 (4covE): again, bad, so leave Snail_rich out

# Try Dropping Snail_Rich + Amph_RA2 (3cov): better than null, but 
# Try Dropping Snail_Rich + Snails_RA2 (3covB): worse, but
# Try Dropping Snail_Rich + Snails_MDS2 (3covC): worse, but
# Try Dropping Snail_Rich + Lat (3covC): much worse, so

# Try Dropping Snail_Rich + Amph_RA2 + Snails_MDS2 (2cov): no diff. from 3cov
# Try Dropping Snail_Rich + Amph_RA2 + Snails_RA2  (2covB): worse
# Try Dropping Snail_Rich + Amph_RA2 + Lat         (2covC): much worse\

# Try Dropping Snail_Rich + Amph_RA2 + Snails_MDS2 + Snails_RA2 (1cov): not as good
# Try Dropping Snail_Rich + Amph_RA2 + Snails_MDS2 + Lat (1cov)
################################################

X_2010_1covB <- data.frame(X_2010_reduced[, -c(1, 2, 3, 4)])
head(X_2010_1covB)
Ncov_2010_1covB <- ncol(X_2010_1covB)

# data:
jags_d_1covB <- list(Y=Y_2010,
                       X=X_2010_1covB,
                       Species=Species_2010,
                       Nspecies=Nspecies_2010,
                       Ncov=Ncov_2010_1covB,
                       Nobs=Nobs_2010,
                       J=J_2010)


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
                    data = jags_d_1covB, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, 
                      variable.names = params, thin=thin)
  return(out)
}

bundle_1covB <- NULL
bundle_1covB <- list(jags.parsamps[[1]][[1]],
                       jags.parsamps[[2]][[1]],
                       jags.parsamps[[3]][[1]])

class(bundle_1covB) <- "mcmc.list"

stopCluster(cl)

caterplot(bundle_1covB, parms="betas", horizontal=F)
caterplot(bundle_1covB, parms="mean.beta.post", horizontal=F)
caterplot(bundle_1covB, parms="sd.beta.post", horizontal=F)


################################################
# Calculate WAIC for step-wise removal models:
################################################
source(file="calc_waic.R")

WAIC_full <- calc_waic(bundle_full, jags_d)
WAIC_full$WAIC
# 2224.971

WAIC_reduced <- calc_waic(bundle_reduced, jags_d_reduced)
WAIC_reduced$WAIC
# 2205.511

###########
WAIC_4cov <- calc_waic(bundle_4cov, jags_d_4cov)
WAIC_4cov$WAIC
# 2203.365

WAIC_4covB <- calc_waic(bundle_4covB, jags_d_4covB)
WAIC_4covB$WAIC
# 2201.818

WAIC_4covC <- calc_waic(bundle_4covC, jags_d_4covC)
WAIC_4covC$WAIC
# 2197.258 * much better

WAIC_4covD <- calc_waic(bundle_4covD, jags_d_4covD)
WAIC_4covD$WAIC
# 2207.695  much worse

WAIC_4covE <- calc_waic(bundle_4covE, jags_d_4covE)
WAIC_4covE$WAIC
# 2203.327 worse

###########
WAIC_3cov <- calc_waic(bundle_3cov, jags_d_3cov)
WAIC_3cov$WAIC
# 2195.976* better than null

WAIC_3covB <- calc_waic(bundle_3covB, jags_d_3covB)
WAIC_3covB$WAIC
# 2198.372 worse

WAIC_3covC <- calc_waic(bundle_3covC, jags_d_3covC)
WAIC_3covC$WAIC
# 2197.466 worse

WAIC_3covD <- calc_waic(bundle_3covD, jags_d_3covD)
WAIC_3covD$WAIC
# 2201.457 worse

###########
WAIC_2cov <- calc_waic(bundle_2cov, jags_d_2cov)
WAIC_2cov$WAIC
# 2195.745 no different from 3cov

WAIC_2covB <- calc_waic(bundle_2covB, jags_d_2covB)
WAIC_2covB$WAIC
# 2198.12 worse 

WAIC_2covC <- calc_waic(bundle_2covC, jags_d_2covC)
WAIC_2covC$WAIC
# 2202.319  much worse 

###########
WAIC_1cov <- calc_waic(bundle_1cov, jags_d_1cov)
WAIC_1cov$WAIC
# 2196.572 not as good

WAIC_1covB <- calc_waic(bundle_1covB, jags_d_1covB)
WAIC_1covB$WAIC
# 2198.516 much worse

###########
WAIC_null <- calc_waic(bundle_null, jags_d_null)
WAIC_null$WAIC
# 2199.982


# BEST MODEL:
# 3 cov: Lat + Snails_MDS2 + Snails_RA2
# 2 cov: Lat + Snails_RA2

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now fit NULL model:
################################################
# data:
jags_d_null <- list(Y=Y_2010,
                    Species=Species_2010,
                    Nspecies=Nspecies_2010,
                    Nobs=Nobs_2010,
                    J=J_2010)

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
# # Calculate WAIC for NULL model:
# ################################################
source(file="calc_waic.R")

WAIC_null <- calc_waic(bundle_null, jags_d_null)
WAIC_null$WAIC
# 2199.982

################################################
# Check random vs. fixed of BEST MODEL:
################################################
mean.beta.df.best <- ggs(bundle_3cov, family="mean.beta.post")
sd.beta.df.best <- ggs(bundle_3cov, family="sd.beta.post")

# Check slopes (Fixed effects)

hdi.mean.best <- array(0, dim=c(Ncov_2010_3cov, 2))

for(i in 1:Ncov_2010_3cov){
  sub <- subset(mean.beta.df.best, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean.best[i, ] <- hdi
}
hdi.mean.best

# Which do not include zero?
# Lat

# Check st.dev (Random Effects)

hdi.sd.best <- array(0, dim=c(Ncov_2010_3cov, 2))

for(i in 1:Ncov_2010_3cov){
  sub <- subset(sd.beta.df.best, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd.best[i, ] <- hdi
}
hdi.sd.best