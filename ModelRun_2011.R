colnames(X_2011) <- colnames(X)

# Standardize each covariate: (value - mean) / 2sd
for(j in 1:15){ # all the non-factor level covariates
  X_2011[, j] <- (X_2011[, j] - mean(X_2011[, j], na.rm=T)) / (2 * sd(X_2011[, j], na.rm=T))
}

# Look for collinearity:
cor(X_2011, use="complete.obs")
x11(height=10, width=15)
pairs(X_2011[,13:23])
pairs(X_2011[,c(6:10, 14:23)])

#Collinearity:
#Lat-long 
#Long-Elev (quadratic)
#FOR-SSG
#Cond-TDS
#OpenW-DOmg
#Snails_RA1 - Snails_MDS2
#Snails_RA1 - Snail_Rich
#Amph_RA1 - Amph_MDS2


# REMOVE: Long, Elev, SSG, TDS, OpenW, Snails_MDS2, Snails_RA1, Amph_MDS2

X_2011 <- X_2011[, -c(2,3,7,10,12,17,19,22)]

# REMOVE: Fish

X_2011 <- X_2011[, -16]

Ncov_2011 <- ncol(X_2011)
################################################
# JAGS model run:
################################################
library(rjags)
library(random)
load.module("lecuyer")
load.module("glm")
# data:
jags_d <- list(Y=Y_2011,
               X=X_2011,
               Species=Species_2011,
               Nspecies=Nspecies_2011,
               Ncov=Ncov_2011,
               Nobs=Nobs_2011,
               J=J_2011)

# parameters:
params <- c("alpha", "betas", "p.detect", 
            "p.include", "sd.beta.post", "mean.beta.post",
            "z")

jinits <- function() {
  list(
    z=ifelse(Y_2011 > 0, 1, 0),
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
# 1735.612

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Check 'best' model:
# ################################################

# Check to see which covariate effects (slopes) != 0

source(file="HDI.R")
hdi.mean <- array(0, dim=c(Ncov_2011, 2))

for(i in 1:Ncov_2011){
  sub <- subset(mean.beta.df, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean[i, ] <- hdi
}
hdi.mean

# Which do not include zero?
# Lat (1), FOR (4), Snail_Rich (10), Amph_RA2 (14), Snails_RA2 (15)

# Check to see which covariate effects are significantly random (st.dev > 0)

hdi.sd <- array(0, dim=c(Ncov_2011, 2))

for(i in 1:Ncov_2011){
  sub <- subset(sd.beta.df, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd[i, ] <- hdi
}
hdi.sd

# Which are greater than zero?
# Aspect (3), veg_s (6), Snail_Rich (10), Amph_RA1 (13), Amph_RA2 (14)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now refit reduced model:
################################################

# X_reduced for reduced covariates:

X_2011_reduced <- X_2011[, c(1,3,4,6, 10, 13, 14, 15)]
Ncov_2011_reduced <- ncol(X_2011_reduced)

# data:
jags_d_reduced <- list(Y=Y_2011,
                       X=X_2011_reduced,
                       Species=Species_2011,
                       Nspecies=Nspecies_2011,
                       Ncov=Ncov_2011_reduced,
                       Nobs=Nobs_2011,
                       J=J_2011)

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

ggs_Rhat(alpha.df.reduced)
ggs_Rhat(beta.df.reduced)
ggs_Rhat(sd.beta.df.reduced)
ggs_Rhat(p.detect.df.reduced)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(bundle_reduced, parms="betas", horizontal=F)
caterplot(bundle_reduced, parms="mean.beta.post", horizontal=F)
caterplot(bundle_reduced, parms="sd.beta.post", horizontal=F)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Calculate WAIC for reduced model:
# ################################################
source(file="calc_waic.R")

WAIC_reduced <- calc_waic(bundle_reduced, jags_d_reduced)
WAIC_reduced$WAIC
# 1656.541

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Check random vs. fixed:
################################################

# Check slopes (Fixed effects)

hdi.mean.reduced <- array(0, dim=c(Ncov_2011_reduced, 2))

for(i in 1:Ncov_2011_reduced){
  sub <- subset(mean.beta.df.reduced, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean.reduced[i, ] <- hdi
}
hdi.mean.reduced

# Which do not include zero?
# Lat (1), Snail_Rich (5), Amph_RA2 (7), Snails_RA2 (8)

# Check st.dev (Random Effects)

hdi.sd.reduced <- array(0, dim=c(Ncov_2011_reduced, 2))

for(i in 1:Ncov_2011_reduced){
  sub <- subset(sd.beta.df.reduced, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd.reduced[i, ] <- hdi
}
hdi.sd.reduced

# Which do not include zero?
# veg_s (4), Snail_Rich (5), Amph_RA1 (6), Snails_RA2 (8)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Try Dropping Aspect      (7cov):
# Try Dropping Lat         (7covB):
# Try Dropping FOR         (7covC):
# Try Dropping veg_s       (7covD): best
# Try Dropping Snail_Rich  (7covE):
# Try Dropping Amph_RA1    (7covF):
# Try Dropping Amph_RA2    (7covG): 
# Try Dropping Snails_RA2  (7covH):

# Try Dropping veg_s + Lat        (6cov): best
# Try Dropping veg_s + Amph_RA2   (6covB):
# Try Dropping veg_s + Amph_RA1   (6covC):
# Try Dropping veg_s + FOR        (6covD):
# Try Dropping veg_s + Snails_RA2 (6covE):
# Try Dropping veg_s + Snail_Rich (6covF):


# Try Dropping veg_s + Lat + Aspect     (5cov):
# Try Dropping veg_s + Lat + Snails_RA2 (5covB): best
# Try Dropping veg_s + Lat + Snail_Rich (5covC):
# Try Dropping veg_s + Lat + Amph_RA1   (5covD):

# Try Dropping veg_s + Lat + Snails_RA2 + For        (4cov): 
# Try Dropping veg_s + Lat + Snails_RA2 + Aspect     (4covB): 
# Try Dropping veg_s + Lat + Snails_RA2 + Snail_Rich (4covC):
# Try Dropping veg_s + Lat + Snails_RA2 + Amph_RA1   (4covD):
# Try Dropping veg_s + Lat + Snails_RA2 + Amph_RA2   (4covE):

# Try Dropping veg_s + Lat + Snails_RA2 + Snail_Rich + Amph_RA1 (3cov)
# Try Dropping veg_s + Lat + Snails_RA2 + Snail_Rich + FOR      (3covB): best
# Try Dropping veg_s + Lat + Snails_RA2 + Snail_Rich + Amph_RA2 (3covC)
# Try Dropping veg_s + Lat + Snails_RA2 + Snail_Rich + Aspect   (3covD)

# Try Dropping veg_s + Lat + Snails_RA2 + Snail_Rich + FOR + Amph_RA1 (2cov): best
# Try Dropping veg_s + Lat + Snails_RA2 + Snail_Rich + FOR + Aspect   (2covB):
# Try Dropping veg_s + Lat + Snails_RA2 + Snail_Rich + FOR + Amph_RA2 (2covC);

# Try with just Aspect    (1cov)
# Try with just Amph_RA2  (1covB):
# Just curious: Try with just Amph_RA1:
################################################

X_2011_1covC <- data.frame(X_2011_reduced[, c(6)])
head(X_2011_1covC)
Ncov_2011_1covC <- ncol(X_2011_1covC)

# data:
jags_d_1covC <- list(Y=Y_2011,
                    X=X_2011_1covC,
                    Species=Species_2011,
                    Nspecies=Nspecies_2011,
                    Ncov=Ncov_2011_1covC,
                    Nobs=Nobs_2011,
                    J=J_2011)


# initialize model:

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

jags.parsamps <- NULL
jags.parsamps <- foreach(i=1:3, .packages=c('rjags','random')) %dopar% {
  #setwd("C:\Users\Joe\Documents\GitHub\CA_Metacoms")
  store<-1500
  nadap<-90000
  nburn<-80000
  thin<-50
  mod <- jags.model(file = "MLM_model.txt", 
                    data = jags_d_1covC, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, 
                      variable.names = params, thin=thin)
  return(out)
}

bundle_1covC <- NULL
bundle_1covC <- list(jags.parsamps[[1]][[1]],
                    jags.parsamps[[2]][[1]],
                    jags.parsamps[[3]][[1]])

class(bundle_1covC) <- "mcmc.list"

stopCluster(cl)

caterplot(bundle_2covB, parms="betas", horizontal=F)
caterplot(bundle_2covB, parms="mean.beta.post", horizontal=F)
caterplot(bundle_2covB, parms="sd.beta.post", horizontal=F)


################################################
# Calculate WAIC for step-wise removal models:
################################################
source(file="calc_waic.R")

WAIC_full <- calc_waic(bundle_full, jags_d)
WAIC_full$WAIC
# 1735.612

WAIC_reduced <- calc_waic(bundle_reduced, jags_d_reduced)
WAIC_reduced$WAIC
# 1656.541

#########
WAIC_7cov <- calc_waic(bundle_7cov, jags_d_7cov)
WAIC_7cov$WAIC
# 1650.577

WAIC_7covB <- calc_waic(bundle_7covB, jags_d_7covB)
WAIC_7covB$WAIC
# 1649.393

WAIC_7covC <- calc_waic(bundle_7covC, jags_d_7covC)
WAIC_7covC$WAIC
# 1654.355

WAIC_7covD <- calc_waic(bundle_7covD, jags_d_7covD)
WAIC_7covD$WAIC
# 1649.254

WAIC_7covE <- calc_waic(bundle_7covE, jags_d_7covE)
WAIC_7covE$WAIC
# 1650.735

WAIC_7covF <- calc_waic(bundle_7covF, jags_d_7covF)
WAIC_7covF$WAIC
# 1649.476

WAIC_7covG <- calc_waic(bundle_7covG, jags_d_7covG)
WAIC_7covG$WAIC
# 1654.774

WAIC_7covH <- calc_waic(bundle_7covH, jags_d_7covH)
WAIC_7covH$WAIC
# 1650.011

#########
WAIC_6cov <- calc_waic(bundle_6cov, jags_d_6cov)
WAIC_6cov$WAIC
# 1642.265 ****

WAIC_6covB <- calc_waic(bundle_6covB, jags_d_6covB)
WAIC_6covB$WAIC
# 1650.567

WAIC_6covC <- calc_waic(bundle_6covC, jags_d_6covC)
WAIC_6covC$WAIC
# 1646.675

WAIC_6covD <- calc_waic(bundle_6covD, jags_d_6covD)
WAIC_6covD$WAIC
# 1647.946

WAIC_6covE <- calc_waic(bundle_6covE, jags_d_6covE)
WAIC_6covE$WAIC
# 1644.341

WAIC_6covF <- calc_waic(bundle_6covF, jags_d_6covF)
WAIC_6covF$WAIC
# 1645.704

#########
WAIC_5cov <- calc_waic(bundle_5cov, jags_d_5cov)
WAIC_5cov$WAIC
# 1643.975 

WAIC_5covB <- calc_waic(bundle_5covB, jags_d_5covB)
WAIC_5covB$WAIC
# 1639.356 ***

WAIC_5covC <- calc_waic(bundle_5covC, jags_d_5covC)
WAIC_5covC$WAIC
# 1641.868

WAIC_5covD <- calc_waic(bundle_5covD, jags_d_5covD)
WAIC_5covD$WAIC
# 1644.152

#########
WAIC_4cov <- calc_waic(bundle_4cov, jags_d_4cov)
WAIC_4cov$WAIC
# 1640.556

WAIC_4covB <- calc_waic(bundle_4covB, jags_d_4covB)
WAIC_4covB$WAIC
# 1642.441

WAIC_4covC <- calc_waic(bundle_4covC, jags_d_4covC)
WAIC_4covC$WAIC
# 1638.396 ***

WAIC_4covD <- calc_waic(bundle_4covD, jags_d_4covD)
WAIC_4covD$WAIC
# 1638.999

WAIC_4covE <- calc_waic(bundle_4covE, jags_d_4covE)
WAIC_4covE$WAIC
# 1644.547

#########
WAIC_3cov <- calc_waic(bundle_3cov, jags_d_3cov)
WAIC_3cov$WAIC
# 1637.76 

WAIC_3covB <- calc_waic(bundle_3covB, jags_d_3covB)
WAIC_3covB$WAIC
# 1637.414 **

WAIC_3covC <- calc_waic(bundle_3covC, jags_d_3covC)
WAIC_3covC$WAIC
# 1640.614

WAIC_3covD <- calc_waic(bundle_3covD, jags_d_3covD)
WAIC_3covD$WAIC
# 1638.132

#########
WAIC_2cov <- calc_waic(bundle_2cov, jags_d_2cov)
WAIC_2cov$WAIC
# 1636.657 **

WAIC_2covB <- calc_waic(bundle_2covB, jags_d_2covB)
WAIC_2covB$WAIC
# 1636.741 **

WAIC_2covC <- calc_waic(bundle_2covC, jags_d_2covC)
WAIC_2covC$WAIC
# 1638.447

#########
WAIC_1cov <- calc_waic(bundle_1cov, jags_d_1cov)
WAIC_1cov$WAIC
# 1636.951

WAIC_1covB <- calc_waic(bundle_1covB, jags_d_1covB)
WAIC_1covB$WAIC
# 1636.445

WAIC_1covC <- calc_waic(bundle_1covC, jags_d_1covC)
WAIC_1covC$WAIC
# 1639.213

#########
WAIC_null <- calc_waic(bundle_null, jags_d_null)
WAIC_null$WAIC
# 1638.877


# BEST MODEL:
# 4 equivalent models:
# - Aspect + Amph_RA2 (2cov)
# - Amph_RA1 + Amph_RA2 (2covB) **
# - Aspect (1cov)
# - Amph_RA2 (1covB)


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now fit a null model (only random intercepts):
################################################

# data:
jags_d_null <- list(Y=Y_2011,
                    Species=Species_2011,
                    Nspecies=Nspecies_2011,
                    Nobs=Nobs_2011,
                    J=J_2011)

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
# 1638.877

# $lppd
# [1] -771.1309
# 
# $p_WAIC
# [1] 48.2028

################################################
# Check random vs. fixed of BEST MODEL:
################################################
mean.beta.df.best <- ggs(bundle_2covB, family="mean.beta.post")
sd.beta.df.best <- ggs(bundle_2covB, family="sd.beta.post")

# Check slopes (Fixed effects)

hdi.mean.best <- array(0, dim=c(Ncov_2011_2covB, 2))

for(i in 1:Ncov_2011_2covB){
  sub <- subset(mean.beta.df.best, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean.best[i, ] <- hdi
}
hdi.mean.best

# Which do not include zero?
# NONE: No fixed effects

# Check st.dev (Random Effects)

hdi.sd.best <- array(0, dim=c(Ncov_2011_2covB, 2))

for(i in 1:Ncov_2011_2covB){
  sub <- subset(sd.beta.df.best, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd.best[i, ] <- hdi
}
hdi.sd.best

