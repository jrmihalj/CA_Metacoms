colnames(X_all) <- colnames(X)
# Standardize each covariate: (value - mean) / 2sd
for(j in 1:15){ # all the non-factor level covariates
  X_all[, j] <- (X_all[, j] - mean(X_all[, j], na.rm=T)) / (2 * sd(X_all[, j], na.rm=T))
}

# Look for collinearity:
cor(X_all, use="complete.obs")
x11(height=10, width=15)
pairs(X_all[, c(14:25)])
pairs(X_all[,c(6:10, 14:23)])

# Collinearity:
# Lat - Long
# Long - Elev
# FOR - SSG
# Cond - TDS
# SnailsRA1 - SnailsMDS2
# SnailsRA2 - Snail Rich

# REMOVE: Long, Elev, SSG, TDS, Snails_MDS2, Snails_RA2

X_all <- X_all[, -c(2,3,7,12,19,22)]

# REMOVE: Fish
X_all <- X_all[, -c(18)]

Ncov_all <- ncol(X_all)
################################################
# JAGS model run:
################################################
library(rjags)
library(random)
load.module("lecuyer")
load.module("glm")
# data:
jags_d <- list(Y=Y_all,
               X=X_all,
               Species=Species_all,
               Nspecies=Nspecies_all,
               Ncov=Ncov_all,
               Nobs=Nobs_all,
               J=J_all)

# parameters:
params <- c("alpha", "betas", "p.detect", 
            "p.include", "sd.beta.post", "mean.beta.post",
            "z")

jinits <- function() {
  list(
    z=ifelse(Y_all > 0, 1, 0),
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
# 6338.863

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Check 'best' model:
# ################################################

# Check to see which covariate effects (slopes) != 0

source(file="HDI.R")
hdi.mean <- array(0, dim=c(Ncov_all, 2))

for(i in 1:Ncov_all){
  sub <- subset(mean.beta.df, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean[i, ] <- hdi
}
hdi.mean

# Which do not include zero?
# veg_s(6) Cond(8) Snail_Rich(11)* Amph_MDS1(12) Amph_RA2(16) Snails_RA2(17)

# Check to see which covariate effects are significantly random (st.dev > 0)

hdi.sd <- array(0, dim=c(Ncov_all, 2))

for(i in 1:Ncov_all){
  sub <- subset(sd.beta.df, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd[i, ] <- hdi
}
hdi.sd

# Which are greater than zero?
# Snail_Rich(11) Amph_RA1(15) Amph_RA2(16)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now refit reduced model:
################################################

# X_reduced for reduced covariates:

X_all_reduced <- X_all[, c(6,8,11,12,15,16,17)]
head(X_all_reduced)
Ncov_all_reduced <- ncol(X_all_reduced)

# data:
jags_d_reduced <- list(Y=Y_all,
                       X=X_all_reduced,
                       Species=Species_all,
                       Nspecies=Nspecies_all,
                       Ncov=Ncov_all_reduced,
                       Nobs=Nobs_all,
                       J=J_all)

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
# 6328.371

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
beta.df.reduced <- ggs(bundle_reduced, family="betas")
sd.beta.df.reduced <- ggs(bundle_reduced, family="sd.beta.post")
mean.beta.df.reduced <- ggs(bundle_reduced, family="mean.beta.post")


caterplot(bundle_reduced, parms="betas", horizontal=F)
caterplot(bundle_reduced, parms="mean.beta.post", horizontal=F)
caterplot(bundle_reduced, parms="sd.beta.post", horizontal=F)

# Check slopes (Fixed effects)
source(file="HDI.R")
hdi.mean.reduced <- array(0, dim=c(Ncov_all_reduced, 2))

for(i in 1:Ncov_all_reduced){
  sub <- subset(mean.beta.df.reduced, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean.reduced[i, ] <- hdi
}
hdi.mean.reduced

# Which do not include zero?
# 1,2,3,4,7

# Check st.dev (Random Effects)

hdi.sd.reduced <- array(0, dim=c(Ncov_all_reduced, 2))

for(i in 1:Ncov_all_reduced){
  sub <- subset(sd.beta.df.reduced, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd.reduced[i, ] <- hdi
}
hdi.sd.reduced

# 3, 5, 6
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Try Dropping Cond*             (6cov): best
# Try Dropping Amph_MDS1:       (6covB): 
# Try Dropping Amph_RA1         (6covC):
# Try Dropping Snail_Rich       (6covD):
# Try Dropping veg_s            (6covE):
# Try Dropping Amph_RA2*        (6covF):
# Try Dropping Snails_RA2       (6covG):


# Try Dropping Cond + Amph_MDS1 (5cov): 
# Try Dropping Cond + veg_s     (5covB):


# Try Dropping Cond + Amph_MDS1 + veg_s (4cov): no different, so
# Try Dropping Cond + Amph_MDS1 + Snails_RA2 (4covB): worse, so
# Try Dropping Cond + Amph_MDS1 + veg_s + Amph_RA2 (3cov): no different, so
# Try Dropping Cond + Amph_MDS1 + veg_s + Snails_RA2 (3covB): worse, so
# Try Dropping Cond + Amph_MDS1 + Amph_RA2 + Amph_RA1 (3covC):
# Try Dropping Cond + Amph_MDS1 + veg_s + Amph_RA2 + Amph_RA1 (2cov): bad, so
# Try Dropping Cond + Amph_MDS1 + veg_s + Amph_RA2 + Snails_RA2 (2covB):not as bad:
# Try Dropping All but Snail_Rich (1 cov): worse,
# Try Dropping All but Snails_RA2 (1 covB): worse again,
# Try Dropping All but Amph_RA1 (1 covC): best so far

################################################

X_all_5covB <- (X_all_reduced[, -c(1,2)])
head(X_all_5covB)
Ncov_all_5covB <- ncol(X_all_5covB)

# data:
jags_d_5covB <- list(Y=Y_all,
                     X=X_all_5covB,
                     Species=Species_all,
                     Nspecies=Nspecies_all,
                     Ncov=Ncov_all_5covB,
                     Nobs=Nobs_all,
                     J=J_all)


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
                    data = jags_d_5covB, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, 
                      variable.names = params, thin=thin)
  return(out)
}

bundle_5covB <- NULL
bundle_5covB <- list(jags.parsamps[[1]][[1]],
                     jags.parsamps[[2]][[1]],
                     jags.parsamps[[3]][[1]])

class(bundle_5covB) <- "mcmc.list"

stopCluster(cl)

caterplot(bundle_5covB, parms="betas", horizontal=F)
caterplot(bundle_5covB, parms="mean.beta.post", horizontal=F)
caterplot(bundle_5covB, parms="sd.beta.post", horizontal=F)
caterplot(bundle_5covB, parms="alpha", horizontal=F)
caterplot(bundle_5covB, parms="p.detect", horizontal=F)


################################################
# Calculate WAIC for step-wise removal models:
################################################
source(file="calc_waic.R")

WAIC_full <- calc_waic(bundle_full, jags_d)
WAIC_full$WAIC
# 6338.863

WAIC_reduced <- calc_waic(bundle_reduced, jags_d_reduced)
WAIC_reduced$WAIC
# 6328.371

##########
WAIC_6cov <- calc_waic(bundle_6cov, jags_d_6cov)
WAIC_6cov$WAIC
# 6325.669 **

WAIC_6covB <- calc_waic(bundle_6covB, jags_d_6covB)
WAIC_6covB$WAIC
# 6327.671 (worse)

WAIC_6covC <- calc_waic(bundle_6covC, jags_d_6covC)
WAIC_6covC$WAIC
# 6330.751 (worse)

WAIC_6covD <- calc_waic(bundle_6covD, jags_d_6covD)
WAIC_6covD$WAIC
# 6329.091 (worse)

WAIC_6covE <- calc_waic(bundle_6covE, jags_d_6covE)
WAIC_6covE$WAIC
# 6329.058 (worse)

WAIC_6covF <- calc_waic(bundle_6covF, jags_d_6covF)
WAIC_6covF$WAIC
# 6326.651

WAIC_6covG <- calc_waic(bundle_6covG, jags_d_6covG)
WAIC_6covG$WAIC
# 6331.006

##########
WAIC_5cov <- calc_waic(bundle_5cov, jags_d_5cov)
WAIC_5cov$WAIC
# 6322.347 ***

WAIC_5covB <- calc_waic(bundle_5covB, jags_d_5covB)
WAIC_5covB$WAIC
# 6323.663 

##########
WAIC_4cov <- calc_waic(bundle_4cov, jags_d_4cov)
WAIC_4cov$WAIC
# 6322.89 ***

WAIC_4covB <- calc_waic(bundle_4covB, jags_d_4covB)
WAIC_4covB$WAIC
# 6323.904 worse...

WAIC_3cov <- calc_waic(bundle_3cov, jags_d_3cov)
WAIC_3cov$WAIC
# 6322.609 ***

WAIC_3covB <- calc_waic(bundle_3covB, jags_d_3covB)
WAIC_3covB$WAIC
# 6324.36 worse...

WAIC_3covB <- calc_waic(bundle_3covB, jags_d_3covB)
WAIC_3covB$WAIC
# 6326.175 worse still...

WAIC_2cov <- calc_waic(bundle_2cov, jags_d_2cov)
WAIC_2cov$WAIC
# 6326.345 a fair bit worse

WAIC_2covB <- calc_waic(bundle_2covB, jags_d_2covB)
WAIC_2covB$WAIC
# 6323.09 Not great

WAIC_1cov <- calc_waic(bundle_1cov, jags_d_1cov)
WAIC_1cov$WAIC
# 6325.948 worse

WAIC_1covB <- calc_waic(bundle_1covB, jags_d_1covB)
WAIC_1covB$WAIC
# 6325.396 worse

WAIC_1covC <- calc_waic(bundle_1covC, jags_d_1covC)
WAIC_1covC$WAIC
# 6320.836

WAIC_null <- calc_waic(bundle_null, jags_d_null)
WAIC_null$WAIC
# 6323.25 ***

# BEST MODEL:
# 


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now fit a null model (only random intercepts):
################################################

# data:
jags_d_null <- list(Y=Y_all,
                    Species=Species_all,
                    Nspecies=Nspecies_all,
                    Nobs=Nobs_all,
                    J=J_all)

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
# 6323.25

caterplot(bundle_null, parms="alpha", horizontal=F)
caterplot(bundle_null, parms="p.detect", horizontal=F)
