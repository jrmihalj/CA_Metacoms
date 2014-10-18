colnames(X_2009) <- colnames(X)

# Standardize each covariate: (value - mean) / 2sd
for(j in 1:15){ # all the non-factor level covariates
  X_2009[, j] <- (X_2009[, j] - mean(X_2009[, j], na.rm=T)) / (2 * sd(X_2009[, j], na.rm=T))
}

# Look for collinearity:
cor(X_2009, use="complete.obs")
x11(height=10, width=13)
pairs(X_2009[,c(6:10, 14:23)])

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

caterplot(bundle_full, parms="mean.beta.post", horizontal=F)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Calculate WAIC for full model:
# ################################################
source(file="calc_waic.R")

WAIC_full <- calc_waic(bundle_full, jags_d)
WAIC_full$WAIC
# 1486.244
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Check 'best' model:
# ################################################

# Check to see which covariate effects (slopes) != 0

source(file="HDI.R")
hdi.mean <- array(0, dim=c(Ncov_2009, 2))

for(i in 1:Ncov_2009){
  sub <- subset(mean.beta.df, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean[i, ] <- hdi
}
hdi.mean

# Which do not include zero?
# Amph_MDS2 (11)

# Check to see which covariate effects are significantly random (st.dev > 0)

hdi.sd <- array(0, dim=c(Ncov_2009, 2))

for(i in 1:Ncov_2009){
  sub <- subset(sd.beta.df, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd[i, ] <- hdi
}
hdi.sd

# Which are greater than zero?
# Aspect (3), Snails_RA1 (15) Snails_RA2 (16)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now refit reduced model:
################################################

# X_reduced for reduced covariates:

X_2009_reduced <- X_2009[, c(3,11,12,15,16)]
Ncov_2009_reduced <- ncol(X_2009_reduced)

# data:
jags_d_reduced <- list(Y=Y_2009,
                    X=X_2009_reduced,
                    Species=Species_2009,
                    Nspecies=Nspecies_2009,
                    Ncov=Ncov_2009_reduced,
                    Nobs=Nobs_2009,
                    J=J_2009)

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
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Calculate WAIC for reduced model:
# ################################################
source(file="calc_waic.R")

WAIC_reduced <- calc_waic(bundle_reduced, jags_d_reduced)
WAIC_reduced$WAIC
# 1459.747 #better

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Check random vs. fixed:
################################################

# Check slopes (Fixed effects)

hdi.mean.reduced <- array(0, dim=c(Ncov_2009_reduced, 2))

for(i in 1:Ncov_2009_reduced){
  sub <- subset(mean.beta.df.reduced, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean.reduced[i, ] <- hdi
}
hdi.mean.reduced

# Which do not include zero?
# Amph_MDS2 (2)

# Check st.dev (Random Effects)

hdi.sd.reduced <- array(0, dim=c(Ncov_2009_reduced, 2))

for(i in 1:Ncov_2009_reduced){
  sub <- subset(sd.beta.df.reduced, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd.reduced[i, ] <- hdi
}
hdi.sd.reduced

# Snails_RA2 only

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Try Dropping Snails_RA1  (4cov): best
# Try Dropping Snails_RA2  (4covB):
# Try Dropping Snails_MDS1 (4covC):
# Try Dropping Amph_MDS2   (4covD):
# Try Dropping Aspect      (4covE):

# Try Dropping Snails_RA1 + Snails_MDS1 (3cov): best** 
# Try Dropping Snails_RA1 + Amph_MDS2   (3covB):

# Try Dropping Snails_RA1 + Snails_MDS1 + Aspect      (2covA):  
# Try Dropping Snails_RA1 + Snails_MDS1 + Snails_RA2  (2covB):

# Try Only Aspect     (1cov): 
# Try only Amph_MDS2  (1covB):
# Try only Snails_RA2 (1covC):
################################################

X_2009_1covC <- data.frame(X_2009_reduced[, c(5)])
head(X_2009_1covC)
Ncov_2009_1covC <- ncol(X_2009_1covC)

# data:
jags_d_1covC <- list(Y=Y_2009,
                    X=X_2009_1covC,
                    Species=Species_2009,
                    Nspecies=Nspecies_2009,
                    Ncov=Ncov_2009_1covC,
                    Nobs=Nobs_2009,
                    J=J_2009)


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

caterplot(bundle_1covC, parms="betas", horizontal=F)
caterplot(bundle_1covC, parms="mean.beta.post", horizontal=F)
caterplot(bundle_1covC, parms="sd.beta.post", horizontal=F)


################################################
# Calculate WAIC for step-wise removal models:
################################################
source(file="calc_waic.R")

WAIC_full <- calc_waic(bundle_full, jags_d)
WAIC_full$WAIC
# 1486.244

WAIC_reduced <- calc_waic(bundle_reduced, jags_d_reduced)
WAIC_reduced$WAIC
# 1459.747

#########
WAIC_4cov <- calc_waic(bundle_4cov, jags_d_4cov)
WAIC_4cov$WAIC
# 1456.439 ***

WAIC_4covB <- calc_waic(bundle_4covB, jags_d_4covB)
WAIC_4covB$WAIC
# 1461.954 bad

WAIC_4covC <- calc_waic(bundle_4covC, jags_d_4covC)
WAIC_4covC$WAIC
# 1459.04 bad

WAIC_4covD <- calc_waic(bundle_4covD, jags_d_4covD)
WAIC_4covD$WAIC
# 1459.923 bad

WAIC_4covE <- calc_waic(bundle_4covE, jags_d_4covE)
WAIC_4covE$WAIC
# 1460.345 bad

#########
WAIC_3cov <- calc_waic(bundle_3cov, jags_d_3cov)
WAIC_3cov$WAIC
# 1454.78 **** best model 

WAIC_3covB <- calc_waic(bundle_3covB, jags_d_3covB)
WAIC_3covB$WAIC
# 1456.6

#########
WAIC_2covA <- calc_waic(bundle_2covA, jags_d_2covA)
WAIC_2covA$WAIC
# 1456.359 bad

WAIC_2covB <- calc_waic(bundle_2covB, jags_d_2covB)
WAIC_2covB$WAIC
# 1458.04 bad

#########
WAIC_1cov <- calc_waic(bundle_1cov, jags_d_1cov)
WAIC_1cov$WAIC
# 1457.327 

WAIC_1covB <- calc_waic(bundle_1covB, jags_d_1covB)
WAIC_1covB$WAIC
# 1458.802 

WAIC_1covC <- calc_waic(bundle_1covC, jags_d_1covC)
WAIC_1covC$WAIC
# 1457.126 


#########
WAIC_null <- calc_waic(bundle_null, jags_d_null)
WAIC_null$WAIC
# 1458.632


# BEST MODEL:
# Apsect, Amph_MDS2, Snails_RA2 (3cov)


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now fit a null model (only random intercepts):
################################################

# data:
jags_d_null <- list(Y=Y_2009,
                       Species=Species_2009,
                       Nspecies=Nspecies_2009,
                       Nobs=Nobs_2009,
                       J=J_2009)

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
# 1458.632



################################################
# Check random vs. fixed of BEST MODEL:
################################################
mean.beta.df.best <- ggs(bundle_3cov, family="mean.beta.post")
sd.beta.df.best <- ggs(bundle_3cov, family="sd.beta.post")

# Check slopes (Fixed effects)

hdi.mean.best <- array(0, dim=c(Ncov_2009_3cov, 2))

for(i in 1:Ncov_2009_3cov){
  sub <- subset(mean.beta.df.best, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean.best[i, ] <- hdi
}
hdi.mean.best

# Which do not include zero?
# NONE: NO SIG. FIXED EFFECTS

# Check st.dev (Random Effects)

hdi.sd.best <- array(0, dim=c(Ncov_2009_3cov, 2))

for(i in 1:Ncov_2009_3cov){
  sub <- subset(sd.beta.df.best, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd.best[i, ] <- hdi
}
hdi.sd.best

# Snails_RA2 only
#(.023 - 1.90)








################################################
# Extract 'linear predictor' of the model: logit(psi)
################################################

linpred.reduced <- NULL

for(i in 1:Nobs){
  sub <- subset(psi.df.reduced, Parameter==paste("psi[", i, "]", sep=""))$value
  sub <- log(sub/(1-sub))
  linpred.reduced[i] <- mean(sub)
}

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now fit the model that has NO RANDOM TERMS (i.e. no species-level variability in beta[j]:
################################################
# (I deleted most of the data, just to save on space)

# initialize model:

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

jags.parsamps <- NULL
jags.parsamps <- foreach(i=1:3, .packages=c('rjags','random')) %dopar% {
  #setwd("C:\Users\Joe\Documents\GitHub\CA_Metacoms")
  store<-1000
  nadap<-20000
  nburn<-50000
  thin<-50
  mod <- jags.model(file = "MLM_model_NoRandom.txt", 
                    data = jags_d_best, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, 
                      variable.names = params_best, thin=thin)
  return(out)
}

bundle_norand <- NULL
bundle_norand <- list(jags.parsamps[[1]][[1]],
                      jags.parsamps[[2]][[1]],
                      jags.parsamps[[3]][[1]])

class(bundle_norand) <- "mcmc.list"

stopCluster(cl)

################################################
# Check Convergence:
################################################
library(ggmcmc)
library(mcmcplots)
alpha.df.norand <- ggs(bundle_norand, family="alpha")
beta.df.norand <- ggs(bundle_norand, family="betas")
p.detect.df.norand <- ggs(bundle_norand, family="p.detect")
psi.df.norand <- ggs(bundle_norand, family="psi")

ggs_Rhat(alpha.df.norand)
ggs_Rhat(beta.df.norand)
ggs_Rhat(p.detect.df.norand)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(bundle_norand, parms="betas", horizontal=F)
caterplot(bundle_norand, parms="alpha", horizontal=F)

################################################
# Extract 'linear predictor' of the model: logit(psi)
################################################

linpred.norand <- NULL

for(i in 1:Nobs){
  sub <- subset(psi.df.norand, Parameter==paste("psi[", i, "]", sep=""))$value
  sub <- log(sub/(1-sub))
  linpred.norand[i] <- mean(sub)
}

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Conduct PCA to ordinate sites, map environmental effects
################################################
MLM.fitted <- array(linpred.best - linpred.norand, c(Nobs/Nspecies, Nspecies))

rownames(MLM.fitted)=c(1:Nsite)
colnames(MLM.fitted)=paste0("Species", 1:Nspecies)
# standardize over spp

MLM.fitted.standard <- MLM.fitted
for(j in 1:Nspecies){
  MLM.fitted.standard[,j] <- (MLM.fitted[,j]-mean(MLM.fitted[,j]))/sd(MLM.fitted[,j])
} 

ss <- cor(MLM.fitted.standard)
U <- svd(ss)
mlm.fit <- MLM.fitted.standard %*% U$v
mlm.fit <- mlm.fit[,1:2]

# environmental variables (only those with significant random effects)
envir.vars <- Xcov[, c(1:4)]
mlm.envir <- NULL
for(j in 1:ncol(envir.vars)){
  mlm.envir <- cbind(mlm.envir, envir.vars[,j]*mlm.fit[,1],envir.vars[,j]*mlm.fit[,2])
}

envir.points <- t(array(colMeans(mlm.envir),c(2,dim(mlm.envir)[2]/2)))

# plot mlm

x11(height=5, width=5)
plot(-mlm.fit,xlab="PC1",ylab="PC2",type="n")
#text(-mlm.fit,label=c(1:(Nobs/Nspecies)),cex=.5)
points(-mlm.fit, pch=19, cex=0.5)

arrow.coordMLM <- cbind(array(0,dim(envir.points)),-envir.points)

arrows(arrow.coordMLM[,1],arrow.coordMLM[,2],arrow.coordMLM[,3],arrow.coordMLM[,4], 
       code=2, col="black", length=0.05, lwd=.8)

text(1.3*-envir.points,label=c("Cov1", "Cov2", "Cov3", "Cov4"),cex=1, font=2)

