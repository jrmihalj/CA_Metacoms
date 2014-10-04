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
params <- c("alpha", "betas", "I", "p.detect", 
            "p.include", "sd.beta.post", "mean.beta.post")

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
  store<-1000
  nadap<-50000
  nburn<-80000
  thin<-50
  mod <- jags.model(file = "MLM_model.txt", 
                    data = jags_d, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, variable.names = params, thin=thin)
  return(out)
}

bundle <- NULL
bundle <- list(jags.parsamps[[1]][[1]],
               jags.parsamps[[2]][[1]],
               jags.parsamps[[3]][[1]])

class(bundle) <- "mcmc.list"

stopCluster(cl)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Check Convergence:
################################################
library(ggmcmc)
library(mcmcplots)
alpha.df <- ggs(bundle, family="alpha")
beta.df <- ggs(bundle, family="betas")
#I.df <- ggs(bundle, family="I")
sd.beta.df <- ggs(bundle, family="sd.beta.post")
mean.beta.df <- ggs(bundle, family="mean.beta.post")
p.detect.df <- ggs(bundle, family="p.detect")
#p.include.df <- ggs(bundle, family="p.include")

ggs_Rhat(alpha.df)
ggs_Rhat(beta.df)
ggs_Rhat(sd.beta.df)
ggs_Rhat(p.detect.df)
ggs_Rhat(mean.beta.df)
#ggs_Rhat(p.include.df)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(bundle, parms="betas", horizontal=F)#, random=50)

caterplot(bundle, parms="mean.beta.post", horizontal=F)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ################################################
# # Check best model:
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
# 1-Lat, 9-Amph_Rich, 10-Snail_Rich, 11-Amph_MDS1, 13-Amph_RA1
# 14-Amph_RA2, 15-Snails_RA2, 16-Fish

# Check to see which covariate effects are significantly random (st.dev > 0)

hdi.sd <- array(0, dim=c(Ncov_2011, 2))

for(i in 1:Ncov_2011){
  sub <- subset(sd.beta.df, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd[i, ] <- hdi
}
hdi.sd

# Which are greater than zero?
# 1,13-16, also: 
# 3-Aspect, 5-area, 6-veg_s, 8-DOmg


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now refit best model:
################################################

# X_best for best covariates:

X_2011_best <- X_2011[, c(1,3,5,6,8,9:11,13:16)]
Ncov_2011_best <- ncol(X_2011_best)

# data:
jags_d_best <- list(Y=Y_2011,
                    X=X_2011_best,
                    Species=Species_2011,
                    Nspecies=Nspecies_2011,
                    Ncov=Ncov_2011_best,
                    Nobs=Nobs_2011,
                    J=J_2011)

# parameters:
params_best <- c("alpha", "betas", "p.detect", 
                 "sd.beta.post", "mean.beta.post", "psi", 
                 "z")

# initialize model:

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

jags.parsamps <- NULL
jags.parsamps <- foreach(i=1:3, .packages=c('rjags','random')) %dopar% {
  #setwd("C:\Users\Joe\Documents\GitHub\CA_Metacoms")
  store<-1000
  nadap<-50000
  nburn<-50000
  thin<-50
  mod <- jags.model(file = "MLM_model.txt", 
                    data = jags_d_best, n.chains = 1, n.adapt=nadap,
                    inits = jinits)
  update(mod, n.iter=nburn)
  out <- coda.samples(mod, n.iter = store*thin, 
                      variable.names = params_best, thin=thin)
  return(out)
}

bundle_best <- NULL
bundle_best <- list(jags.parsamps[[1]][[1]],
                    jags.parsamps[[2]][[1]],
                    jags.parsamps[[3]][[1]])

class(bundle_best) <- "mcmc.list"

stopCluster(cl)

################################################
# Check Convergence:
################################################
library(ggmcmc)
library(mcmcplots)
alpha.df.best <- ggs(bundle_best, family="alpha")
beta.df.best <- ggs(bundle_best, family="betas")
sd.beta.df.best <- ggs(bundle_best, family="sd.beta.post")
mean.beta.df.best <- ggs(bundle_best, family="mean.beta.post")
p.detect.df.best <- ggs(bundle_best, family="p.detect")
psi.df.best <- ggs(bundle_best, family="psi")
z.df.best <- ggs(bundle_best, family="z")

ggs_Rhat(alpha.df.best)
ggs_Rhat(beta.df.best)
ggs_Rhat(sd.beta.df.best)
ggs_Rhat(p.detect.df.best)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(bundle_best, parms="betas", horizontal=F)
caterplot(bundle_best, parms="mean.beta.post", horizontal=F)
caterplot(bundle_best, parms="sd.beta.post", horizontal=F)

################################################
# Check random vs. fixed:
################################################

# Check slopes (Fixed effects)

hdi.mean.best <- array(0, dim=c(Ncov_2011_best, 2))

for(i in 1:Ncov_2011_best){
  sub <- subset(mean.beta.df.best, Parameter==paste("mean.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.mean.best[i, ] <- hdi
}
hdi.mean.best

# Which do not include zero?
# 4-veg_s, 6-Amph_Rich, 7-Snail_Rich, 8-Amph_MDS1, 
# 10-Amph_RA2, 11-Snails_RA2, 12-Fish

# Check st.dev (Random Effects)

hdi.sd.best <- array(0, dim=c(Ncov_2011_best, 2))

for(i in 1:Ncov_2011_best){
  sub <- subset(sd.beta.df.best, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd.best[i, ] <- hdi
}
hdi.sd.best

#1,3,4,9-12

################################################
# Extract 'linear predictor' of the model: logit(psi)
################################################

linpred.best <- NULL

for(i in 1:Nobs){
  sub <- subset(psi.df.best, Parameter==paste("psi[", i, "]", sep=""))$value
  sub <- log(sub/(1-sub))
  linpred.best[i] <- mean(sub)
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

