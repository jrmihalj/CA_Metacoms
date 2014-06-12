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
               Ncov=Ncov,
               Nobs=Nobs_2009,
               J=J_2009)

# parameters:
params <- c("alpha", "betas", "I", "p.detect", 
            "p.include", "sd.beta.post")

jinits <- function() {
  list(
    z=ifelse(Y > 0, 1, 0),
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
  setwd("~/GitHub/MLM_EcologyInSilico")
  store<-1000
  nadap<-20000
  nburn<-50000
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
I.df <- ggs(bundle, family="I")
sd.beta.df <- ggs(bundle, family="sd.beta.post")
p.detect.df <- ggs(bundle, family="p.detect")
p.include.df <- ggs(bundle, family="p.include")

ggs_Rhat(alpha.df)
ggs_Rhat(beta.df)
ggs_Rhat(sd.beta.df)
ggs_Rhat(p.detect.df)
ggs_Rhat(p.include.df)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(bundle, parms="betas", horizontal=F, random=50)

caterplot(bundle, parms="sd.beta.post", horizontal=F)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Check best model:
################################################

Ipost <- NULL
Ipost <- array(dim=c(store*3, Ncov)) # indicator variable array
for (i in 1:Ncov){
  string <- paste("I[", i, "]", sep="")
  Ipost[, i] <- subset(I.df, Parameter==string)$value
}

# what are the unique models that have nonzero posterior probability?
uniquemods <- unique(Ipost, MARGIN=1)
# how many do we have?
nmods <- dim(uniquemods)[1]
nmods
model.probabilities <- rep(NA, nmods)
for (i in 1:nmods){
  TFs <- apply(Ipost, 1, function(x) all(x == uniquemods[i,]))
  model.probabilities[i] <- sum(TFs) / (store*3)
}

sum(model.probabilities)
ordered.mod.probs <- model.probabilities[order(-model.probabilities)]
ordered.mods <- uniquemods[order(-model.probabilities), ]

ordered.mod.probs
ordered.mods

# Best model probability: 0.986, next 0.0057
# Includes variables: 1-5

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Now refit best model:
################################################
# (I deleted most of the data, just to save on space)

# X_best for best covariates:

X_best <- X[, c(1:5)]
Ncov <- ncol(X_best)

# data:
jags_d_best <- list(Y=Y,
                    X=X_best,
                    Species=Species,
                    Nspecies=Nspecies,
                    Ncov=Ncov,
                    Nobs=Nobs,
                    J=J)

# parameters:
params_best <- c("alpha", "betas", "p.detect", 
                 "sd.beta.post", "psi")
jinits <- function() {
  list(
    z=ifelse(Y > 0, 1, 0),
    .RNG.name=c("base::Super-Duper"),
    .RNG.seed=as.numeric(randomNumbers(n = 1, min = 1, max = 1e+06, col=1))
  )
}
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
  mod <- jags.model(file = "MLM_model_Best.txt", 
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
p.detect.df.best <- ggs(bundle_best, family="p.detect")
psi.df.best <- ggs(bundle_best, family="psi")

ggs_Rhat(alpha.df.best)
ggs_Rhat(beta.df.best)
ggs_Rhat(sd.beta.df.best)
ggs_Rhat(p.detect.df.best)

quartz(height=4, width=11)
x11(height=4, width=11)
caterplot(bundle_best, parms="betas", horizontal=F)
caterpoints(as.vector(Beta)[1:(Nspecies*Ncov)], horizontal=F)
caterplot(bundle_best, parms="sd.beta.post", horizontal=F)

################################################
# Check random vs. fixed:
################################################
# If 95 HDI of st.dev. of beta[j] overlaps zero, then fixed effect.
# (i.e. no significant variability in effect among species)

source(file="HDI.R")
hdi.sd <- array(0, dim=c(Ncov, 2))

for(i in 1:Ncov){
  sub <- subset(sd.beta.df.best, Parameter==paste("sd.beta.post[",i,"]",sep=""))$value
  hdi <- HDI(sub) #HDI of st.dev. for each covariate
  
  hdi.sd[i, ] <- hdi
}
hdi.sd
# The model estimated that only covariate 5 has a fixed effect


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
