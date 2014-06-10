# Testing MLM model with correction for detection probability:
# Data from: Jackson et al. 2012

################################################
# load data
################################################
# set wd
setwd("~/Documents/Thesis Research/CA metacoms/CAmetacomsCode/MLM_Var.Select")
# use file "MLM_dataset.csv"
h <- read.csv("MLM_dataset.csv")
# remove sites with no species
h <- h[-3,]

#Normalize predictor variables:
for(i in 5:25){h[,i] <- (h[,i] - mean(h[,i]))/sd(h[,i])} 

hh <- cbind(h[,1:26],factor('ACRA')) 
names(hh)[26] <- 'PRESENCE' 

names(hh)[27] <- 'SPP' 
levels(hh$SPP) <- c("ACRA","ARTR","CATH","CHMA","DILA","GALA","GEMA","LISU","POBI","SACA","THDI","TRGR","UVGR","VIRO")

T <- dim(h)[1] 
T #[1] 54
t <- T 
t  #[1] 54
names(hh)

for(i in 27:39){
  hh[(t+1):(t+T),] <- cbind(h[,1:25], h[,i], names(h)[i])
  t <- t+T
}
t #756

hh$Elevation2 <- hh$Elevation^2
hh$LITU2 <- hh$LITU^2
hh$Ca2 <- hh$Ca^2
hh$P2 <- hh$P^2
hh$Herb2=hh$Herb^2

################################################
# end load data
################################################

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Format data for model
################################################
# All data for Bayesian model:
X <- hh[, c(5:25, 28:32)] #include all covariates
Y <- hh$PRESENCE
Species <- as.factor(as.numeric(hh$SPP))
Nspecies <- length(levels(Species))
Ncov <- ncol(X)
Nobs <- length(Y)
J <- rep(5, times=Nobs)
  
# Convert Y (pres/abs) to binomial based on simulated species detection probs.
ps <- runif(Nspecies, 0.25, 1)
pdetect <- rep(ps, each=54) # each species has 54 sites
for(i in 1:Nobs){
  if(Y[i] ==1){
    Y[i] <- rbinom(1, 5, pdetect[i])
  }
}

# Some of the X variables are highly correlated:
# Keep CA, Elevation
# Remove: OM, N, Mg, pH, K, Litu2, Ca2

X <- X[, -c(15,20,17,14,19,23)]
Ncov <- ncol(X)
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# JAGS model run:
################################################
library(rjags)
library(random)
load.module("lecuyer")
load.module("glm")
# data:
jags_d <- list(Y=Y,
               X=X,
               Species=Species,
               Nspecies=Nspecies,
               Ncov=Ncov,
               Nobs=Nobs,
               J=J)

# zinit <- NULL
# zinit <- ifelse(Y > 0, 1, 0)


# set wd:
# setwd("~/Documents/Thesis Research/CA metacoms/CAmetacomsCode/MLM_Var.Select")

# parameters:
params <- c("alpha", "betas", "I", "tau.beta", "p.detect", "p.include")

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
  nadap<-50000
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

# # initialize model:
# mod <- NULL
# mod <- jags.model(file = "MLM_model.txt", 
#                   data = jags_d, n.chains = 3, n.adapt=nadap,
#                   inits = list(z=zinit))
# update(mod, n.iter=nburn)
# 
# out <- NULL
# out <- coda.samples(mod, n.iter = store*thin, variable.names = params, thin=thin)

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
tau.beta.df <- ggs(bundle, family="tau.beta")
p.detect.df <- ggs(bundle, family="p.detect")
p.include.df <- ggs(bundle, family="p.include")

ggs_Rhat(alpha.df)
ggs_Rhat(beta.df)
ggs_Rhat(tau.beta.df)
ggs_Rhat(p.detect.df)
ggs_Rhat(p.include.df)


x11(height=4, width=11)
caterplot(bundle, parms="betas", horizontal=F, random=50)

caterplot(bundle, parms="p.detect", horizontal=F)#, val.lim=c(-1, 10))
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

ordered.mod.probs[1:10]
ordered.mods[1:5, ]

################################################
# Check random vs. fixed:
################################################

median.tau <- NULL
for(i in 1:Ncov){
  sub <- subset(tau.beta.df, Parameter==paste("tau.beta[", i, "]", sep=""))$value
  median.tau[i] <- median(sub)
}
median.tau
