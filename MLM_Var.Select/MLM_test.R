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

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# JAGS model run:
################################################
library(rjags)
# data:
jags_d <- list(Y=Y,
               X=X,
               Species=Species,
               Nspecies=Nspecies,
               Ncov=Ncov,
               Nobs=Nobs,
               J=J)

zinit <- NULL
zinit <- ifelse(Y > 0, 1, 0)

store<-10
nadap<-1000
nburn<-5
thin<-5

# set wd:
setwd("~/Documents/Thesis Research/CA metacoms/CAmetacomsCode/MLM_Var.Select")

# parameters:
params <- c("alpha", "betas", "I", "tau.beta")

# initialize model:
mod <- NULL
mod <- jags.model(file = "MLM_model.txt", 
                  data = jags_d, n.chains = 3, n.adapt=nadap,
                  inits = list(z=zinit))
update(mod, n.iter=nburn)

out <- NULL
out <- coda.samples(mod, n.iter = store*thin, variable.names = params, thin=thin)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

################################################
# Check Convergence:
################################################
library(ggmcmc)
library(mcmcplots)
alpha.df <- ggs(out, family="alpha")
beta.df <- ggs(out, family="betas")
I.df <- ggs(out, family="I")
tau.beta.df <- ggs(out, family="tau.beta")

ggs_Rhat(alpha.df)
ggs_Rhat(beta.df)
