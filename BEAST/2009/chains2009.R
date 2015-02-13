ch1 <- NULL; ch2 <- NULL; ch3 <- NULL;
ch1 <- read.coda(output.file = , 
                 index.file = , quiet=T)
ch2 <- read.coda(output.file = paste(uniqueID, "chain2.txt", sep=""), 
                 index.file = paste(uniqueID, "index.txt", sep=""), quiet=T)
ch3 <- read.coda(output.file = paste(uniqueID, "chain3.txt", sep=""), 
                 index.file = paste(uniqueID, "index.txt", sep=""), quiet=T)

out <- NULL
out <- mcmc.list(list(ch1, ch2, ch3))


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