model {
  #-------------------------------------------------------------------------------
  #######################################
  ############ ASSIGN PRIORS ############
  #######################################
  #-------------------------------------------------------------------------------
  ## Metacommunity-wide average of parameters: psi, p, phi, gamma:
  
  psiMean ~ dbeta(1,1) # Initial occurrence probability
  pMean ~ dbeta(1,1) # Detection probability
  
  ## Manual logit transformation of mean parameter values:
  
  lpsiMean <- log(psiMean) - log(1-psiMean)
  lpMean <- log(pMean) - log(1-pMean) 
  
  #-------------------------------------------------------------------------------
  ## Precision estimates for each metacommunity average:
  
  lpsiSD ~ dunif(0,10)
  lpSD ~ dunif(0,10)
  
  lpsiPrec <- pow(lpsiSD,-2)
  lpPrec <- pow(lpSD,-2)
  
  #-------------------------------------------------------------------------------
  ## Estimates for covariate effects:
  ## Stochastic search variable selection priors:
  
  for(t in 1:4){ # Different covariate effects per year:
    for(i in 1:n){
      bSD[i,t] ~ dunif(0,50)
      bPrec_in[i,t] <- pow(bSD[i,t],-2)
      bPrec[i,t,1] <- bPrec_in[i,t]         #coeff. effectively zero
      bPrec[i,t,2] <- bPrec_in[i,t]/1000    #non-zero coeff.
    }
  }
  
  p_ind[1] <- 1/2
  p_ind[2] <- 1 - p_ind[1]
  
  ## Missing data priors:
  #-------------------------------------------------------------------------------
  for(cov in 1:11){
    sd.x[cov] ~ dunif(0,10)
    Prec.x[cov] <- pow(sd.x[cov],-2)
  }
  for(cov in 12:22){
    bern.p[cov] ~ dbeta(1,1)
  }
  
  p.hydro[1] <- 1/3; p.hydro[2] <- 1/3; p.hydro[3] <- 1/3; 
  #-------------------------------------------------------------------------------
  #######################################
  ########## Likelihood Model ###########
  #######################################
  #-------------------------------------------------------------------------------
  
  # Deal with the missing covariates:
  for(cov in 1:11){
    for(k in 1:K){
      x[k, cov] ~ dnorm(0, Prec.x[cov])
    }
  }
  for(cov in 12:22){
    for(k in 1:K){
      x[k, cov] ~ dbern(bern.p[cov])
    }
  }
  for(k in 1:K){
    x[k, 23] ~ dcat(p.hydro[])
  }
  
  # Occupancy states:
  for (i in 1:n) { # n=number of species
    
    ## Species-specific covariate effects on occurrence:
    for(t in 1:4){
      for (cov in 1:ncovs) {
        b_indA[i, cov, t] ~ dcat(p_ind[]) # returns 1 or 2
        b_incl[i, cov, t] <- b_indA[i, cov, t] - 1   # returns 0 or 1
        b[i,cov,t] ~ dnorm(0, bPrec[i, t, b_indA[i, cov, t]])T(-12,12)
      }
    }
    
    for(cov in 1:ncovs){
      for(k in 1:83){
        temp.b[i, k, cov] <- b[i, cov, 1] * x[k, cov] # speed up matrix calculations
      }
      for(k in 84:188){
        temp.b[i, k, cov] <- b[i, cov, 2] * x[k, cov] # speed up matrix calculations
      }
      for(k in 189:255){
        temp.b[i, k, cov] <- b[i, cov, 3] * x[k, cov] # speed up matrix calculations
      }
      for(k in 256:289){
        temp.b[i, k, cov] <- b[i, cov, 4] * x[k, cov] # speed up matrix calculations
      }
    }

    
    for(t in 1:4){
      b0[i,t] ~ dnorm(lpsiMean, lpsiPrec)T(-12,12) # Species-specific baseline occurrence prob. 
    }
    lp[i] ~ dnorm(lpMean, lpPrec)T(-12,12) # Species-specific detection prob at t=1
    p[i] <- 1/(1+exp(-lp[i])) # Manual anti-logit
    
    for (k in 1:83) {
      lpsi[i,k] <- b0[i,1] + sum(temp.b[i, k, ])
      psi[i,k] <- 1/(1 + exp(-lpsi[i,k]))
      z[i,k] ~ dbern(psi[i,k])
      y[i,k] ~ dbin(p[i]*z[i,k], J[k]) # J=number of hosts sampled from each site
    }
    for (k in 84:188) {
      lpsi[i,k] <- b0[i,2] + sum(temp.b[i, k, ])
      psi[i,k] <- 1/(1 + exp(-lpsi[i,k]))
      z[i,k] ~ dbern(psi[i,k])
      y[i,k] ~ dbin(p[i]*z[i,k], J[k]) # J=number of hosts sampled from each site
    }
    for (k in 189:255) {
      lpsi[i,k] <- b0[i,3] + sum(temp.b[i, k, ])
      psi[i,k] <- 1/(1 + exp(-lpsi[i,k]))
      z[i,k] ~ dbern(psi[i,k])
      y[i,k] ~ dbin(p[i]*z[i,k], J[k]) # J=number of hosts sampled from each site
    }
    for (k in 256:289) {
      lpsi[i,k] <- b0[i,4] + sum(temp.b[i, k, ])
      psi[i,k] <- 1/(1 + exp(-lpsi[i,k]))
      z[i,k] ~ dbern(psi[i,k])
      y[i,k] ~ dbin(p[i]*z[i,k], J[k]) # J=number of hosts sampled from each site
    }
  
  }  
# Close Model:
}