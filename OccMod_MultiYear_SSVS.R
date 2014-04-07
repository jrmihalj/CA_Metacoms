model {
  #-------------------------------------------------------------------------------
  #######################################
  ############ ASSIGN PRIORS ############
  #######################################
  #-------------------------------------------------------------------------------
  ## Metacommunity-wide average of parameters: psi, p, phi, gamma:
  
  psiMean ~ dbeta(1,1) # Initial occurrence probability
  pMean ~ dbeta(1,1) # Detection probability
  phiMean ~ dbeta(1,1) # Persistence probability
  gamMean ~ dbeta(1,1) # Colonization probability
  
  ## Manual logit transformation of mean parameter values:
  
  lpsiMean <- log(psiMean) - log(1-psiMean)
  lpMean <- log(pMean) - log(1-pMean) 
  lphiMean <- log(phiMean) - log(1-phiMean) 
  lgamMean <- log(gamMean) - log(1-gamMean) 
  
  #-------------------------------------------------------------------------------
  ## Precision estimates for each metacommunity average:
  
  lpsiSD ~ dunif(0,100)
  lpSD ~ dunif(0,100)
  lphiSD ~ dunif(0,100)
  lgamSD ~ dunif(0,100)
  
  lpsiPrec <- pow(lpsiSD,-2)
  lpPrec <- pow(lpSD,-2)
  lphiPrec <- pow(lphiSD,-2)
  lgamPrec <- pow(lgamSD,-2)
  
  #-------------------------------------------------------------------------------
  ## Estimates for covariate effects:
  ## Stochastic search variable selection priors:
  
  bSD ~ dunif(0,100)
  bPrec_in <- pow(bSD,-2)
  bPrec[1] <- bPrec_in         #coeff. effectively zero
  bPrec[2] <- bPrec_in/1000    #non-zero coeff.
  
  
  cSD ~ dunif(0,100)  
  cPrec_in <- pow(cSD,-2)
  cPrec[1] <- cPrec_in         #coeff. effectively zero
  cPrec[2] <- cPrec_in/1000    #non-zero coeff.
  
  dSD ~ dunif(0,100)
  dPrec_in <- pow(dSD,-2)
  dPrec[1] <- dPrec_in         #coeff. effectively zero
  dPrec[2] <- dPrec_in/1000    #non-zero coeff.
  
  p_ind[1] <- 1/2
  p_ind[2] <- 1 - p_ind[1]
  
  
  #-------------------------------------------------------------------------------
  #######################################
  ########## Likelihood Model ###########
  #######################################
  #-------------------------------------------------------------------------------
  
  for (i in 1:n) { # n=number of species
    
    ## Initial occupancy state at t=1
    
    b0[i] ~ dnorm(lpsiMean, lpsiPrec)T(-12,12) # Species-specific baseline of occurrence prob.
    
    ## Species-specific covariate effects on initial occurrence:
    for (cov in 1:ncovs) {
      b_indA[i, cov] ~ dcat(p_ind[]) # returns 1 or 2
      b_incl[i, cov] <- b_indA[i, cov] - 1   # returns 0 or 1
      b[i,cov] ~ dnorm(0, bPrec[b_indA[i, cov]])
      
      for(k in 1:K){
        temp.b[i, k, cov, 1] <- b[i, cov] * x[k, cov, 1]
      }
    }
    
    lp[i,1] ~ dnorm(lpMean, lpPrec)T(-12,12) # Species-specific detection prob at t=1
    p[i,1] <- 1/(1+exp(-lp[i,1])) # Manual anti-logit
    
    for (k in 1:K) {
      lpsi[i,k,1] <- b0[i] + sum(temp.b[i, k, ,1]) # Year-specific covariate values
      psi[i,k,1] <- 1/(1 + exp(-lpsi[i,k,1]))
      mu.z[i,k,1] <- psi[i,k,1]
      z[i,k,1] ~ dbern(mu.z[i,k,1])
      mu.y[i,k,1] <- p[i,1]*z[i,k,1]
      y[i,k,1] ~ dbin(mu.y[i,k,1], J[k,1]) # J=number of hosts sampled from each site
    }
    
    ## Model of changes in occupancy state for t=2, ..., T
    
    # First, draw covariate effects:
    # - These are species-specific
    # - I assume covariate effects don't change across years
    
    for (cov in 1:ncovs) {
      c_indA[i, cov] ~ dcat(p_ind[]) # returns 1 or 2
      c_incl[i, cov] <- c_indA[i, cov] - 1   # returns 0 or 1
      c[i,cov] ~ dnorm(0, cPrec[c_indA[i, cov]])
      
      d_indA[i, cov] ~ dcat(p_ind[]) # returns 1 or 2
      d_incl[i, cov] <- d_indA[i, cov] - 1   # returns 0 or 1
      d[i,cov] ~ dnorm(0, dPrec[d_indA[i, cov]])
      
      for(t in 1:(T-1)){
        for(k in 1:K){
          temp.c[i, t, k, cov] <- c[i, cov] * x[k, cov, t]
          temp.d[i, t, k, cov] <- d[i, cov] * x[k, cov, t]
        }
      }
    }
    
    for (t in 1:(T-1)) {
      lp[i,t+1] ~ dnorm(lpMean, lpPrec)T(-12,12) 
      p[i,t+1] <- 1/(1+exp(-lp[i,t+1])) # Species- and year-specific detection prob.
      
      c0[i,t] ~ dnorm(lgamMean, lgamPrec)T(-12,12) # Species- and year-specific baseline 
      # colonization probability 
      d0[i,t] ~ dnorm(lphiMean, lphiPrec)T(-12,12) # Species- and year-specific baseline 
      # persistence probability
      
      for (k in 1:K) { # K=number of sites
        lgam[i,k,t] <- c0[i,t] + sum(temp.c[i, t, k, ])
        gam[i,k,t] <-  1/(1+exp(-lgam[i,k,t]))
        lphi[i,k,t] <- d0[i,t] + sum(temp.d[i, t, k, ])
        phi[i,k,t] <-  1/(1+exp(-lphi[i,k,t]))
        
        psi[i,k,t+1] <- phi[i,k,t] * psi[i,k,t] + gam[i,k,t] * (1 - psi[i,k,t])
        mu.z[i,k,t+1] <- (phi[i,k,t] * z[i,k,t] + gam[i,k,t] * (1-z[i,k,t]))
        z[i,k,t+1] ~ dbern(mu.z[i,k,t+1])
        mu.y[i,k,t+1] <- p[i,t+1] * z[i,k,t+1]
        y[i,k,t+1] ~ dbin(mu.y[i,k,t+1], J[k,t+1])
      } 
    }
  } 
  
} # Model closed