# Store Z values from the best model's posterior
# Assemble into Z_post matrices
# Calculate metacommunity structure metrics for each Z_post
# Assign metacommunity structure for each Z_post


# Author: JR Mihaljevic
# Date: Oct 2014

################################################################################################
################################################################################################

library(ggmcmc)

#Store z_ij values from best model run:
# Best model: 2covB (among others)

post.z <- NULL
post.z <- ggs(bundle_2covB, family="z")



######################################
#### SAMPLE FROM THE Z_POSTERIOR  ####
######################################

iter <- 500 # number of samples to draw from the posterior 

# Each chain has 1500 observations
# Choose which iterations will be used from the posterior
samples <- sample(c(1:1500), iter, replace=F)

# Create storage for the output
z.post <- array(0, dim=c(Nsite_2011, Nspecies_2011, iter))

for(i in 1:iter){
  # Randomly choose from which chain the sample will originate
  chain <- NULL
  chain <- sample(c(1:3), 1)
  
  # Create a subset vector of the values for all z[n, k]:
  subset <- NULL
  subset <- subset(post.z, Chain==chain & Iteration==samples[i])$value
  # Store this vector as a matrix
  z.mat <- NULL
  z.mat <- matrix(subset, nrow=Nsite_2011, ncol=Nspecies_2011, byrow=F)
  # Add matrix to the array
  z.post[, , i] <- z.mat
}

# Check to see if any sites ever have zero parasite species
# If so, remove the site entirely, for ease of next step:

for(i in 1:iter){
  if(any(rowSums(z.post[,,i])==0)){
    zeros <- NULL
    zeros <- which(rowSums(z.post[,,i])==0)
    z.post <- z.post[-zeros, , ]
  }
}

# How many sites are left?
nrow(z.post) # 61 (none removed)
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

####################################################
######## CALCULATE METACOMMUNITY METRICS  ##########
####################################################

# Ordinate all the z matrices:

library(metacom)
# Store the Ordinated Matrices:
z.ord <- array(0, dim=c(nrow(z.post), Nspecies_2011, iter))
# Store the statistics for Coherence and Turnover.
Coher <- array(0, dim=c(iter, 5))
colnames(Coher) <- c("Emb", "z", "pval", "sim.mean", "sim.sd")
Turn <- array(0, dim=c(iter, 5))
colnames(Turn) <- c("Repl", "z", "pval", "sim.mean", "sim.sd")
Bound <- array(0, dim=c(iter, 3))
colnames(Bound) <- c("index", "pval", "df")

pb <- txtProgressBar(min = 0, max = iter, style = 3)
for(i in 1:iter){
  meta <- NULL
  meta <- Metacommunity(z.post[, , i], method="r1", sims=1000, scores=1)
  z.ord[, , i] <- meta[[1]]
  Coher[i, ] <- as.numeric(as.character(meta[[2]][1:5, ]))
  Turn[i, ] <- as.numeric(as.character(meta[[3]][1:5, ]))
  for(j in 1:3){
    Bound[i, j] <- meta[[4]][1,j]
  }
  setTxtProgressBar(pb, i)
}
# Determine Structure for each:
Structure <- NULL

for(i in 1:iter){
  if(Coher[i, 3] > 0.05){
    Structure[i] <- "Random"
  }
  
  if(Coher[i, 3] <= 0.05 & Coher[i, 1] > Coher[i, 4]){
    Structure[i] <- "Checkerboard"
  }
  
  if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
       Turn[i, 3] > 0.05 & Turn[i, 1] < Turn[i, 4]){
    Structure[i] <- "Quasi-Nested"
  }
  
  if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
       Turn[i, 3] > 0.05 & Turn[i, 1] > Turn[i, 4] &
       Bound[i, 1] > 1 & Bound[i, 2] <= 0.05){
    Structure[i] <- "Quasi-Clementsian"
  }
  
  if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
       Turn[i, 3] > 0.05 & Turn[i, 1] > Turn[i, 4] &
       Bound[i, 1] < 1 & Bound[i, 2] <= 0.05){
    Structure[i] <- "Quasi-EvenSpaced"
  }
  
  if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
       Turn[i, 3] > 0.05 & Turn[i, 1] > Turn[i, 4] &
       Bound[i, 2] > 0.05){
    Structure[i] <- "Quasi-Gleasonian"
  }
  
  if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
       Turn[i, 3] <= 0.05 & Turn[i, 1] < Turn[i, 4]){
    Structure[i] <- "Nested"
  }
  
  if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
       Turn[i, 3] <= 0.05 & Turn[i, 1] > Turn[i, 4] &
       Bound[i, 1] > 1 & Bound[i, 2] <= 0.05){
    Structure[i] <- "Clementsian"
  }
  
  if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
       Turn[i, 3] <= 0.05 & Turn[i, 1] > Turn[i, 4] &
       Bound[i, 1] < 1 & Bound[i, 2] <= 0.05){
    Structure[i] <- "EvenSpaced"
  }
  
  if(Coher[i, 3] <= 0.05 & Coher[i, 1] < Coher[i, 4] & 
       Turn[i, 3] <= 0.05 & Turn[i, 1] > Turn[i, 4] &
       Bound[i, 2] > 0.05){
    Structure[i] <- "Gleasonian"
  }
}

###################################
#### POSTERIORS OF EMS METRICS ####
###################################
library(ggplot2)

Coher.plot <- ggplot(data.frame(Coher), aes(x=z))+
  geom_histogram(binwidth=0.3)+
  geom_histogram(binwidth=0.3, data=data.frame(Coher[which(Coher[, 3]>0.05), ]), color="grey")+
  labs(x="", y="")+
  theme_classic()+
  geom_vline(xintercept=1.96, linetype=5)+ # Significance cut-off
  scale_y_continuous(limits=c(0, 120), breaks=c(0, 60, 120))

quartz(height=3.5, width=3.5)
print(Coher.plot)

Turn.plot <- ggplot(data.frame(Turn), aes(x=z))+
  geom_histogram(binwidth=0.3)+
  geom_histogram(binwidth=0.3, data=data.frame(Turn[which(Coher[, 3]>0.05), ]), color="grey")+
  labs(x="", y="")+
  theme_classic()+
  geom_vline(xintercept=c(-1.96), linetype=5)+
  scale_x_continuous(breaks=c(-9, -3, 3))+
  scale_y_continuous(limits=c(0, 90), breaks=c(0, 45, 90))

quartz(height=3.5, width=3.5)
print(Turn.plot)

Bound.plot <- ggplot(data.frame(Bound), aes(x=index))+
  geom_histogram(binwidth=0.02)+
  geom_histogram(binwidth=0.02, data=data.frame(Bound[which(Coher[, 3]>0.05), ]), color="grey")+
  labs(x="", y="")+
  theme_classic()+
  scale_x_continuous(breaks=c(1, 1.3, 1.6))+
  scale_y_continuous(limits=c(0, 90), breaks=c(0, 45, 90))

quartz(height=3.5, width=3.5)
print(Bound.plot)

####################################
#### DISTRIBUTION OF STRUCTURES ####
####################################

head(Structure)
unique(Structure)

Structure2 <- Structure

# Change the order to something more logical/ideal for plotting
Structure2 <- replace(Structure2, which(Structure2=="Quasi-Clementsian"), "b_Quasi-Clem")
Structure2 <- replace(Structure2, which(Structure2=="Clementsian"), "a_Clem")
Structure2 <- replace(Structure2, which(Structure2=="Gleasonian"), "c_Gleas")
Structure2 <- replace(Structure2, which(Structure2=="Random"), "f_Rand")
Structure2 <- replace(Structure2, which(Structure2=="Quasi-Nested"), "e_Quasi-Nest")
Structure2 <- replace(Structure2, which(Structure2=="Quasi-Gleasonian"), "d_Quasi_Gleas")

Str.Counts <- as.data.frame(table(Structure2))
Str.Counts$ID <- rep("A", nrow(Str.Counts)) # Need an ID var for the bar plot

Str.Counts


frac <- ggplot(Str.Counts, aes(x=ID, y=Freq, fill=Structure2))+
  geom_bar(stat="identity")+
  theme_classic()+
  labs(x="", y="")+
  scale_fill_grey(breaks=paste(Str.Counts$Structure),
                  labels=c("Clementsian: 65.3% (Z structure)", "Quasi-Clementsian: 6.7%", "Gleasonian: 3.9%",
                           "Quasi-Gleasonian: 0.1%", "Quasi-Nested: 0.3%", 
                           "Random: 23.7% (Y structure)"))+
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_text(angle=90, hjust=0.5), 
        axis.text.x=element_blank(), legend.title=element_blank())

quartz(height=5, width=4)
print(frac)

#####################################
####   HEAT MAP OF Z-POSTERIOR   ####
#####################################

Z_heat <- mat.or.vec(nr=nrow(z.post), nc=Nspecies_2011)
for(i in 1:iter){
  Z_heat <- Z_heat + z.post[,,i]
}
Z_heat <- Z_heat/iter # Now each cell is a probability of occurrence averaged over z.post
rownames(Z_heat) <- as.factor(1:nrow(z.post))
colnames(Z_heat) <- as.factor(1:Nspecies_2011)
# Conduct a CCA and order the sites based on the first axes scores:
library(vegan)

Z_heat.CCA <- decorana(Z_heat, ira=0) # CHECK CONVERGENCE, ETC.
Z_heat.CCA.sites <- Z_heat.CCA$rproj[, 1]
Z_heat.CCA.spp <- Z_heat.CCA$cproj[, 1]

colnames(Z_heat) <- rownames(Yobs_2011)
Z_heat_ord <- Z_heat[order(Z_heat.CCA.sites, decreasing = FALSE), 
                     order(Z_heat.CCA.spp, decreasing = FALSE)]

x11(height=6, width=4)
#quartz(height=6, width=4)
print(Matrix_HeatMap_NMDS(Z_heat_ord, xlab="", ylab=""))
#################################################################################
#### LOOK AT PROB DETECTION AND SPECIES-SPECIFIC COVARIATE EFFECT POSTERIORS ####
#################################################################################
library(rjags)
library(ggmcmc)

mod2 <- jags.model(file = "OccMod_SingleYear.txt", 
                   data = jags_d, n.chains = 3, n.adapt=1000,
                   inits = list(z=zinit))
update(mod2, n.iter=5000)
out2 <- coda.samples(mod2, n.iter = 10000, variable.names = c("p", "b"), thin=10)

post.p <- NULL
post.p <- ggs(out2, family="p")
Rhat_p <- get_Rhat(post.p)
head(post.p)

post.b <- NULL
post.b <- ggs(out2, family="b")
Rhat_b <- get_Rhat(post.b)
head(post.b)

###########################      
#######  DET. PROB  ####### 
###########################  
library(mcmcplots)

quartz(height=5, width=7)
caterplot(out2, parms="p", col="black", val.lim=c(0, 1),
          labels=paste(LETTERS[1:N], sep=""), )
caterpoints(p0)


##############################      
#######  COV. EFFECTS  ####### 
############################## 

quartz(height=5, width=7)
caterplot(out2, parms="b", col="black")
caterpoints(b.spp)
