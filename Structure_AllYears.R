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
# Best model: 2cov (among others)

post.z <- NULL
post.z <- ggs(bundle_2cov, family="z")



######################################
#### SAMPLE FROM THE Z_POSTERIOR  ####
######################################

iter <- 500 # number of samples to draw from the posterior 

# Each chain has 1500 observations
# Choose which iterations will be used from the posterior
samples <- sample(c(1:1500), iter, replace=F)

# Create storage for the output
z.post <- array(0, dim=c(Nsite_all, Nspecies_all, iter))

for(i in 1:iter){
  # Randomly choose from which chain the sample will originate
  chain <- NULL
  chain <- sample(c(1:3), 1)
  
  # Create a subset vector of the values for all z[n, k]:
  subset <- NULL
  subset <- subset(post.z, Chain==chain & Iteration==samples[i])$value
  # Store this vector as a matrix
  z.mat <- NULL
  z.mat <- matrix(subset, nrow=Nsite_all, ncol=Nspecies_all, byrow=F)
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
nrow(z.post) # 261 (5 removed)

# Remove PinW, and Clin because almost never found...
z.post <- z.post[,-c(2,9),]
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

####################################################
######## CALCULATE METACOMMUNITY METRICS  ##########
####################################################

# Ordinate all the z matrices:

library(metacom)
# Store the Ordinated Matrices:
z.ord <- array(0, dim=c(nrow(z.post), ncol(z.post), iter))
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
  meta <- Metacommunity(z.post[, , i], method="r1", sims=1000, scores=2, 
                        allow.empty=T)
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
Structure2 <- replace(Structure2, which(Structure2=="Nested"), "d_Nest")
Structure2 <- replace(Structure2, which(Structure2=="Random"), "f_Rand")
Structure2 <- replace(Structure2, which(Structure2=="Quasi-Nested"), "c_QNest")
Structure2 <- replace(Structure2, which(Structure2=="Checkerboard"), "e_Check")
Structure2 <- replace(Structure2, which(Structure2=="Clementsian"), "b_Clem")
Structure2 <- replace(Structure2, which(Structure2=="Quasi-Clementsian"), "a_QClem")

Str.Counts <- as.data.frame(table(Structure2))
Str.Counts$ID <- rep("A", nrow(Str.Counts)) # Need an ID var for the bar plot

Str.Counts


frac <- ggplot(Str.Counts, aes(x=ID, y=Freq, fill=Structure2))+
  geom_bar(stat="identity")+
  theme_classic()+
  labs(x="", y="")+
  scale_fill_grey(labels=c("Quasi-Clementsian (45.4%)", "Clementsian (6.8%)",
                           "Quasi-Nested (32.0%)", "Nested (1.6%)",
                           "Checkerboard (0.2%)","Random (14.0%)"))+
  scale_y_continuous(breaks=c(0,250,500))+
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_text(angle=90, hjust=0.5), 
        axis.text.x=element_blank(), legend.title=element_blank())

quartz(height=5, width=4)
print(frac)

#####################################
####   HEAT MAP OF Z-POSTERIOR   ####
#####################################

Z_heat <- mat.or.vec(nr=nrow(z.post), nc=ncol(z.post))
for(i in 1:iter){
  Z_heat <- Z_heat + z.post[,,i]
}
Z_heat <- Z_heat/iter # Now each cell is a probability of occurrence averaged over z.post
rownames(Z_heat) <- as.factor(1:nrow(z.post))
colnames(Z_heat) <- as.factor(1:ncol(z.post))
# Conduct a CCA and order the sites based on the first axes scores:
library(vegan)

Z_heat.CCA <- decorana(Z_heat, ira=0) # CHECK CONVERGENCE, ETC.
Z_heat.CCA.sites <- Z_heat.CCA$rproj[, 2]
Z_heat.CCA.spp <- Z_heat.CCA$cproj[, 2]

colnames(Z_heat) <- c("Alar", "Fib", "Glob", "Echi", "Mano", "Nyct", "Opal", "Rib")
Z_heat_ord <- Z_heat[order(Z_heat.CCA.sites, decreasing = FALSE), 
                     order(Z_heat.CCA.spp, decreasing = FALSE)]

#x11(height=6, width=4)
quartz(height=6, width=3)
print(Matrix_HeatMap_NMDS(Z_heat_ord, xlab="", ylab=""))
#################################################################################
#### LOOK AT PROB DETECTION AND SPECIES-SPECIFIC COVARIATE EFFECT POSTERIORS ####
#################################################################################

library(mcmcplots)

plotThese <- paste("p.detect[", c(1,3:8,10), "]", sep="")
quartz(height=4, width=5)
caterplot(bundle_2cov, parms=plotThese, col="black", val.lim=c(0,1),
          labels=c("Alar", "Fib", "Glob", "Echi", "Mano", "Nyct", "Opal", "Rib"),
          horizontal=T)


betalabs <- c(paste(c("Alar", "Fib", "Glob", "Echi", "Mano", "Nyct", "Opal", "Rib"),1,sep=", "),
              paste(c("Alar", "Fib", "Glob", "Echi", "Mano", "Nyct", "Opal", "Rib"),2,sep=", "))

plotThese2 <- c(paste("betas[", c(1,3:8,10), ",", 1, "]", sep=""),
                paste("betas[", c(1,3:8,10), ",", 2, "]", sep=""))
quartz(height=5, width=7)
caterplot(bundle_2cov, parms=plotThese2, col="black", labels=betalabs)

quartz(height=4, width=6)
caterplot(bundle_2cov, parms="mean.beta.post", col="black",
          labels=c("Amph_RA1", "Amph_RA2"), las=0,
          horizontal=F)

quartz(height=4, width=6)
caterplot(bundle_2cov, parms="sd.beta.post", col="black",
          labels=c("Amph_RA1", "Amph_RA2"), las=0,
          horizontal=F)
