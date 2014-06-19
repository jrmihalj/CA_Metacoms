# Using non-metric multi-dimensional scaling (NMDS) to better capture
# amphibian and snail host community structure:
# Based on abund:

host_abund <- read.csv(file.choose(), header=T) #HostAbund.csv
head(host_abund)

amphs <- host_abund[, c(1:9)]
snails <- host_abund[, c(1:3, 10:14)]

# Sort by year to match Y.obs parasite occurrences:
amphs <- amphs[order(amphs$Year), ]
snails <- snails[order(snails$Year), ]

# Subset only the sites to be included:
amphs.included <- data.frame()
snails.included <- data.frame()

for(i in 1:ncol(Y.obs)){
  sub1 <- NULL; sub2 <- NULL
  sub1 <- subset(amphs, Site_yr==colnames(Y.obs)[i])
  sub2 <- subset(snails, Site_yr==colnames(Y.obs)[i])
  amphs.included <- rbind(amphs.included, sub1)
  snails.included <- rbind(snails.included, sub2)
}

rownames(amphs.included) <- amphs.included$Site_yr
rownames(snails.included) <- snails.included$Site_yr

# Remove unused columns, then convert to matrix:
amphs.included <- amphs.included[, c(4:9)]
snails.included <- snails.included[, c(4:8)]

amphs.mat <- as.matrix(amphs.included)
snails.mat <- as.matrix(snails.included)

# Make a list with a matrix per year of data:
Amph.Mats <- list()
Amph.Mats[[1]] <- amphs.mat[1:77, ]
Amph.Mats[[2]] <- amphs.mat[78:175, ]
Amph.Mats[[3]] <- amphs.mat[176:236, ]
Amph.Mats[[4]] <- amphs.mat[237:266, ]
Amph.Mats[[5]] <- amphs.mat

Snails.Mats <- list()
Snails.Mats[[1]] <- snails.mat[1:77, ]
Snails.Mats[[2]] <- snails.mat[78:175, ]
Snails.Mats[[3]] <- snails.mat[176:236, ]
Snails.Mats[[4]] <- snails.mat[237:266, ]
Snails.Mats[[5]] <- snails.mat

### AMPHIBIANS ###
# Make sure to remove rows that have no data
library(vegan)
which(rowSums(Amph.Mats[[4]])==0)

Amph.Mats[[1]] <- Amph.Mats[[1]][-which(rowSums(Amph.Mats[[1]])==0), ]
Amph1 <- metaMDS(Amph.Mats[[1]], k=2, trymax=100)
print(Amph1)
Amph1$points

Amph.Mats[[2]] <- Amph.Mats[[2]][-which(rowSums(Amph.Mats[[2]])==0), ]
Amph2 <- metaMDS(Amph.Mats[[2]], k=2, trymax=100)

Amph.Mats[[3]] <- Amph.Mats[[3]][-which(rowSums(Amph.Mats[[3]])==0), ]
Amph3 <- metaMDS(Amph.Mats[[3]], k=2, trymax=300)

Amph.Mats[[4]] <- Amph.Mats[[4]][-which(rowSums(Amph.Mats[[4]])==0), ]
Amph4 <- metaMDS(Amph.Mats[[4]], k=2, trymax=100)

Amph.Mats[[5]] <- Amph.Mats[[5]][-which(rowSums(Amph.Mats[[5]])==0), ]
Amph5_all <- metaMDS(Amph.Mats[[5]], k=3, trymax=200) # Need 3 dimensions here


# Combine into one data.frame
Amph_NMDS <- data.frame(rbind(Amph1$points, Amph2$points[,1:2], 
                              Amph3$points, Amph4$points))
colnames(Amph_NMDS) <- c("Amph_MDS1", "Amph_MDS2")

Amph_NMDS_all <- data.frame(Amph5_all$points[,1:2]) # Just first two axes
colnames(Amph_NMDS_all) <-  c("Amph_MDS1", "Amph_MDS2")

### SNAILS ###
which(rowSums(Snails.Mats[[4]])==0)

Snails.Mats[[1]] <- Snails.Mats[[1]][-which(rowSums(Snails.Mats[[1]])==0), ]
Snails1 <- metaMDS(Snails.Mats[[1]], k=2, trymax=100)
print(Snails1)
Snails1$points

Snails.Mats[[2]] <- Snails.Mats[[2]][-which(rowSums(Snails.Mats[[2]])==0), ]
Snails2 <- metaMDS(Snails.Mats[[2]], k=2, trymax=100)

Snails.Mats[[3]] <- Snails.Mats[[3]][-which(rowSums(Snails.Mats[[3]])==0), ]
Snails3 <- metaMDS(Snails.Mats[[3]], k=2, trymax=100)

Snails.Mats[[4]] <- Snails.Mats[[4]][-which(rowSums(Snails.Mats[[4]])==0), ]
Snails4 <- metaMDS(Snails.Mats[[4]], k=2, trymax=100)

Snails.Mats[[5]] <- Snails.Mats[[5]][-which(rowSums(Snails.Mats[[5]])==0), ]
Snails5_all <- metaMDS(Snails.Mats[[5]], k=2, trymax=100) # Only need 2 here


# Combine into one data.frame
Snails_NMDS <- data.frame(rbind(Snails1$points, Snails2$points, 
                              Snails3$points, Snails4$points))
colnames(Snails_NMDS) <- c("Snails_MDS1", "Snails_MDS2")

Snails_NMDS_all <- data.frame(Snails5_all$points)
colnames(Snails_NMDS_all) <- c("Snails_MDS1", "Snails_MDS2")


#############################################
#############################################

# Fill in the blanks (the sites that had no data)

site.yr <- colnames(Y.obs)
Snails_NMDS$site.yr <- rownames(Snails_NMDS)
Amph_NMDS$site.yr <- rownames(Amph_NMDS)

Snails_NMDS_fixed <- data.frame(site.yr)
Amph_NMDS_fixed <- data.frame(site.yr)

Snails_NMDS_fixed$Snails_MDS1 <- rep(0, nrow(Snails_NMDS_fixed))
Snails_NMDS_fixed$Snails_MDS2 <- rep(0, nrow(Snails_NMDS_fixed))
Amph_NMDS_fixed$Amph_MDS1 <- rep(0, nrow(Amph_NMDS_fixed))
Amph_NMDS_fixed$Amph_MDS2 <- rep(0, nrow(Amph_NMDS_fixed))

for(i in 1:ncol(Y.obs)){
  sub1 <- NULL; sub2 <- NULL
  sub1 <- subset(Snails_NMDS, site.yr==colnames(Y.obs)[i])
  sub2 <- subset(Amph_NMDS, site.yr==colnames(Y.obs)[i])
  if(nrow(sub1)==0){
    Snails_NMDS_fixed[i, 2:3] <- c(NA, NA)
  }else{
    Snails_NMDS_fixed[i, 2:3] <- Snails_NMDS[which(Snails_NMDS$site.yr==colnames(Y.obs)[i]), 1:2]
  }
  
  if(nrow(sub2)==0){
    Amph_NMDS_fixed[i, 2:3] <- c(NA, NA)
  }else{
    Amph_NMDS_fixed[i, 2:3] <- Amph_NMDS[which(Amph_NMDS$site.yr==colnames(Y.obs)[i]), 1:2]
  }
}

#############################################
#############################################

# Fill in the blanks (the sites that had no data)...ALL YEARS

site.yr <- colnames(Y.obs)
Snails_NMDS_all$site.yr <- rownames(Snails_NMDS_all)
Amph_NMDS_all$site.yr <- rownames(Amph_NMDS_all)

Snails_NMDS_all_fixed <- data.frame(site.yr)
Amph_NMDS_all_fixed <- data.frame(site.yr)

Snails_NMDS_all_fixed$Snails_MDS1 <- rep(0, nrow(Snails_NMDS_all_fixed))
Snails_NMDS_all_fixed$Snails_MDS2 <- rep(0, nrow(Snails_NMDS_all_fixed))
Amph_NMDS_all_fixed$Amph_MDS1 <- rep(0, nrow(Amph_NMDS_all_fixed))
Amph_NMDS_all_fixed$Amph_MDS2 <- rep(0, nrow(Amph_NMDS_all_fixed))

for(i in 1:ncol(Y.obs)){
  sub1 <- NULL; sub2 <- NULL
  sub1 <- subset(Snails_NMDS_all, site.yr==colnames(Y.obs)[i])
  sub2 <- subset(Amph_NMDS_all, site.yr==colnames(Y.obs)[i])
  if(nrow(sub1)==0){
    Snails_NMDS_all_fixed[i, 2:3] <- c(NA, NA)
  }else{
    Snails_NMDS_all_fixed[i, 2:3] <- Snails_NMDS_all[which(Snails_NMDS_all$site.yr==colnames(Y.obs)[i]), 1:2]
  }
  
  if(nrow(sub2)==0){
    Amph_NMDS_all_fixed[i, 2:3] <- c(NA, NA)
  }else{
    Amph_NMDS_all_fixed[i, 2:3] <- Amph_NMDS_all[which(Amph_NMDS_all$site.yr==colnames(Y.obs)[i]), 1:2]
  }
}

