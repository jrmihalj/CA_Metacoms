# Use basic reciprocal averaging to capture
# amphibian and snail host community structure.
# Based on pres/abs:

host_data <- read.csv(file.choose(), header=T) #HostData.csv
head(host_data)

# Make a site_yr column:
host_data$site.yr <- paste(host_data$site, "_", host_data$year, sep="")

amphs <- host_data[, c(1:8, 16)]
snails <- host_data[, c(1:2, 10:14, 16)]

# Sort by year to match Y.obs parasite occurrences:
amphs <- amphs[order(amphs$year), ]
snails <- snails[order(snails$year), ]

# Subset only the sites to be included:
amphs.included <- data.frame()
snails.included <- data.frame()

for(i in 1:ncol(Y.obs)){
  sub1 <- NULL; sub2 <- NULL
  sub1 <- subset(amphs, site.yr==colnames(Y.obs)[i])
  sub2 <- subset(snails, site.yr==colnames(Y.obs)[i])
  amphs.included <- rbind(amphs.included, sub1)
  snails.included <- rbind(snails.included, sub2)
}

rownames(amphs.included) <- amphs.included$site.yr
rownames(snails.included) <- snails.included$site.yr

# Remove unused columns, then convert to matrix:
amphs.included <- amphs.included[, c(3:8)]
snails.included <- snails.included[, c(3:7)]

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
# Conduct RDA:
library(vegan)
Amph1 <- decorana(Amph.Mats[[1]], ira=1)
print(Amph1)
Amph1$rproj
Amph2 <- decorana(Amph.Mats[[2]], ira=1)
Amph3 <- decorana(Amph.Mats[[3]], ira=1)
Amph4 <- decorana(Amph.Mats[[4]], ira=1)
Amph5_all <- decorana(Amph.Mats[[5]], ira=1)

# Combine into one data.frame
Amph_RA <- data.frame(rbind(Amph1$rproj[,1:2], Amph2$rproj[,1:2], 
                              Amph3$rproj[,1:2], Amph4$rproj[,1:2]))
colnames(Amph_RA) <- c("Amph_RA1", "Amph_RA2")

Amph_RA_all <- data.frame(Amph5_all$rproj[,1:2])
colnames(Amph_RA_all) <- c("Amph_RA1", "Amph_RA2")


### SNAILS ###

Snails.Mats[[1]] <- Snails.Mats[[1]][c(-which(rowSums(Snails.Mats[[1]])==0)), ]
Snails1 <- decorana(Snails.Mats[[1]], ira=1)
print(Snails1)
Snails1$rproj

Snails.Mats[[2]] <- Snails.Mats[[2]][-which(rowSums(Snails.Mats[[2]])==0), ]
Snails2 <- decorana(Snails.Mats[[2]], ira=1)

Snails.Mats[[3]] <- Snails.Mats[[3]][-which(rowSums(Snails.Mats[[3]])==0), ]
Snails3 <- decorana(Snails.Mats[[3]], ira=1)

Snails.Mats[[4]] <- Snails.Mats[[4]][-which(rowSums(Snails.Mats[[4]])==0), ]
Snails4 <- decorana(Snails.Mats[[4]], ira=1)

Snails.Mats[[5]] <- Snails.Mats[[5]][-which(rowSums(Snails.Mats[[5]])==0), ]
Snails5 <- decorana(Snails.Mats[[5]], ira=1)

# Combine into one data.frame
Snails_RA <- data.frame(rbind(Snails1$rproj[,1:2], Snails2$rproj[,1:2], 
                                Snails3$rproj[,1:2], Snails4$rproj[,1:2]))
colnames(Snails_RA) <- c("Snails_RA1", "Snails_RA2")
Snails_RA$site.yr <- rownames(Snails_RA)

Snails_RA_all <- data.frame(Snails5$rproj[,1:2])
colnames(Snails_RA_all) <- c("Snails_RA1", "Snails_RA2")
Snails_RA_all$site.yr <- rownames(Snails_RA_all)

# Fill in the blanks for Snails (the sites that had no snails present)

site.yr <- colnames(Y.obs)
Snails_RA_fixed <- data.frame(site.yr)
Snails_RA_fixed$Snails_RA1 <- rep(0, nrow(Snails_RA_fixed))
Snails_RA_fixed$Snails_RA2 <- rep(0, nrow(Snails_RA_fixed))

for(i in 1:ncol(Y.obs)){
  sub <- NULL
  sub <- subset(Snails_RA, site.yr==colnames(Y.obs)[i])
  if(nrow(sub)==0){
    Snails_RA_fixed[i, 2:3] <- c(NA, NA)
  }else{
    Snails_RA_fixed[i, 2:3] <- Snails_RA[which(Snails_RA$site.yr==colnames(Y.obs)[i]), 1:2]
  }
}

# Fill in the blanks for Snails (all years)

site.yr <- colnames(Y.obs)
Snails_RA_all_fixed <- data.frame(site.yr)
Snails_RA_all_fixed$Snails_RA1 <- rep(0, nrow(Snails_RA_all_fixed))
Snails_RA_all_fixed$Snails_RA2 <- rep(0, nrow(Snails_RA_all_fixed))

for(i in 1:ncol(Y.obs)){
  sub <- NULL
  sub <- subset(Snails_RA_all, site.yr==colnames(Y.obs)[i])
  if(nrow(sub)==0){
    Snails_RA_all_fixed[i, 2:3] <- c(NA, NA)
  }else{
    Snails_RA_all_fixed[i, 2:3] <- Snails_RA_all[which(Snails_RA_all$site.yr==colnames(Y.obs)[i]), 1:2]
  }
}

