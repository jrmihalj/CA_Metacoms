# Using non-metric multi-dimensional scaling (NMDS) to better capture
# amphibian and snail host community structure:

# First I'll use pres/abs:

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

for(i in 1:length(unique(PSRE.included$site.yr))){
  sub1 <- NULL; sub2 <- NULL
  sub1 <- subset(amphs, site.yr==unique(PSRE.included$site.yr)[i])
  sub2 <- subset(snails, site.yr==unique(PSRE.included$site.yr)[i])
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
Amph.Mats[[1]] <- amphs.mat[1:83, ]
Amph.Mats[[2]] <- amphs.mat[84:188, ]
Amph.Mats[[3]] <- amphs.mat[189:255, ]
Amph.Mats[[4]] <- amphs.mat[256:289, ]

Snails.Mats <- list()
Snails.Mats[[1]] <- snails.mat[1:83, ]
Snails.Mats[[2]] <- snails.mat[84:188, ]
Snails.Mats[[3]] <- snails.mat[189:255, ]
Snails.Mats[[4]] <- snails.mat[256:289, ]

### AMPHIBIANS ###
Amph1 <- metaMDS(Amph.Mats[[1]], k=2, trymax=100)
print(Amph1)
Amph1$points
Amph2 <- metaMDS(Amph.Mats[[2]], k=2, trymax=100)
Amph3 <- metaMDS(Amph.Mats[[3]], k=2, trymax=100)
Amph4 <- metaMDS(Amph.Mats[[4]], k=2, trymax=100)

# Combine into one data.frame
Amph_NMDS <- data.frame(rbind(Amph1$points, Amph2$points, 
                              Amph3$points, Amph4$points))
colnames(Amph_NMDS) <- c("Amph_MDS1", "Amph_MDS2")

### SNAILS ###

# CA-CORTE_2009 has missing values; must remove
which(rownames(Snails.Mats[[1]])=="CA-CORTE_2009") # Row 19
Snails.Mats[[1]] <- Snails.Mats[[1]][-19, ]
# Remove "Beaver_2009" because no snails present:
Snails.Mats[[1]] <- Snails.Mats[[1]][-3, ]

Snails1 <- metaMDS(Snails.Mats[[1]], k=2, trymax=100)
print(Snails1)
Snails1$points

# Some have zero sums:
which(rowSums(Snails.Mats[[2]])==0)
# Remove
Snails.Mats[[2]] <- Snails.Mats[[2]][-c(4,40,42,55,82,83,91,99,100), ]
Snails2 <- metaMDS(Snails.Mats[[2]], k=2, trymax=100)

# Some have zero sums:
which(rowSums(Snails.Mats[[3]])==0)
# Remove
Snails.Mats[[3]] <- Snails.Mats[[3]][-c(21,46,47,56), ]
Snails3 <- metaMDS(Snails.Mats[[3]], k=2, trymax=100)

# Some have zero sums:
which(rowSums(Snails.Mats[[4]])==0)
# Remove
Snails.Mats[[4]] <- Snails.Mats[[4]][-c(3,5,20,21,29), ]
Snails4 <- metaMDS(Snails.Mats[[4]], k=2, trymax=100)

# Combine into one data.frame
Snails_NMDS <- data.frame(rbind(Snails1$points, Snails2$points, 
                              Snails3$points, Snails4$points))
colnames(Snails_NMDS) <- c("Snails_MDS1", "Snails_MDS2")
Snails_NMDS$site.yr <- rownames(Snails_NMDS)

# Fill in the blanks for Snails (the sites that had no snails present)

site.yr <- unique(PSRE.included$site.yr)
Snails_NMDS_fixed <- data.frame(site.yr)
Snails_NMDS_fixed$Snails_MDS1 <- rep(0, nrow(Snails_NMDS_fixed))
Snails_NMDS_fixed$Snails_MDS2 <- rep(0, nrow(Snails_NMDS_fixed))

for(i in 1:length(unique(PSRE.included$site.yr))){
  sub <- NULL
  sub <- subset(Snails_NMDS, site.yr==unique(PSRE.included$site.yr)[i])
  if(nrow(sub)==0){
    Snails_NMDS_fixed[i, 2:3] <- c(NA, NA)
  }else{
    Snails_NMDS_fixed[i, 2:3] <- Snails_NMDS[which(Snails_NMDS$site.yr==unique(PSRE.included$site.yr)[i]), 1:2]
  }
}

