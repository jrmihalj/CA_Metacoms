# Make a heat map/column of Snail Richness and Canopy Cover
# for 2009:
sites_2009 <- rownames(meta_2009$Comm)
snail_rich <- NULL
canopy_cov <- NULL
Group <- c(rep("A", 18), rep("B", 62))
ID <- rep(1:length(sites_2009))

for(i in 1:length(sites_2009)){
  snail_rich[i] <- X[which(rownames(X)==sites_2009[i]), 11] #Snail Richness value
  canopy_cov[i] <- X[which(rownames(X)==sites_2009[i]), 5] #Canopy Cover value
}

df_snailrich <- data.frame(sites_2009, ID, snail_rich, Group)
df_canopy <- data.frame(sites_2009, ID, canopy_cov, Group)

gradient_2009 <- ggplot(df_canopy, aes(y=ID, x=Group))+
  geom_tile(aes(fill=canopy_cov), color=element_blank())+
  labs(x="Canopy Cover", y="")+
  scale_fill_gradient(low="#F0F0F0", high="black",
                      na.value="transparent")+
  theme_classic()+
  theme(axis.text=element_blank(), axis.ticks=element_blank())+
  theme(legend.position="none")+
  scale_y_reverse()

quartz(height=6, width=3)
print(gradient_2009)

# Ordinated Y matrix for 2009:
quartz(height=6, width=3)
print(Matrix_Plot(meta_2009$Comm, xlab="Parasite Species", ylab="Sites"))

##############################################################################################################
##############################################################################################################


# Make a heat map/column of Snail Richness and TotalN
# for 2010:
sites_2010 <- rownames(meta_2010$Comm)
snail_rich_10 <- NULL
totalN <- NULL
Group <- c(rep("A", 12), rep("B", 82))
ID <- rep(1:length(sites_2010))

for(i in 1:length(sites_2010)){
  snail_rich_10[i] <- X[which(rownames(X)==sites_2010[i]), 11] #Snail Richness value
  totalN[i] <- X[which(rownames(X)==sites_2010[i]), 8] #TotalN value
}

df_snailrich_10 <- data.frame(ID, snail_rich_10, Group)
df_totalN <- data.frame(ID, totalN, Group)

gradient_2010 <- ggplot(df_totalN, aes(y=ID, x=Group))+
  geom_tile(aes(fill=totalN), color=element_blank())+
  labs(x="Total N", y="")+
  scale_fill_gradient(low="#F0F0F0", high="black",
                      na.value="transparent")+
  theme_classic()+
  theme(axis.text=element_blank(), axis.ticks=element_blank())+
  theme(legend.position="none")+
  scale_y_reverse()

quartz(height=6, width=3)
print(gradient_2010)

# Ordinated Y matrix for 2009:
quartz(height=6, width=3)
print(Matrix_Plot(meta_2010$Comm, xlab="Parasite Species", ylab="Sites"))
