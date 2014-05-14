# California PSRE parasite metacommunities
# Using occupancy modeling in conjunction with the EMS framework

# Script Author: JR Mihaljevic
# Project co-authors: Bethany J. Hoye, Pieter T.J. Johnson

# Plotting curves for species-specific probability of occurrence

##############################################################################################################
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################

# Establish some useful functions:
Logit <- function(x){
  log(x) - log(1 - x)
}

AntiLogit <- function(x){
  exp(x) / (exp(x) + 1)
}

##############################################################################################################
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################

# 2009:
# First, snail richness:
SR_2009 <- seq(from=min(X[, 11], na.rm=T), to=max(X[, 11], na.rm=T), length.out=500)

# Mano:
lpOcc_2009_Mano <- median(subset(post.b0, Parameter=="b0[6,1]")$value) +
                    median(subset(post.b_Mano, Parameter=="b[6,11,1]")$value)*SR_2009
pOcc_2009_Mano <- AntiLogit(lpOcc_2009_Mano)

Mano_2009_SR_df <- data.frame(SR_2009, pOcc_2009_Mano)

plot_Mano_2009_SR <- ggplot(Mano_2009_SR_df, aes(x=SR_2009, y=pOcc_2009_Mano))+
  geom_line()+
  theme_classic()+
  labs(x="Snail Richness", y="P(Occurrence)")+
  scale_y_continuous(breaks=c(0.1,0.3,0.5), labels=c("0.10", "0.30", "0.50"))

quartz(height=2.5, width=2.5)
print(plot_Mano_2009_SR)
#------------------------------------------------------------------------------------------------------------#

# Canopy Cover
CC_2009 <- seq(from=min(X[, 5], na.rm=T), to=max(X[, 5], na.rm=T), length.out=500)

# Mano:
lpOcc_2009_Mano <- median(subset(post.b0, Parameter=="b0[6,1]")$value) +
  median(subset(post.b_Mano, Parameter=="b[6,5,1]")$value)*CC_2009
pOcc_2009_Mano <- AntiLogit(lpOcc_2009_Mano)

Mano_2009_CC_df <- data.frame(CC_2009, pOcc_2009_Mano)

# Alaria:
lpOcc_2009_Alaria <- median(subset(post.b0, Parameter=="b0[1,1]")$value) +
  median(subset(post.b_Alaria, Parameter=="b[1,5,1]")$value)*CC_2009
pOcc_2009_Alaria <- AntiLogit(lpOcc_2009_Alaria)

Alaria_2009_CC_df <- data.frame(CC_2009, pOcc_2009_Alaria)

plot_2009_CC <- ggplot(Mano_2009_CC_df, aes(x=CC_2009, y=pOcc_2009_Mano))+
  geom_line()+
  geom_line(data=Alaria_2009_CC_df, aes(y=pOcc_2009_Alaria), linetype=2)+
  theme_classic()+
  scale_y_continuous(breaks=c(0.15, 0.25, 0.35))+
  labs(x="Canopy Cover", y="P(Occurrence)")

quartz(height=2.5, width=2.5)
print(plot_2009_CC)

##############################################################################################################
#------------------------------------------------------------------------------------------------------------#
##############################################################################################################

# 2010:
# First, snail richness:
SR_2010 <- seq(from=min(X[, 11], na.rm=T), to=max(X[, 11], na.rm=T), length.out=500)

# Mano:
lpOcc_2010_Mano <- median(subset(post.b0, Parameter=="b0[6,2]")$value) +
                    median(subset(post.b_Mano, Parameter=="b[6,11,2]")$value)*SR_2010
pOcc_2010_Mano <- AntiLogit(lpOcc_2010_Mano)

Mano_2010_SR_df <- data.frame(SR_2010, pOcc_2010_Mano)

plot_Mano_2010_SR <- ggplot(Mano_2010_SR_df, aes(x=SR_2010, y=pOcc_2010_Mano))+
  geom_line()+
  theme_classic()+
  labs(x="Snail Richness", y="P(Occurrence)")+
  scale_y_continuous(breaks=c(0.1,0.3,0.5), labels=c("0.10", "0.30", "0.50"))

quartz(height=2.5, width=2.5)
print(plot_Mano_2010_SR)
#------------------------------------------------------------------------------------------------------------#

# TotalN:
TN_2010 <- seq(from=min(X[, 8], na.rm=T), to=max(X[, 8], na.rm=T), length.out=500)

# Mano:
lpOcc_2010_Mano <- median(subset(post.b0, Parameter=="b0[6,2]")$value) +
  median(subset(post.b_Mano, Parameter=="b[6,8,2]")$value)*TN_2010
pOcc_2010_Mano <- AntiLogit(lpOcc_2010_Mano)

Mano_2010_TN_df <- data.frame(TN_2010, pOcc_2010_Mano)

# Alaria:
lpOcc_2010_Glob <- median(subset(post.b0, Parameter=="b0[4,2]")$value) +
  median(subset(post.b_Glob, Parameter=="b[4,8,2]")$value)*TN_2010
pOcc_2010_Glob <- AntiLogit(lpOcc_2010_Glob)

Glob_2010_TN_df <- data.frame(TN_2010, pOcc_2010_Glob)

plot_2010_TN <- ggplot(Mano_2010_TN_df, aes(x=TN_2010, y=pOcc_2010_Mano))+
  geom_line()+
  geom_line(data=Glob_2010_TN_df, aes(y=pOcc_2010_Glob), linetype=3)+
  theme_classic()+
  labs(x="Total N", y="P(Occurrence)")+
  scale_y_continuous(breaks=c(0.1,0.4,0.7), labels=c("0.10", "0.40", "0.70"))

quartz(height=2.5, width=2.5)
print(plot_2010_TN)