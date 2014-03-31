# CA metacommunity project
# Code to relate ordinated site-scores of PSRE individuals, amphibian hosts and snail hosts
# Year: 2012
# September 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# Merge the site-scores from the 2010 Trema_only ordination with assessment data

alldata_indiv_2012 <- Match_Rows(sitedata=Indiv_Site_Ass_Data, 
                                       sitescores=sitescores_PSRE_2012)
head(alldata_indiv_2012)

#Rename the "value" variable to be more specific
alldata_indiv_2012$Ord_Host_value <- alldata_indiv_2012$value
#Remove the "value" column to avoid confusion
alldata_indiv_2012 <- alldata_indiv_2012[, -53]

#################### AMPH ###################
# For the amphibian site-scores, need to match sites to site.yr in the original amph dataset
head(sitescores_amph_2012)
sitescores_amph_2012$Ord_Amph_value <- sitescores_amph_2012$value # rename to something more useful
sitescores_amph_2012 <- sitescores_amph_2012[, -2] #remove 'values' variable

# The rownames correspond to the 'hosts' dataset
# Create a new column for site.yr
sitescores_amph_2012$site.yr <- vector(mode="character", length=nrow(sitescores_amph_2012))

# Match the site.yr to the 'hosts' dataset based on AmphSiteID
for (i in 1:nrow(sitescores_amph_2012)){
  sitescores_amph_2012$site.yr[i] <- hosts$site.yr[rownames(hosts)==sitescores_amph_2012$AmphSiteID[i]]
}

# Create a new column in the alldata_indiv_2012 dataset for Amph ordination score
alldata_indiv_2012$Ord_Amph_Value <- vector(mode="numeric", length=nrow(alldata_indiv_2012))

# Match up to the amph site scores
NAs <- NULL
for (i in 1:nrow(alldata_indiv_2012)){
  row <- NULL
  row <- which(sitescores_amph_2012$site.yr %in% alldata_indiv_2012[i, 2] == TRUE)
  if (sum(row) != 0){
    alldata_indiv_2012$Ord_Amph_Value[i] <- sitescores_amph_2012$Ord_Amph_value[row]
  } else {
    # Store which rows don't have an accompanying 'Assmt' in 'assess'
    NAs <- c(NAs, i)
  }   
}
NAs # No NAs so each row in alldata_indiv_2012 should have an amphib site score
head(alldata_indiv_2012) # yup

#################### SNAILS ###################
# For the snail site-scores, need to match sites to site.yr in the original snails dataset
head(sitescores_snails_2012)
sitescores_snails_2012$Ord_snails_value <- sitescores_snails_2012$value # rename to something more useful
sitescores_snails_2012 <- sitescores_snails_2012[, -2] #remove 'values' variable

# The rownames correspond to the 'snails' dataset
# Create a new column for site.yr
sitescores_snails_2012$site.yr <- vector(mode="character", length=nrow(sitescores_snails_2012))

# Match the site.yr to the 'snails' dataset based on snailsSiteID
for (i in 1:nrow(sitescores_snails_2012)){
  sitescores_snails_2012$site.yr[i] <- snails$site.yr[rownames(snails)==sitescores_snails_2012$snailsSiteID[i]]
}

# Create a new column in the alldata_indiv_2012 dataset for snails ordination score
alldata_indiv_2012$Ord_snails_Value <- vector(mode="numeric", length=nrow(alldata_indiv_2012))

# Match up to the snails site scores
NAs <- NULL
for (i in 1:nrow(alldata_indiv_2012)){
  row <- NULL
  row <- which(sitescores_snails_2012$site.yr %in% alldata_indiv_2012[i, 2] == TRUE)
  if (sum(row) != 0){
    alldata_indiv_2012$Ord_snails_Value[i] <- sitescores_snails_2012$Ord_snails_value[row]
  } else {
    # Store which rows don't have an accompanying 'Assmt' in 'assess'
    NAs <- c(NAs, i)
  }   
}
NAs # No NAs
head(alldata_indiv_2012) 
alldata_indiv_2012$Ord_snails_Value

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Look at the correlation between ordination scores
# Outliers above 1 and less than -1
plot(Ord_Host_value ~ Ord_Amph_Value, 
     data=subset(alldata_indiv_2012, Ord_Host_value < 1 & Ord_Host_value > -1))

plot(Ord_Host_value ~ Ord_snails_Value, 
     data=subset(alldata_indiv_2012, Ord_Host_value < 1 & Ord_Host_value > -1))

plot(Ord_snails_Value ~ Ord_Amph_Value, data=alldata_indiv_2012)

# all data
m1 <- lm(Ord_Host_value ~ Ord_Amph_Value*Ord_snails_Value, 
         data=subset(alldata_indiv_2012, Ord_Host_value < 1 & Ord_Host_value > -1))
summary(m1)

m2 <- update(m1, .~. -Ord_Amph_Value:Ord_snails_Value)
anova(m1, m2) # It's ok to drop the interaction

m3 <- lm(Ord_Host_value ~ Ord_Amph_Value + Ord_snails_Value, 
         data=subset(alldata_indiv_2012, Ord_Host_value < 1 & Ord_Host_value > -1))
summary(m3)
# Significant effect of snail structure, but R2 sucks...

# Plot against the residuals
library(car)
crPlots(m3, terms= ~.)


# Add a random component of site.yr
# When random factor is included, no effects
library(nlme)
m1r <- lme(Ord_Host_value ~ Ord_Amph_Value+Ord_snails_Value, random= ~1 | site.yr,
           data=subset(alldata_indiv_2012, Ord_Host_value < 1 & Ord_Host_value > -1), 
           method="REML",
           na.action="na.omit")
summary(m1r)
