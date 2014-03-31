# CA metacommunity project
# Code to relate ordinated site-scores of PSRE individuals, amphibian hosts and snail hosts
# Year: 2010
# September 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# This is the data set that contains all of the environmental assessment for 2010 sites and hosts
head(alldata_indiv_2010)

#Rename the "value" variable to be more specific
alldata_indiv_2010$Ord_Host_value <- alldata_indiv_2010$value
#Remove the "value" column to avoid confusion
alldata_indiv_2010 <- alldata_indiv_2010[, -53]

#################### AMPH ###################
# For the amphibian site-scores, need to match sites to site.yr in the original amph dataset
head(sitescores_amph_2010)
sitescores_amph_2010$Ord_Amph_value <- sitescores_amph_2010$value # rename to something more useful
sitescores_amph_2010 <- sitescores_amph_2010[, -2] #remove 'values' variable

# The rownames correspond to the 'hosts' dataset
# Create a new column for site.yr
sitescores_amph_2010$site.yr <- vector(mode="character", length=nrow(sitescores_amph_2010))

# Match the site.yr to the 'hosts' dataset based on AmphSiteID
for (i in 1:nrow(sitescores_amph_2010)){
  sitescores_amph_2010$site.yr[i] <- hosts$site.yr[rownames(hosts)==sitescores_amph_2010$AmphSiteID[i]]
}

# Create a new column in the alldata_indiv_2010 dataset for Amph ordination score
alldata_indiv_2010$Ord_Amph_Value <- vector(mode="numeric", length=nrow(alldata_indiv_2010))

# Match up to the amph site scores
NAs <- NULL
for (i in 1:nrow(alldata_indiv_2010)){
  row <- NULL
  row <- which(sitescores_amph_2010$site.yr %in% alldata_indiv_2010[i, 2] == TRUE)
  if (sum(row) != 0){
    alldata_indiv_2010$Ord_Amph_Value[i] <- sitescores_amph_2010$Ord_Amph_value[row]
  } else {
    # Store which rows don't have an accompanying 'Assmt' in 'assess'
    NAs <- c(NAs, i)
  }   
}
NAs # No NAs so each row in alldata_indiv_2010 should have an amphib site score
head(alldata_indiv_2010) # yup

#################### SNAILS ###################
# For the snail site-scores, need to match sites to site.yr in the original snails dataset
head(sitescores_snails_2010)
sitescores_snails_2010$Ord_snails_value <- sitescores_snails_2010$value # rename to something more useful
sitescores_snails_2010 <- sitescores_snails_2010[, -2] #remove 'values' variable

# The rownames correspond to the 'snails' dataset
# Create a new column for site.yr
sitescores_snails_2010$site.yr <- vector(mode="character", length=nrow(sitescores_snails_2010))

# Match the site.yr to the 'snails' dataset based on snailsSiteID
for (i in 1:nrow(sitescores_snails_2010)){
  sitescores_snails_2010$site.yr[i] <- snails$site.yr[rownames(snails)==sitescores_snails_2010$snailsSiteID[i]]
}

# Create a new column in the alldata_indiv_2010 dataset for snails ordination score
alldata_indiv_2010$Ord_snails_Value <- vector(mode="numeric", length=nrow(alldata_indiv_2010))

# Match up to the snails site scores
NAs <- NULL
for (i in 1:nrow(alldata_indiv_2010)){
  row <- NULL
  row <- which(sitescores_snails_2010$site.yr %in% alldata_indiv_2010[i, 2] == TRUE)
  if (sum(row) != 0){
    alldata_indiv_2010$Ord_snails_Value[i] <- sitescores_snails_2010$Ord_snails_value[row]
  } else {
    # Store which rows don't have an accompanying 'Assmt' in 'assess'
    NAs <- c(NAs, i)
  }   
}
NAs 
length(NAs) #57 hosts do not have a snail site score because those sites had no snails
head(alldata_indiv_2010) 
alldata_indiv_2010$Ord_snails_Value

# Need to replace the zeros w/ NAs
alldata_indiv_2010$Ord_snails_Value[NAs] <- "NA"

# That seems to have changed the column to a character vector, so reverse that
alldata_indiv_2010$Ord_snails_Value <- as.numeric(as.character(alldata_indiv_2010$Ord_snails_Value))
# That worked

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Look at the correlation between ordination scores

plot(Ord_Host_value ~ Ord_Amph_Value, 
     data=alldata_indiv_2010[alldata_indiv_2010$Ord_Host_value < 5 & alldata_indiv_2010$Comp==1,])

plot(Ord_Host_value ~ Ord_Amph_Value, 
     data=alldata_indiv_2010[alldata_indiv_2010$Ord_Host_value < 5 & alldata_indiv_2010$Comp==2,])


plot(Ord_Host_value ~ Ord_snails_Value, 
     data=alldata_indiv_2010[alldata_indiv_2010$Ord_Host_value < 5 & alldata_indiv_2010$Comp==1,])

plot(Ord_Host_value ~ Ord_snails_Value, 
     data=alldata_indiv_2010[alldata_indiv_2010$Ord_Host_value < 5 & alldata_indiv_2010$Comp==2,])

# Negative relationship between snail and amphibian metacommunity scores
plot(Ord_snails_Value ~ Ord_Amph_Value, data=alldata_indiv_2010)

# all data
m1 <- lm(Ord_Host_value ~ Ord_Amph_Value*Ord_snails_Value, 
           data=alldata_indiv_2010[alldata_indiv_2010$Ord_Host_value < 5 ,])
summary(m1)

m2 <- update(m1, .~. -Ord_Amph_Value:Ord_snails_Value)
anova(m1, m2) # It's ok to drop the interaction

m3 <- lm(Ord_Host_value ~ Ord_Amph_Value + Ord_snails_Value, 
         data=alldata_indiv_2010[alldata_indiv_2010$Ord_Host_value < 5 ,])
summary(m3)
# Significant effect of amphibian structure

# Plot against the residuals
library(car)
crPlots(m3, terms= ~.)

# Do this with just compartments
# No patterns there...
m1a <- lm(Ord_Host_value ~ Ord_Amph_Value + Ord_snails_Value, 
         data=alldata_indiv_2010[alldata_indiv_2010$Ord_Host_value < 5 & alldata_indiv_2010$Comp==1,])
summary(m1a)

m1b <- lm(Ord_Host_value ~ Ord_Amph_Value + Ord_snails_Value, 
          data=alldata_indiv_2010[alldata_indiv_2010$Ord_Host_value < 5 & alldata_indiv_2010$Comp==2,])
summary(m1b)


# Add a random component of site.yr
# When random factor is included, no effects
library(nlme)
m1r <- lme(Ord_Host_value ~ Ord_Amph_Value+Ord_snails_Value, random= ~1 | site.yr,
         data=alldata_indiv_2010[alldata_indiv_2010$Ord_Host_value < 5 ,], method="REML",
           na.action="na.omit")
summary(m1r)
