# CA metacommunity project
# Code to relate ordinated site-scores of PSRE individuals, amphibian hosts and snail hosts
# Using the scores from the TREMA_only Ordination.
# Year: 2012
# September 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Merge the site-scores from the 2012 Trema_only ordination with assessment data

alldata_indiv_2012_trema <- Match_Rows(sitedata=Indiv_Site_Ass_Data, 
                                       sitescores=sitescores_PSRE_2012_trema)
head(alldata_indiv_2012_trema)

#Rename the "value" variable to be more specific
alldata_indiv_2012_trema$Ord_Host_value <- alldata_indiv_2012_trema$value
#Remove the "value" column to avoid confusion
alldata_indiv_2012_trema <- alldata_indiv_2012_trema[, -53]

#################### AMPH ###################
# Create a new column in the alldata_indiv_2012 dataset for Amph ordination score
alldata_indiv_2012_trema$Ord_Amph_Value <- vector(mode="numeric", length=nrow(alldata_indiv_2012_trema))

# Match up to the amph site scores
NAs <- NULL
for (i in 1:nrow(alldata_indiv_2012_trema)){
  row <- NULL
  row <- which(sitescores_amph_2012$site.yr %in% alldata_indiv_2012_trema[i, 2] == TRUE)
  if (sum(row) != 0){
    alldata_indiv_2012_trema$Ord_Amph_Value[i] <- sitescores_amph_2012$Ord_Amph_value[row]
  } else {
    # Store which rows don't have an accompanying 'Assmt' in 'assess'
    NAs <- c(NAs, i)
  }   
}
NAs # No NAs so each row in alldata_indiv_2012_trema should have an amphib site score
head(alldata_indiv_2012_trema) # yup

#################### SNAILS ###################
alldata_indiv_2012_trema$Ord_snails_Value <- vector(mode="numeric", length=nrow(alldata_indiv_2012_trema))

# Match up to the snails site scores
NAs <- NULL
for (i in 1:nrow(alldata_indiv_2012_trema)){
  row <- NULL
  row <- which(sitescores_snails_2012$site.yr %in% alldata_indiv_2012_trema[i, 2] == TRUE)
  if (sum(row) != 0){
    alldata_indiv_2012_trema$Ord_snails_Value[i] <- sitescores_snails_2012$Ord_snails_value[row]
  } else {
    # Store which rows don't have an accompanying 'Assmt' in 'assess'
    NAs <- c(NAs, i)
  }   
}
NAs # No NA's

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Look at the correlation between ordination scores
# Outliers above 1 and less than -1
plot(Ord_Host_value ~ Ord_Amph_Value, 
     data=subset(alldata_indiv_2012_trema, Ord_Host_value < 1 & Ord_Host_value > -1))

plot(Ord_Host_value ~ Ord_snails_Value, 
     data=subset(alldata_indiv_2012_trema, Ord_Host_value < 1 & Ord_Host_value > -1))

plot(Ord_snails_Value ~ Ord_Amph_Value, data=alldata_indiv_2012_trema)

# all data
m1 <- lm(Ord_Host_value ~ Ord_Amph_Value*Ord_snails_Value, 
         data=subset(alldata_indiv_2012_trema, Ord_Host_value < 1 & Ord_Host_value > -1))
summary(m1)

m2 <- update(m1, .~. -Ord_Amph_Value:Ord_snails_Value)
anova(m1, m2) # It's ok to drop the interaction

m3 <- lm(Ord_Host_value ~ Ord_Amph_Value + Ord_snails_Value, 
         data=subset(alldata_indiv_2012_trema, Ord_Host_value < 1 & Ord_Host_value > -1))
summary(m3)
# Significant effect of snail structure, but R2 sucks...

# Plot against the residuals
library(car)
crPlots(m3, terms= ~.)


# Add a random component of site.yr
# When random factor is included, no effects
library(nlme)
m1r <- lme(Ord_Host_value ~ Ord_Amph_Value+Ord_snails_Value, random= ~1 | site.yr,
           data=subset(alldata_indiv_2012_trema, Ord_Host_value < 1 & Ord_Host_value > -1), 
           method="REML",
           na.action="na.omit")
summary(m1r)