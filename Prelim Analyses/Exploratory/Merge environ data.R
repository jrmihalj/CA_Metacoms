# This code will be used to match up individual-level and site-level data (e.g. site-level)
# environmental assessment data. This combined data will be used in analyses to associate
# metacommunity structure with environmental covariates for each paritioned dataset. 
# August 2013
# Author: JR Mihaljevic

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Start with the individual-level, full dataset
# I am going to create a data frame with all of the host-level data and the environmental
# assessment data. This includes non-PSRE hosts because this matches the ID numbers on the 
# ordinated incidence matrices. 

head(indiv_allhosts)
# I just want the first eight columns, which give me relevant, indiv-level data

indiv_allhosts2 <- indiv_allhosts[, c(1:8)]

# Import the site-level environmental data
# file name: 'Site_envir_data.csv'
assess <- read.csv(file.choose(), header=T, na.string="")
head(assess)
ncol(assess)

# Add columns and rows to indiv_allhosts2 that will encompass the 'assess' data
site.dataframe <- mat.or.vec(nr=nrow(indiv_allhosts2), nc=ncol(assess))
indiv_allhosts2 <- cbind(indiv_allhosts2, site.dataframe)
colnames(indiv_allhosts2)[9:51] <- colnames(assess)
head(indiv_allhosts2)

# Check each row in 'indiv_allhosts2' for the row in 'assess' that is relevant,
# then paste it into 'inidv_allhosts2'

NAs <- NULL
for ( i in 1:nrow(indiv_allhosts2)){
  # Store the rows in 'assess' for which the assessment data is applicable
  row <- NULL
  row <- which(assess$Assmt %in% indiv_allhosts2[i, 2] == TRUE)
  if (sum(row) != 0){
    indiv_allhosts2[i , 9:51] <- assess[row ,]  
  } else {
    # Store which rows don't have an accompanying 'Assmt' in 'assess'
    NAs <- c(NAs, i)
  }   
}
# Which rows have NAs?
missing <- indiv_allhosts2[NAs, ]
head(missing)


# Seems like for these missing ones, the assessment date was just a few days off from the 
# date the animals were collected. Therefore, I will add the data based on 'site.yr',
# instead of 'Assmt' name. 

for ( i in 1:length(NAs)){
  # Store the rows in 'assess' for which the assessment data is applicable
  row <- NULL
  row <- which(assess$Site_yr %in% indiv_allhosts2[NAs[i], 1] == TRUE)
  if (sum(row) != 0){
    indiv_allhosts2[i , 9:51] <- assess[row ,]  
  } else {
    # Store which rows don't have an accompanying 'Assmt' in 'assess'
    NAs <- c(NAs, i)
  }   
}


# Rename data.frame to something more useful
Indiv_Site_Ass_Data <- indiv_allhosts2
Indiv_Site_Ass_Data$HostID <- as.numeric(as.character(rownames(Indiv_Site_Ass_Data)))
summary(Indiv_Site_Ass_Data)

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Now make a function that will match up the rows in the ordinated matrices to the 
# associated rows in this newly created data frame: 'Indiv_Site_Ass_data'

Match_Rows <- function(sitescores, sitedata){
  
  sub <- subset(sitedata, HostID %in% sitescores$HostID)
  
  merged <- merge(sub, sitescores, by="HostID")
  return(merged) 
}

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# 2010 Inidividual-level PSRE data
alldata_indiv_2010 <- Match_Rows(sitedata=Indiv_Site_Ass_Data, 
                                  sitescores=sitescores_PSRE_2010)
head(alldata_indiv_2010)
# The 'value' is the site score

# subset out the sites that have really large influence?
wo_large <- subset(alldata_indiv_2010, value < 5)


# Pairs plot
pairs(~value+no_malf+SVL+Lat+Long+Elev+Slope+Aspect+hydro+AGG+DEV+FOR+SSG+WET, data=wo_large)
# no_malf looks pretty interesting
pairs(~value+area+depth+pH+tempW+Troph+veg_s+fish, data=wo_large)
pairs(~value+Cond+TDS+Salin+DO+DOmg+DOC+NH4+DON+PO4+DOP+NO3, data=wo_large)
pairs(~value+TATO+RADR+BUBO+AMCA+RACA+PSRE+PHySA+LyMN+HELI+GYRA+RADIX, data=wo_large)
# TATO and PSRE look interesting

###### need to look at total density/richness of hosts at a site! ######

wo_large$no_malf <- as.numeric(as.character(wo_large$no_malf))
plot(value~no_malf, type="p", data=wo_large)
#Definitely looks like there's a pattern with number of malformations. 

wo_large$per_malf <- as.numeric(as.character(wo_large$per_malf))
#Slight pattern with percent malformation
plot(value~per_malf, data=wo_large)

boxplot(value~Malf, data=wo_large)
#Definitely seems like a pattern here

plot(value~TATO, data=wo_large)
#Maybe a pattern here, maybe driven by an outlier?

plot(value~PSRE, data=wo_large)
#Potential pattern here as well. 


# Look at the separate compartments
# Compartment 1(A)
compA_all <- subset(alldata_indiv_2010, HostID %in% comp1_2010$First_comp)
hist(compA_all$value[compA_all$value<5])
# Positive and negative values

pairs(~value+TATO+RADR+BUBO+AMCA+RACA+PSRE+PHySA+LyMN+HELI+GYRA+RADIX, data=compA_all[compA_all$value<5 ,])
plot(value~as.numeric(as.character(no_malf)), data=compA_all[compA_all$value<5 ,])
plot(value~PSRE, data=compA_all[compA_all$value<5 ,])

# Compartment 2(B)
compB_all <- subset(alldata_indiv_2010, HostID %in% comp2_2010$Second_comp)
hist(compB_all$value[compB_all$value<5])
#Basically all negative values here. 

plot(value~as.numeric(as.character(no_malf)), data=compB_all[compB_all$value<5 ,])
plot(value~PSRE, data=compB_all[compB_all$value<5 ,])

# Seems like in compartment B, there's no correlation with malformation, 
# but in CompA there is a strong one. 
# Need to figure out which variables differ between compA and compB



