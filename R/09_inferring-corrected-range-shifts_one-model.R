## script to correct range shift values for methdological variation 
## method and code developed by RB, adapted by NM
## based on script used for the analysis published in Lenoir et al. 2020 and available here: https://figshare.com/articles/dataset/BioShifts_a_global_geodatabase_of_climate-induced_species_redistribution_over_land_and_sea/7413365?file=22057815

## WHAT THIS SCRIPT DOES:
## to correct observed range shifts for differences in methodology, we fit linear mixed effect models to the raw range shift values including methodological variables as fixed effects and taxonomy as random effects
## we then extract the residuals of this model, which represent the leftover variation in range shifts after accounting for methodology

## we fit one model per combination of:
## - range shift position (leading edge, trailing edge, centroid)
## - gradient (elevation, latitude)
## - realm (terrestrial, marine)

## model selection process prevents no-convergence and singularities by systematically adapting the model structure according to the issues encountered 

library(MuMIn) #version 1.43.6
library(lme4) #version 1.1.21
library(optimx) #version 2018-7.10
library(tidyverse)

#############################
####   data preparation  ####
#############################
rs_data = read.table("data-raw/bioshifts-download/Lenoir_et_al/Analysis/Table_S1.csv",sep=";",
                     h=T,dec=".",stringsAsFactors = FALSE) 

#it's the data used in Lenoir et al. 2020; so it's a subset of Bioshift v1 available here: https://figshare.com/articles/dataset/BioShifts_a_global_geodatabase_of_climate-induced_species_redistribution_over_land_and_sea/7413365?file=22057815

#Kingdom=	Taxonomic kingdom of the species for which the range shift has been estimated (according to NCBI; September 2019) 
#Phylum=	Taxonomic phylum of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Order=	Taxonomic order of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Class=	Taxonomic class of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Family=	Taxonomic family of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Genus=	Taxonomic genus of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Species=	Validated species name of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Hemisphere=	Hemisphere of the study: North or South
#Ecosystem=	Ecosystem of the study: Marine or terrestrial (note that freshwater species have been classified as terrestrial)
#Gradient=	Spatial gradient along which the species range shift has been estimated: Elevation or Latitude
#Position=	"Range shift parameter for which the range shift has been estimated: Leading edge, Trailing edge, Centroid"
#ShiftR=	Range shift estimate standardized by the length of the study period over which the range shift has been estimated
#Unit=	Unit of the range shift estimate: m/year or km/year
#EleVeloT= elevational temperature velocity (m/year)
#LatVeloT= latitudinal temperature velocity (km/year)
#HFI = human footprint index (see Lenoir et al. 2020)
#BaselineT = baseline temperature (°C)
#LifeForm = lifeform of the species (endo, ecto, crypto, plant)
#Start = year of the beginning of the range shift monitoring
#Area = area of the polygon of the study (km2)
#Ntaxa = number of taxa monitored in the study
#PrAb=	Type of the data used to estimate the range shift: Occurrence (OCCUR) or Abundance (ABUND)
#Sampling=	"Characteristics of the time periods over which the range shift has been estimated: Continuous (CONT; yearly data or less), Irregular (IRR; multiple periods irregularly distributed), Multiple (MULT; multiple periods regularly distributed), Two (TWO; two periods)"
#Grain=	"Spatial grain of the data used to estimate the range shift: fine (FINE; spatial resolution lower than 10 km), coarse (COARSE; spatial resolution greater than 100 km), or medium (MEDIUM; intermediate situations)"
#Signif=	Whether the significance of the range shift estimate has been tested in the original publication (note that we do not report the result of such test in the database): YES or NO
#Quality=	"The quality of the approach used to estimate the range shift: LOW (no pre-processing of the raw data), BALANCED (data cleaning or resampling procedures were carried out to quantify the range shift on a balanced dataset), MODELED (range shifts were quantified species ranges modelled using  species distribution models calibrated independently for each time period); RESURVEYED (range shifts were quantified from paired designs such as permanent plots)"
#Start=	The first year of the temporal period over which the range shift has been estimated
#Reference=	The reference to the original publication
#Source=	Unique identifier of the spatial polygons delinating the study areas (provided in the geodatabase)

rs_data = rs_data[order(rs_data$n_cl),]

## transform continuous methodological variables into qualitative variables
## this has two benefits: 
## (1) all variables are one type  
## (2) qualitative is better than quantitative in case the effect is not linear
q1 = quantile(rs_data$Start, probs=c(0,0.25,0.5,0.75,1))
rs_data$StartF = cut(rs_data$Start, breaks=q1, include.lowest=T)
q1 = quantile(rs_data$Area, probs=c(0,0.25,0.5,0.75,1))
rs_data$AreaF = cut(rs_data$Area, breaks=q1, include.lowest=T)
q1 = quantile(rs_data$Ntaxa, probs=c(0,0.25,0.5,0.75,1))
rs_data$NtaxaF = cut(rs_data$Ntaxa, breaks=q1, include.lowest=T)

rs_data$Sampling = ifelse(rs_data$Sampling %in% c("IRR","MULT"),"MULT", rs_data$Sampling)

rs_data$IDn=1:nrow(rs_data)

## subset the data by gradient
rs_lat <- which(rs_data$Gradient == "Latitudinal")
rs_ele <- which(rs_data$Gradient == "Elevation")

## select variable to keep for each gradient 
chosen_varlat = c("ShiftR", "Position", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain",
                  "Quality", "Signif",
                  "Source","Class","Family","Genus","Species","Start","Area","Ntaxa","IDn")
chosen_varele = c("ShiftR", "Position", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain",
                  "Quality", "Signif",
                  "Source","Class","Family","Genus","Species","Start","Area","Ntaxa","IDn")

##########################################################
#### model fitting: random effect structure selection ####
##########################################################
## we decided to fit one model for latitudinal shifts, and a separate model for elevation shifts 

####################################
####   model fitting: latitude  ####
####################################
## Disparities between Classes and data selection
data <- rs_data[rs_lat, chosen_varlat] # create a data.frame with only variables we are interested in
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq       Shift
# 1     Actinopterygii  739  1.89190757
# 2           Amphibia  211 -0.20626687
# 4           Anthozoa   24  3.50618558
# 6          Arachnida  444  6.37105308
# 8         Asteroidea   20  1.13439103
# 9               Aves 4898  0.97659182
# 10 Bacillariophyceae   12  1.00701313
# 12          Bivalvia   71  2.16322556
# 13         Bryopsida  288 -0.46128360
# 15       Cephalopoda   19  0.09525340
# 16         Chilopoda   14  3.16326531
# 17    Chondrichthyes   83  0.45228464
# 22       Dinophyceae   23  3.49354276
# 23         Diplopoda   12  3.84920635
# 25        Echinoidea   11  3.53634994
# 26     Equisetopsida   34  0.04194176
# 28   Florideophyceae   70  1.35215098
# 30        Gastropoda  134  1.98410681
# 32           Insecta 4237  1.74941735
# 33   Lecanoromycetes   30  0.87507527
# 34        Liliopsida 1112 -0.09897762
# 35    Lycopodiopsida   18  0.08206480
# 36     Magnoliopsida 3617 -0.08973201
# 37      Malacostraca  160  1.09149536
# 38          Mammalia   14  5.99689091
# 39       Maxillopoda   55  4.84413708
# 43      Phaeophyceae   59  2.12435568
# 44         Pinopsida   98 -0.04036576
# 45        Polychaeta  143  3.71917179
# 47    Polypodiopsida  105  0.09068342
# 50          Reptilia   90  0.63731996
# 52      Sphagnopsida   20 -0.76504839
# 54       Ulvophyceae   18  1.87225530

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 16882 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF) ## colinear with Grain
table(data$StartF)
table(data$NtaxaF)



data = rs_data
data$Area = log(data$Area)

## reorder categorical factors in logical order, transform to numeric
data$Grain = as.numeric(factor(data$Grain, levels = c("FINE", "MEDIUM", "COARSE")))
data$AreaF = as.numeric(factor(data$AreaF))
data$NtaxaF = as.numeric(factor(data$NtaxaF))
data$StartF = as.numeric(factor(data$StartF))

meth <- data[,names(data) %in% chosen_varlat] #remove "n" as can capture diversity gradients
meth <- na.omit(meth)
meth[,names(meth) %in% c("Position","PrAb",
                         "Sampling",
                         "Quality", "Signif",
                         "Source","Class","Family",
                         "Genus","Species")] <- apply(meth[,names(meth) %in% c("Position","PrAb",
                                                                        "Sampling", 
                                                                        "Quality", "Signif",
                                                                        "Source","Class","Family",
                                                                        "Genus","Species")],2,
                                               function(x) as.numeric(factor(x, levels = unique(x))))
cormat <- cor(meth)
corrplot(cormat)


