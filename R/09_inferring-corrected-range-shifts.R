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
rs_data = read.table("data-raw/bioshifts-download/Lenoir_et_al/Analysis/Table_S1.csv",sep=";",h=T,dec=".",stringsAsFactors = FALSE) 
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

## subset the data by Gradient x Position x Ecosystem x Hemisphere
# Optimum - core of the distribution (O)
rsO.LT_posit <- which(rs_data$Position == "Centroid" & rs_data$Gradient == "Latitudinal" & 
                        rs_data$Ecosystem=="Terrestrial")
rsO.LM_posit <- which(rs_data$Position == "Centroid" & rs_data$Gradient == "Latitudinal" & 
                        rs_data$Ecosystem =="Marine")
rsO.ET_posit <- which(rs_data$Position == "Centroid" & rs_data$Gradient == "Elevation" & 
                        rSdata$Ecosystem =="Terrestrial")

# Margins - leading edge/trailing edge (LE & TE)
rsMrg.LT_posit <- which((rs_data$Position == "Trailing edge" | rs_data$Position == "Leading edge")
                        & rs_data$Gradient == "Latitudinal" & rs_data$Ecosystem == "Terrestrial")
rsMrg.LM_posit <- which((rs_data$Position == "Trailing edge" | rs_data$Position == "Leading edge")
                        & rs_data$Gradient == "Latitudinal" & rs_data$Ecosystem == "Marine")
rsMrg.ET_posit <- which((rs_data$Position == "Trailing edge" | rs_data$Position == "Leading edge")
                        & rs_data$Gradient=="Elevation" & rs_data$Ecosystem == "Terrestrial")

## select variable to keep for each gradient 
chosen_varlat = c("ShiftR", "Position", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif",
                  "Source","Class","Family","Genus","Species","Start","Area","Ntaxa","IDn")
chosen_varele = c("ShiftR", "Position", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif",
                "Source","Class","Family","Genus","Species","Start","Area","Ntaxa","IDn")

##########################################################
#### model fitting: random effect structure selection ####
##########################################################
## we decided to fit one model per combination of position (trailing edge, leading edge, centeroid), gradient (latitudinal, elevational) and realm (terrestrial, marine)

#function to select random effect
testSingAl=function(x,f,data,n=1:length(x)){
  #function testing for random effect structure
  #x=set of variables to test as random effect
  #f=formula of the fixed structure of the model
  #data= data base with all variables
  #n= levels of variable combination to test (start from 1 when signle variable are tested to the size of the set of variables to test in combination)
  b=1
  if(is.null(x)==F){
    for(i in n){
      v1=combn(x,i)
      if(i>1){
        v2=apply(t(v1),1,paste,sep="",collapse="+")
      }else{
        v2=t(v1)[,1]
      }
      for(j in 1:length(v2)){
        print(paste(i," : ",j,sep=""))
        f1=paste(f,"+",v2[j],sep="")
        lme1=lmer(as.formula(f1),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        a=summary(lme1)$optinfo$conv$lme4$messages
        if(is.null(a)==F){
          a=paste(a,sep="",collapse="_")
        }else{
          a=NA
        }
        res=data.frame(test=v2[j],R2m=r.squaredGLMM(lme1)[[1]],R2c=r.squaredGLMM(lme1)[[2]],aic=AIC(lme1)[[1]],aicc=AICc(lme1)[[1]],singular=isSingular(lme1)[[1]],nV=i,source=grep("Source",v2[j])[1],warning=a[[1]])
        if(b==1){
          resOK=res
        }else{
          resOK=rbind(resOK,res)
        }
        b=b+1
      }
    }
  }else{
    f1=f
    lme1=lmer(as.formula(f1),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    a=summary(lme1)$optinfo$conv$lme4$messages
    if(is.null(a)==F){
      a=paste(a,sep="",collapse="_")
    }else{
      a=NA
    }
    resOK=data.frame(test=f,R2m=r.squaredGLMM(lme1)[[1]],R2c=r.squaredGLMM(lme1)[[2]],aic=AIC(lme1)[[1]],aicc=AICc(lme1)[[1]],singular=isSingular(lme1)[[1]],warning=a[[1]])
    
  }
  
  return(resOK)
  #output:
  #test= random structure
  #R2m= marginal R2 of the model
  #R2c= conditional R2 of the model
  #aic= Akaike information criterion of the model
  #aicc= Akaike information criterion with small-sample correction
  #singular= output of the model singularity test (FALSE= no signularity; TRUE= singularity issue)
  #nV= number of variables tested as random effect
  #warning= inform for warning during the model fit (such as convergence issue)
}

#################################################################
####   model fitting: leading edge x latitude x terrestrial  ####
#################################################################
## Disparities between Classes and data selection
data <- rs_data[rsMrg.LT_posit, chosen_varlat] # create a data.frame with the data of realm, position and gradient you want to analyse
data = subset(data, Position == "Leading edge") # here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq       Shift
# 1  Actinopterygii   19  0.76429882
# 4       Arachnida  444  6.37105308
# 5            Aves  967  1.54629122
# 6       Bryopsida  288 -0.46128360
# 7       Chilopoda   14  3.16326531
# 8       Diplopoda   12  3.84920635
# 10        Insecta 3510  1.77091325
# 11     Liliopsida  132  0.13691220
# 13  Magnoliopsida  611 -0.28077401
# 14   Malacostraca   37  2.11678530
# 16      Pinopsida   31 -0.08336775
# 19       Reptilia   54  2.24699868
# 20   Sphagnopsida   20 -0.76504839

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 6139 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=1
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF"
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),]
  if(nrow(tx1)>0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
  if(isSingular(mod1)==T){#in case of no model selected considering our criteria (convergence and non signularity), a more simple random model structure is tested: (1|Family/Genus)
    f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
          f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Source)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
          }
        }
      } 
      
    }
  }
}

## look at the statistics of the selected model
R2m = r.squaredGLMM(mod1)[[1]] # marginal R2, ie part of variation explain by fixed effect (so here it's methodological factor)
R2c = r.squaredGLMM(mod1)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect (so here it's methodological factor)
aic = AIC(mod1)[[1]]
aicc = AICc(mod1)[[1]] # correction AIC useful for low sampling size

## compare the linear model and LMM coeff estimation
mod0 = lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff, xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~ShiftR, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ ShiftR, data)
summary(lm0)

hist(data$ShiftR, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$ShiftR, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$ShiftR, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$ShiftR, data$Species, mean)
x2 = tapply(data$corrected_shift, data$Species, mean) # lower mean value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2,data)
summary(lm0) # at the species level, we observe a high positive relationship 
## this means the corrected shifts did not change the pattern of range shift => good

## save the correct shift values:
corr_shift_LE_LT = data 

#################################################################
####     model fitting: leading edge x latitude x marine     ####
#################################################################
## Disparities between Classes and data selection
data <- rs_data[rsMrg.LM_posit, chosen_varlat] # create a data.frame with the data of realm, position and gradient you want to analyse
data = subset(data, Position == "Leading edge") # here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq    Shift
# 1   Actinopterygii  100 2.952814
# 2         Anthozoa   22 3.703240
# 4       Asteroidea   17 1.240541
# 5         Bivalvia   23 2.639126
# 11 Florideophyceae   25 1.368674
# 12      Gastropoda   72 1.254346
# 15    Malacostraca   22 3.409560
# 16     Maxillopoda   17 8.417573
# 19    Phaeophyceae   31 2.173821
# 20      Polychaeta   50 4.078642
# 22     Ulvophyceae   14 1.859915

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 393 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=1
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF"
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),]
  if(nrow(tx1)>0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
  if(isSingular(mod1)==T){#in case of no model selected considering our criteria (convergence and non signularity), a more simple random model structure is tested: (1|Family/Genus)
    f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
          f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Source)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
          }
        }
      } 
      
    }
  }
}

## look at the statistics of the selected model
R2m = r.squaredGLMM(mod1)[[1]] # marginal R2, ie part of variation explain by fixed effect (so here it's methodological factor)
R2c = r.squaredGLMM(mod1)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect (so here it's methodological factor)
aic = AIC(mod1)[[1]]
aicc = AICc(mod1)[[1]] # correction AIC useful for low sampling size

## compare the linear model and LMM coeff estimation
mod0 = lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff, xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~ShiftR, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ ShiftR, data)
summary(lm0)

hist(data$ShiftR, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$ShiftR, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$ShiftR, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$ShiftR, data$Species, mean)
x2 = tapply(data$corrected_shift, data$Species, mean) # lower mean value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2,data)
summary(lm0) # at the species level, we observe a high positive relationship 
## this means the corrected shifts did not change the pattern of range shift => good

## save the correct shift values:
corr_shift_LE_LM = data 

#################################################################
####  model fitting: leading edge x elevation x terrestrial  ####
#################################################################
## Disparities between Classes and data selection
data <- rs_data[rsMrg.ET_posit, chosen_varele] # create a data.frame with the data of realm, position and gradient you want to analyse
data = subset(data, Position == "Leading edge") # here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq      Shift
# 1   Actinopterygii   32  5.5664286
# 2         Amphibia   39  2.2947972
# 4             Aves  731  1.8563040
# 5        Bryopsida   24  0.1311728
# 8          Insecta  674  3.1261346
# 10 Lecanoromycetes   33 -0.1043771
# 11      Liliopsida  547  2.1884622
# 12  Lycopodiopsida   15  1.7307897
# 13   Magnoliopsida 2011  2.4199034
# 14        Mammalia  160  0.4930064
# 15       Pinopsida   62  1.1373279
# 16  Polypodiopsida   70  1.4301883
# 18        Reptilia   11  6.6115702

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 393 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=1
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF"
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),]
  if(nrow(tx1)>0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
  if(isSingular(mod1)==T){#in case of no model selected considering our criteria (convergence and non signularity), a more simple random model structure is tested: (1|Family/Genus)
    f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
          f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Source)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
          }
        }
      } 
      
    }
  }
}

## look at the statistics of the selected model
R2m = r.squaredGLMM(mod1)[[1]] # marginal R2, ie part of variation explain by fixed effect (so here it's methodological factor)
R2c = r.squaredGLMM(mod1)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect (so here it's methodological factor)
aic = AIC(mod1)[[1]]
aicc = AICc(mod1)[[1]] # correction AIC useful for low sampling size

## compare the linear model and LMM coeff estimation
mod0 = lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff, xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~ShiftR, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ ShiftR, data)
summary(lm0)

hist(data$ShiftR, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$ShiftR, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$ShiftR, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$ShiftR, data$Species, mean)
x2 = tapply(data$corrected_shift, data$Species, mean) # lower mean value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2,data)
summary(lm0) # at the species level, we observe a high positive relationship 
## this means the corrected shifts did not change the pattern of range shift => good

## save the correct shift values:
corr_shift_LE_ET = data 

#################################################################
####  model fitting: trailing edge x latitude x terrestrial  ####
#################################################################
## Disparities between Classes and data selection
data <- rs_data[rsMrg.LT_posit, chosen_varlat] # create a data.frame with the data of realm, position and gradient you want to analyse
data = subset(data, Position == "Trailing edge") # here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq      Shift
# 1           Aves  295 -0.7262751
# 4        Insecta  246  2.3350293
# 5     Liliopsida  132 -0.1617514
# 7  Magnoliopsida  584 -0.2034522
# 9      Pinopsida   18 -0.2047738
# 12      Reptilia   18 -3.5604597

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 393 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=1
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF"
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),]
  if(nrow(tx1)>0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
  if(isSingular(mod1)==T){#in case of no model selected considering our criteria (convergence and non signularity), a more simple random model structure is tested: (1|Family/Genus)
    f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
          f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Source)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
          }
        }
      } 
      
    }
  }
}

## look at the statistics of the selected model
R2m = r.squaredGLMM(mod1)[[1]] # marginal R2, ie part of variation explain by fixed effect (so here it's methodological factor)
R2c = r.squaredGLMM(mod1)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect (so here it's methodological factor)
aic = AIC(mod1)[[1]]
aicc = AICc(mod1)[[1]] # correction AIC useful for low sampling size

## compare the linear model and LMM coeff estimation
mod0 = lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff, xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~ShiftR, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ ShiftR, data)
summary(lm0)

hist(data$ShiftR, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$ShiftR, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$ShiftR, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$ShiftR, data$Species, mean)
x2 = tapply(data$corrected_shift, data$Species, mean) # lower mean value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2,data)
summary(lm0) # at the species level, we observe a high positive relationship 
## this means the corrected shifts did not change the pattern of range shift => good

## save the correct shift values:
corr_shift_TE_LT = data 

#################################################################
####    model fitting: trailing edge x latitude x marine     ####
#################################################################
## Disparities between Classes and data selection
data <- rs_data[rsMrg.ET_posit, chosen_varlat] # create a data.frame with the data of realm, position and gradient you want to analyse
data = subset(data, Position == "Trailing edge") # here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq      Shift
# 1   Actinopterygii   32  5.5664286
# 2         Amphibia   39  2.2947972
# 4             Aves  731  1.8563040
# 5        Bryopsida   24  0.1311728
# 8          Insecta  674  3.1261346
# 10 Lecanoromycetes   33 -0.1043771
# 11      Liliopsida  547  2.1884622
# 12  Lycopodiopsida   15  1.7307897
# 13   Magnoliopsida 2011  2.4199034
# 14        Mammalia  160  0.4930064
# 15       Pinopsida   62  1.1373279
# 16  Polypodiopsida   70  1.4301883
# 18        Reptilia   11  6.6115702

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 393 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=1
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF"
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),]
  if(nrow(tx1)>0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
  if(isSingular(mod1)==T){#in case of no model selected considering our criteria (convergence and non signularity), a more simple random model structure is tested: (1|Family/Genus)
    f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
          f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Source)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
          }
        }
      } 
      
    }
  }
}

## look at the statistics of the selected model
R2m = r.squaredGLMM(mod1)[[1]] # marginal R2, ie part of variation explain by fixed effect (so here it's methodological factor)
R2c = r.squaredGLMM(mod1)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect (so here it's methodological factor)
aic = AIC(mod1)[[1]]
aicc = AICc(mod1)[[1]] # correction AIC useful for low sampling size

## compare the linear model and LMM coeff estimation
mod0 = lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff, xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~ShiftR, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ ShiftR, data)
summary(lm0)

hist(data$ShiftR, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$ShiftR, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$ShiftR, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$ShiftR, data$Species, mean)
x2 = tapply(data$corrected_shift, data$Species, mean) # lower mean value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2,data)
summary(lm0) # at the species level, we observe a high positive relationship 
## this means the corrected shifts did not change the pattern of range shift => good

## save the correct shift values:
corr_shift_TE_ET = data 

#################################################################
####  model fitting: trailing edge x elevation x terrestrial ####
#################################################################
## Disparities between Classes and data selection
data <- rs_data[rsMrg.ET_posit, chosen_varele] # create a data.frame with the data of realm, position and gradient you want to analyse
data = subset(data, Position == "Trailing edge") # here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq      Shift
# 1  Actinopterygii   32  0.2311607
# 2        Amphibia   39  4.7396499
# 3            Aves  435  0.8177141
# 4         Insecta  235  3.0181040
# 5      Liliopsida  132  0.6292490
# 7   Magnoliopsida  506  1.3031581
# 8        Mammalia  162  1.2930445
# 9       Pinopsida   18  0.7943908
# 10 Polypodiopsida   14  2.2312189
# 12       Reptilia   11 -2.5206612

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 1584 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=1
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF"
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),]
  if(nrow(tx1)>0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
  if(isSingular(mod1)==T){#in case of no model selected considering our criteria (convergence and non signularity), a more simple random model structure is tested: (1|Family/Genus)
    f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
          f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Source)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
          }
        }
      } 
      
    }
  }
}

## look at the statistics of the selected model
R2m = r.squaredGLMM(mod1)[[1]] # marginal R2, ie part of variation explain by fixed effect (so here it's methodological factor)
R2c = r.squaredGLMM(mod1)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect (so here it's methodological factor)
aic = AIC(mod1)[[1]]
aicc = AICc(mod1)[[1]] # correction AIC useful for low sampling size

## compare the linear model and LMM coeff estimation
mod0 = lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[which(!is.na(mod0$coeff))], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~ShiftR, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ ShiftR, data)
summary(lm0)

hist(data$ShiftR, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$ShiftR, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$ShiftR, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$ShiftR, data$Species, mean)
x2 = tapply(data$corrected_shift, data$Species, mean) # lower mean value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2,data)
summary(lm0) # at the species level, we observe a high positive relationship 
## this means the corrected shifts did not change the pattern of range shift => good

## save the correct shift values:
corr_shift_TE_ET = data 

#################################################################
####    model fitting: optimum x latitude x terrestrial      ####
#################################################################
## Disparities between Classes and data selection
data <- rs_data[rsO.LT_posit, chosen_varlat] # create a data.frame with the data of realm, position and gradient you want to analyse
data = subset(data, Position == "Centroid") # here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq        Shift
# 1         Amphibia  201 -0.008455518
# 2             Aves 3620  0.949498673
# 6    Equisetopsida   28  0.147932035
# 7          Insecta  481  1.293053249
# 8  Lecanoromycetes   30  0.875075269
# 9       Liliopsida  845 -0.061776778
# 10  Lycopodiopsida   16  0.385623996
# 11   Magnoliopsida 2422 -0.014117112
# 13       Pinopsida   49  0.047234364
# 14  Polypodiopsida   95 -0.038819981
# 17        Reptilia   18  0.006063453

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 7805 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=1
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF"
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),]
  if(nrow(tx1)>0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
  if(isSingular(mod1)==T){#in case of no model selected considering our criteria (convergence and non signularity), a more simple random model structure is tested: (1|Family/Genus)
    f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
          f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Source)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
          }
        }
      } 
      
    }
  }
}

## look at the statistics of the selected model
R2m = r.squaredGLMM(mod1)[[1]] # marginal R2, ie part of variation explain by fixed effect (so here it's methodological factor)
R2c = r.squaredGLMM(mod1)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect (so here it's methodological factor)
aic = AIC(mod1)[[1]]
aicc = AICc(mod1)[[1]] # correction AIC useful for low sampling size

## compare the linear model and LMM coeff estimation
mod0 = lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[which(!is.na(mod0$coeff))], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~ShiftR, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ ShiftR, data)
summary(lm0)

hist(data$ShiftR, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$ShiftR, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$ShiftR, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$ShiftR, data$Species, mean)
x2 = tapply(data$corrected_shift, data$Species, mean) # lower mean value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2,data)
summary(lm0) # at the species level, we observe a high positive relationship 
## this means the corrected shifts did not change the pattern of range shift => good

## save the correct shift values:
corr_shift_O_LT = data 

#################################################################
####    model fitting: optimum x latitude x marine           ####
#################################################################
## Disparities between Classes and data selection
data <- rs_data[rsO.LM_posit, chosen_varlat] # create a data.frame with the data of realm, position and gradient you want to analyse
data = subset(data, Position == "Centroid") # here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq     Shift
# 1     Actinopterygii  582 1.6981903
# 4               Aves   16 4.0718179
# 5  Bacillariophyceae   12 1.0070131
# 6           Bivalvia   31 3.0135520
# 7        Cephalopoda   19 0.0952534
# 8     Chondrichthyes   69 0.7254884
# 10       Dinophyceae   23 3.4935428
# 14        Gastropoda   45 1.9653404
# 16      Malacostraca   83 0.1875997
# 17       Maxillopoda   26 2.2096605
# 21        Polychaeta   46 5.0198741

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 952 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=1
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF"
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),]
  if(nrow(tx1)>0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
  if(isSingular(mod1)==T){#in case of no model selected considering our criteria (convergence and non signularity), a more simple random model structure is tested: (1|Family/Genus)
    f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
          f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Source)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
          }
        }
      } 
      
    }
  }
}

## look at the statistics of the selected model
R2m = r.squaredGLMM(mod1)[[1]] # marginal R2, ie part of variation explain by fixed effect (so here it's methodological factor)
R2c = r.squaredGLMM(mod1)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect (so here it's methodological factor)
aic = AIC(mod1)[[1]]
aicc = AICc(mod1)[[1]] # correction AIC useful for low sampling size

## compare the linear model and LMM coeff estimation
mod0 = lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[which(!is.na(mod0$coeff))], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~ShiftR, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ ShiftR, data)
summary(lm0)

hist(data$ShiftR, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$ShiftR, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$ShiftR, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$ShiftR, data$Species, mean)
x2 = tapply(data$corrected_shift, data$Species, mean) # lower mean value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2,data)
summary(lm0) # at the species level, we observe a high positive relationship 
## this means the corrected shifts did not change the pattern of range shift => good

## save the correct shift values:
corr_shift_O_LM = data 

#################################################################
####    model fitting: optimum x elevation x terrestrial     ####
#################################################################
## Disparities between Classes and data selection
data <- rs_data[rsO.ET_posit, chosen_varele] # create a data.frame with the data of realm, position and gradient you want to analyse
data = subset(data, Position == "Centroid") # here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq      Shift
# 1     Actinopterygii   32 1.24505952
# 2           Amphibia  119 1.90017825
# 4               Aves 1061 0.69903773
# 5          Bryopsida   62 0.82023064
# 8            Insecta  509 1.82380285
# 9  Jungermanniopsida   25 0.29662904
# 10   Lecanoromycetes   32 0.57754630
# 11        Liliopsida  895 0.87418069
# 12    Lycopodiopsida   23 1.42287067
# 13     Magnoliopsida 4548 0.59142932
# 14          Mammalia   33 0.88353323
# 16         Pinopsida  101 0.03245349
# 17    Polypodiopsida   98 0.47605321
# 19          Reptilia   11 2.05785124

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 7549 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=1
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Source)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Source)
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }else{
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF"
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),]
  if(nrow(tx1)>0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }
  if(isSingular(mod1)==T){#in case of no model selected considering our criteria (convergence and non signularity), a more simple random model structure is tested: (1|Family/Genus)
    f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
          f1="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Source)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
          }
        }
      } 
      
    }
  }
}

## look at the statistics of the selected model
R2m = r.squaredGLMM(mod1)[[1]] # marginal R2, ie part of variation explain by fixed effect (so here it's methodological factor)
R2c = r.squaredGLMM(mod1)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect (so here it's methodological factor)
aic = AIC(mod1)[[1]]
aicc = AICc(mod1)[[1]] # correction AIC useful for low sampling size

## compare the linear model and LMM coeff estimation
mod0 = lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[which(!is.na(mod0$coeff))], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~ShiftR, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ ShiftR, data)
summary(lm0)

hist(data$ShiftR, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$ShiftR, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$ShiftR, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$ShiftR, data$Species, mean)
x2 = tapply(data$corrected_shift, data$Species, mean) # lower mean value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2,data)
summary(lm0) # at the species level, we observe a high positive relationship 
## this means the corrected shifts did not change the pattern of range shift => good

## save the correct shift values:
corr_shift_O_ET = data 

#####################################
####   combing and saving data   ####
#####################################
## combine all of the corrected range shift data into one database
corr_rs_data <- rbind(corr_shift_LE_ET, corr_shift_LE_LT) %>%
  rbind(., corr_shift_LE_LM) %>%
  rbind(., corr_shift_O_LT) %>%
  rbind(., corr_shift_O_LM) %>%
  rbind(., corr_shift_O_ET) %>%
  rbind(., corr_shift_TE_LT) %>%
  rbind(., corr_shift_TE_LM) %>%
  rbind(., corr_shift_TE_ET) 

dim(corr_rs_data)
dim(rs_data)

## make sure no data is duplicated
length(unique(corr_rs_data$IDn))
length(unique(rs_data$IDn))

## save it:
write.csv(corr_rs_data, "data-processed/corrected-range-shifts.csv", row.names = FALSE)


