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
source("R/taxonomic-harmonization/clean_taxa_functions.R")

#############################
####   data preparation  ####
#############################
## read in current Bioshifts v1
v1 = read.table("data-raw/bioshiftsv1/Shifts2018_checkedtaxo.txt",
                header = T,
                encoding="latin1") 

## read in list of all bioshifts species 
sp <- read_csv("data-raw/splist.csv")

## clean names to make them match reported names in species list
v1$reported_name = Clean_Names(gsub("_", " ", v1$Publi_name), return_gen_sps = F)

spv1 <- filter(sp, v1 == 1) %>%
  select(reported_name, scientificName) %>%
  unique()

## yay! all are there
which(!v1$reported_name %in% spv1$reported_name)

## join to add scientific name column
v1 = left_join(v1, spv1)

## clean Bioshifts 
#correct and reorganize methodological variables
v1$Data <- ifelse(v1$Data %in% c("occurence-based", "occurrence-based"), "OCCUR" , 
                      ifelse(v1$Data %in% c("abundance-based"), "ABUND", 
                             v1$Data)) 
v1$Sampling <- ifelse(v1$Sampling == "MULTIPLE(continuous)", "MULT" , 
                      ifelse(v1$Sampling == "CONTINUOUS", "CONT", 
                             v1$Sampling))
v1$Grain_size <- ifelse(v1$Grain_size %in% c("large","very_large"),"COARSE",
                        ifelse(v1$Grain_size == "small", "FINE", 
                               ifelse(v1$Grain_size == "moderate", "MEDIUM",
                                      v1$Grain_size))) 
v1$Uncertainty_Distribution <- ifelse(v1$Uncertainty_Distribution %in% c("RESAMPLING","RESAMPLING(same)"),"RESAMPLING",
                                      ifelse(v1$Uncertainty_Distribution %in% c("MODEL","MODEL+RESAMPLING(same)",
                                                                                "RESAMPLING+MODEL"),"MODEL",
                                             ifelse(v1$Uncertainty_Distribution %in% c("DETECTABILITY",
                                                                                       "RESAMPLING(sam)+DETECTABILITY"), "DETECTABILITY",
                                                    v1$Uncertainty_Distribution))) 
v1$Position <- ifelse(v1$Param == "O", "Centroid",
                      ifelse(v1$Param == "LE", "Leading edge",
                             ifelse(v1$Param == "TE", "Trailing edge", 
                                    v1$Param))) 
v1$Gradient <- ifelse(v1$Type == "ELE", "Elevation",
                      ifelse(v1$Type == "LAT", "Latitudinal",
                             v1$Type)) 
v1$Ecosystem <- ifelse(v1$ECO == "T", "Terrestrial",
                      ifelse(v1$ECO == "M", "Marine",
                             v1$ECO)) 
v1$Signif <- ifelse(v1$Sign.trend %in% c("Y", "N"), "YES", "NO")

rs_data = v1

## transform continuous methodological variables into qualitative variables
## this has two benefits: 
## (1) all variables are one type  
## (2) qualitative is better than quantitative in case the effect is not linear
q1 = quantile(rs_data$START, probs=c(0,0.25,0.5,0.75,1))
rs_data$StartF = cut(rs_data$START, breaks=q1, include.lowest=T)
q1 = quantile(rs_data$ID.area, probs=c(0,0.25,0.5,0.75,1))
rs_data$AreaF = cut(rs_data$ID.area, breaks=q1, include.lowest=T)
q1 = quantile(rs_data$N, probs=c(0,0.25,0.5,0.75,1))
rs_data$NtaxaF = cut(rs_data$N, breaks=q1, include.lowest=T)

rs_data$Sampling = ifelse(rs_data$Sampling %in% c("IRR","MULT"),"MULT", rs_data$Sampling)

rs_data$IDn=1:nrow(rs_data)

## save the data 
write.csv(rs_data, "data-processed/bioshifts-v1_before-inferring.csv", row.names = FALSE)

## subset the data by Gradient x Position x Ecosystem x Hemisphere
# Optimum - core of the distribution (O)
rsO.LT_posit <- which(rs_data$Position == "Centroid" & rs_data$Gradient == "Latitudinal" & 
                        rs_data$Ecosystem=="Terrestrial")
rsO.LM_posit <- which(rs_data$Position == "Centroid" & rs_data$Gradient == "Latitudinal" & 
                        rs_data$Ecosystem =="Marine")
rsO.ET_posit <- which(rs_data$Position == "Centroid" & rs_data$Gradient == "Elevation" & 
                        rs_data$Ecosystem =="Terrestrial")

# Margins - leading edge/trailing edge (LE & TE)
rsMrg.LT_posit <- which((rs_data$Position == "Trailing edge" | rs_data$Position == "Leading edge")
                        & rs_data$Gradient == "Latitudinal" & rs_data$Ecosystem == "Terrestrial")
rsMrg.LM_posit <- which((rs_data$Position == "Trailing edge" | rs_data$Position == "Leading edge")
                        & rs_data$Gradient == "Latitudinal" & rs_data$Ecosystem == "Marine")
rsMrg.ET_posit <- which((rs_data$Position == "Trailing edge" | rs_data$Position == "Leading edge")
                        & rs_data$Gradient=="Elevation" & rs_data$Ecosystem == "Terrestrial")

## select variable to keep for each gradient 
chosen_varlat = c("SHIFT", "Position", "NtaxaF", "StartF", "AreaF", "Data", "Sampling", "Grain_size", "Signif",
                  "Uncertainty_Distribution", "Article_ID","Class","Family","Genus","Species",
                  "START","ID.area","N","IDn", "scientificName")
chosen_varele =  c("SHIFT", "Position", "NtaxaF", "StartF", "AreaF", "Data", "Sampling", "Grain_size", "Signif",
                   "Uncertainty_Distribution", "Article_ID","Class","Family","Genus","Species",
                   "START","ID.area","N","IDn", "scientificName")

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
        res=data.frame(test=v2[j],R2m=r.squaredGLMM(lme1)[[1]],R2c=r.squaredGLMM(lme1)[[2]],aic=AIC(lme1)[[1]],aicc=AICc(lme1)[[1]],singular=isSingular(lme1)[[1]],nV=i,source=grep("Article_ID",v2[j])[1],warning=a[[1]])
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
Class$Shift <- c(tapply(data$SHIFT, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq       Shift
# 1  Actinopterygii   19  0.76429882
# 4       Arachnida  448  6.36432160
# 5            Aves  914  1.35179215
# 6       Bryopsida  289 -0.46275477
# 7       Chilopoda   14  3.16326531
# 8       Diplopoda   12  3.84920635
# 10        Insecta 3529  1.76729448
# 11     Liliopsida  127  0.14226756
# 13  Magnoliopsida  495 -0.32104058
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
dim(data) # 5988 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$Data) 
table(data$Sampling)
table(data$Grain_size)
table(data$Uncertainty_Distribution)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Article_ID)
  f1=paste(f,"+(1|Article_ID)",sep="")
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF"
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
    f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Article_ID)
          f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Article_ID)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
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
mod0 = lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[which(!is.na(mod0$coeff))], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~SHIFT, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ SHIFT, data)
summary(lm0)

hist(data$SHIFT, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$SHIFT, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$SHIFT, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$SHIFT, data$Species, mean)
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
Class$Shift <- c(tapply(data$SHIFT, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq     Shift
# 1   Actinopterygii  102 2.8741791
# 2         Anthozoa   22 3.7032404
# 3       Asteroidea   17 1.2405405
# 10 Florideophyceae   25 1.3686736
# 11      Gastropoda   69 0.4491539
# 15     Maxillopoda   17 8.4175732
# 17    Phaeophyceae   31 2.1738209
# 18      Polychaeta   22 0.5071850
# 20     Ulvophyceae   14 1.8599154

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 314 obs x 20 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$Data) 
table(data$Sampling)
table(data$Grain_size)
table(data$Uncertainty_Distribution)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Article_ID)
  f1=paste(f,"+(1|Article_ID)",sep="")
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF"
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
    f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Article_ID)
          f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Article_ID)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
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
mod0 = lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[!is.na(mod0$coeff)], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~SHIFT, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ SHIFT, data)
summary(lm0)

hist(data$SHIFT, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$SHIFT, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$SHIFT, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$SHIFT, data$Species, mean)
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
Class$Shift <- c(tapply(data$SHIFT, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq      Shift
# 1   Actinopterygii   32  5.5664286
# 2         Amphibia   40  1.7828818
# 4             Aves  732  1.8453112
# 5        Bryopsida   25  0.1170370
# 8          Insecta  677  3.1240094
# 10 Lecanoromycetes   33 -0.1043771
# 11      Liliopsida  550  2.1832822
# 12  Lycopodiopsida   15  1.7307897
# 13   Magnoliopsida 2024  2.4219791
# 14        Mammalia  160  0.4930064
# 15       Pinopsida   62  1.1373279
# 16  Polypodiopsida   72  1.4369948
# 18        Reptilia   11  6.6115702

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 4393 obs x 20 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$Data) 
table(data$Sampling)
table(data$Grain_size)
table(data$Uncertainty_Distribution)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Article_ID)
  f1=paste(f,"+(1|Article_ID)",sep="")
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF"
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
    f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Article_ID)
          f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Article_ID)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
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
mod0 = lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[!is.na(mod0$coeff)], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~SHIFT, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ SHIFT, data)
summary(lm0)

hist(data$SHIFT, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$SHIFT, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$SHIFT, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$SHIFT, data$Species, mean)
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
Class$Shift <- c(tapply(data$SHIFT, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq      Shift
# 1           Aves  298 -0.7182417
# 4        Insecta  247  2.3242374
# 5     Liliopsida  127 -0.1668823
# 7  Magnoliopsida  468 -0.2563265
# 9      Pinopsida   18 -0.2047738
# 12      Reptilia   18 -3.5604597

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 1176 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$Data) 
table(data$Sampling)
table(data$Grain_size)
table(data$Uncertainty_Distribution)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Article_ID)
  f1=paste(f,"+(1|Article_ID)",sep="")
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF"
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
    f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Article_ID)
          f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Article_ID)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
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
mod0 = lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[!is.na(mod0$coeff)], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~SHIFT, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ SHIFT, data)
summary(lm0)

hist(data$SHIFT, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$SHIFT, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$SHIFT, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$SHIFT, data$Species, mean)
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
data <- rs_data[rsMrg.LM_posit, chosen_varlat] # create a data.frame with the data of realm, position and gradient you want to analyse
data = subset(data, Position == "Trailing edge") # here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$SHIFT, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq    Shift
# 1   Actinopterygii   40 2.621704
# 8  Florideophyceae   45 1.342972
# 9       Gastropoda   13 5.985066
# 12     Maxillopoda   12 5.489802
# 14    Phaeophyceae   28 2.069591
# 15      Polychaeta   19 1.438449

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 157 obs x 20 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$Data) 
table(data$Sampling)
table(data$Grain_size)
table(data$Uncertainty_Distribution)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Article_ID)
  f1=paste(f,"+(1|Article_ID)",sep="")
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF"
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
    f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Article_ID)
          f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Article_ID)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
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
mod0 = lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff, xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~SHIFT, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ SHIFT, data)
summary(lm0)

hist(data$SHIFT, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$SHIFT, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$SHIFT, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$SHIFT, data$Species, mean)
x2 = tapply(data$corrected_shift, data$Species, mean) # lower mean value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2,data)
summary(lm0) # at the species level, we observe a high positive relationship 
## this means the corrected shifts did not change the pattern of range shift => good

## save the correct shift values:
corr_shift_TE_LM = data 

#################################################################
####  model fitting: trailing edge x elevation x terrestrial ####
#################################################################
## Disparities between Classes and data selection
data <- rs_data[rsMrg.ET_posit, chosen_varele] # create a data.frame with the data of realm, position and gradient you want to analyse
data = subset(data, Position == "Trailing edge") # here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$SHIFT, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq      Shift
# 1  Actinopterygii   32  0.2311607
# 2        Amphibia   40  6.5529768
# 3            Aves  435  0.8177141
# 4         Insecta  236  3.0536689
# 5      Liliopsida  133  0.6096402
# 7   Magnoliopsida  512  1.3057994
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
table(data$Data) 
table(data$Sampling)
table(data$Grain_size)
table(data$Uncertainty_Distribution)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Article_ID)
  f1=paste(f,"+(1|Article_ID)",sep="")
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF"
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
    f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Article_ID)
          f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Article_ID)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
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
mod0 = lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[which(!is.na(mod0$coeff))], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~SHIFT, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ SHIFT, data)
summary(lm0)

hist(data$SHIFT, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$SHIFT, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$SHIFT, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$SHIFT, data$Species, mean)
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
Class$Shift <- c(tapply(data$SHIFT, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq        Shift
# 1        Amphibia  101 -0.016486982
# 2            Aves 2972  0.969938244
# 7         Insecta  486  1.283282910
# 8      Liliopsida  180 -0.123742400
# 10  Magnoliopsida  627 -0.126242918
# 12      Pinopsida   38  0.175642912
# 13 Polypodiopsida   12  0.152322327
# 16       Reptilia   18  0.006063453

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 4424 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$Data) 
table(data$Sampling)
table(data$Grain_size)
table(data$Uncertainty_Distribution)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Article_ID)
  f1=paste(f,"+(1|Article_ID)",sep="")
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF"
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
    f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Article_ID)
          f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Article_ID)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
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
mod0 = lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[which(!is.na(mod0$coeff))], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~SHIFT, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ SHIFT, data)
summary(lm0)

hist(data$SHIFT, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$SHIFT, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$SHIFT, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$SHIFT, data$Species, mean)
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
Class$Shift <- c(tapply(data$SHIFT, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq      Shift
# 1     Actinopterygii  479  1.9901463
# 3               Aves   16  4.0718179
# 4  Bacillariophyceae   12  1.0070131
# 5           Bivalvia   17  3.0241107
# 6        Cephalopoda   16  0.2612516
# 7     Chondrichthyes   43  1.1166431
# 9        Dinophyceae   23  3.4935428
# 12        Gastropoda   41  1.4729754
# 14      Malacostraca   61 -0.3353653
# 15       Maxillopoda   26  2.2096605
# 19        Polychaeta   18  4.5912593

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 714 obs x 20 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$Data) 
table(data$Sampling)
table(data$Grain_size)
table(data$Uncertainty_Distribution)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Article_ID)
  f1=paste(f,"+(1|Article_ID)",sep="")
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF"
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
    f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Article_ID)
          f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Article_ID)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
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
mod0 = lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[which(!is.na(mod0$coeff))], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~SHIFT, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ SHIFT, data)
summary(lm0)

hist(data$SHIFT, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$SHIFT, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$SHIFT, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$SHIFT, data$Species, mean)
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
Class$Shift <- c(tapply(data$SHIFT, data$Class, mean))

## here we use a criteria of number of observation by class in order to have sufficient replicate among classes 
## select only classes with > 10 observations 
Class=Class[which(Class$Freq>10), ]
Class
# Class Freq      Shift
# 1     Actinopterygii   32 1.24505952
# 2           Amphibia  121 2.11769207
# 4               Aves 1063 0.69437395
# 5          Bryopsida   63 0.77017407
# 8            Insecta  510 1.82375618
# 9  Jungermanniopsida   27 0.27980605
# 10   Lecanoromycetes   32 0.57754630
# 11        Liliopsida  902 0.86323249
# 12    Lycopodiopsida   23 1.42287067
# 13     Magnoliopsida 4590 0.59348935
# 14          Mammalia   33 0.88353323
# 16         Pinopsida  101 0.03245349
# 17    Polypodiopsida  103 0.32935287
# 19          Reptilia   11 2.05785124

data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 7538 obs x 19 vars

## select the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$Data) 
table(data$Sampling)
table(data$Grain_size)
table(data$Uncertainty_Distribution)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|StartF)","(1|Signif)") #considered as random effect as it can drive bias in model that we want to control but it's not a truly methodological factor
f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
#setwd(rep_out)
#write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Article_ID)
  f1=paste(f,"+(1|Article_ID)",sep="")
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  #setwd(rep_out)
  #write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    f1=paste(f,"+(1|Article_ID)",sep="") #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|1/Genus) + (1|Article_ID)
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
  f="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF"
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
    f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Class)
        f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Class)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Article_ID)
          f1="SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF + (1|Article_ID)"
          mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
          if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
            print("Impossible to fit the linear mixed-effect model due to singularity issue")
            mod1=lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
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
mod0 = lm(SHIFT ~ Data + Sampling + Grain_size + Uncertainty_Distribution + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff[which(!is.na(mod0$coeff))], xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

## extract the residuals of the selected LMM (ie mod0) to represent the corrected range shift estimation
data$corrected_shift = summary(mod1)$residuals
plot(corrected_shift~SHIFT, data)
abline(a=0, b=1, col=2) ## add one to one line
lm0 = lm(corrected_shift ~ SHIFT, data)
summary(lm0)

hist(data$SHIFT, col = rgb(1,0,0,0.5))
hist(data$corrected_shift, col = rgb(0,0,1,0.5), add=T)
## by using residuals, we decrease the variation in the raw range shift observations
## we create a more homogenous variable as all the variation due to methodology has been substracted from the shift

## at the class level:
x1 = tapply(data$SHIFT, data$Class, mean)
x2 = tapply(data$corrected_shift, data$Class, mean) # lower mean value than true obs 
plot(x1 ~ x2)
lm0 = lm(x1 ~ x2,data)
summary(lm0) # it changes many things

x1 = tapply(data$SHIFT, data$Class, var)
x2 = tapply(data$corrected_shift, data$Class, var) # lower variance value than in true obs 
plot(x1 ~ x2)
lm0=lm(x1 ~ x2, data)
summary(lm0) # but high positive correlation meaning that relative variation is conserved

## at the species level:
x1 = tapply(data$SHIFT, data$Species, mean)
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
write.csv(corr_rs_data, "data-processed/corrected-range-shifts_new-bsv1.csv", row.names = FALSE)


