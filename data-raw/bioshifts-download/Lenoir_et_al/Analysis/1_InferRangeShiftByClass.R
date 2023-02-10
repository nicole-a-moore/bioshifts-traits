#work conducted on R version 3.6.1 (2019-07-05); Platform: x86_64-pc-linux-gnu (64-bit); Running under: Ubuntu 16.04.5 LTS
#authors: Bertrand Romain, Lise Comte, Bourgaud Luana and Jonathan Lenoir
#AIM= Inferring range shift at the taxonimic class level

#required packages
#library(nlme) #version 3.1.140
library(MuMIn) #version 1.43.6
library(lme4) #version 1.1.21
library(optimx) #version 2018-7.10

#repository/file location
rep_data="....../Script" #unzip Scrip.zip available at figshare
rep_out="......." #inform the location where all results have to be saved


#############################
#### DATASET preparation ####
#############################
setwd(rep_data)
rSdata = read.table("Table_S1.csv",sep=";",h=T,dec=".",stringsAsFactors = FALSE) #exported from Table_S1.xlsx with separator=";" and decimal="."
rSdata = rSdata[order(rSdata$n_cl),]

#transforming continuous method variables in qualitative variables
q1=quantile(rSdata$Start,probs=c(0,0.25,0.5,0.75,1))
rSdata$StartF=cut(rSdata$Start,breaks=q1,include.lowest=T)
q1=quantile(rSdata$Area,probs=c(0,0.25,0.5,0.75,1))
rSdata$AreaF=cut(rSdata$Area,breaks=q1,include.lowest=T)
q1=quantile(rSdata$Ntaxa,probs=c(0,0.25,0.5,0.75,1))
rSdata$NtaxaF=cut(rSdata$Ntaxa,breaks=q1,include.lowest=T)

rSdata$Sampling = ifelse(rSdata$Sampling %in% c("IRR","MULT"),"MULT", rSdata$Sampling)

#data subsetting: Gradient x Position x Ecosystem x Hemisphere
#Core of the distribution (O)
rSO.LTS_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Terrestrial" & rSdata$Hemisphere=="South")
rSO.LTN_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Terrestrial" & rSdata$Hemisphere=="North")
rSO.LMS_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Marine" & rSdata$Hemisphere=="South")
rSO.LMN_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Marine" & rSdata$Hemisphere=="North")
rSO.ETS_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Elevation"& rSdata$Ecosystem=="Terrestrial" & rSdata$Hemisphere=="South")
rSO.ETN_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Elevation"& rSdata$Ecosystem=="Terrestrial" & rSdata$Hemisphere=="North")

#Margins (LE & TE)
rSMrg.LTS_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Terrestrial" & rSdata$Hemisphere=="South")
rSMrg.LTN_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Terrestrial" & rSdata$Hemisphere=="North")
rSMrg.LMS_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Marine" & rSdata$Hemisphere=="South")
rSMrg.LMN_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Marine" & rSdata$Hemisphere=="North")
rSMrg.ETS_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Elevation"& rSdata$Ecosystem=="Terrestrial" & rSdata$Hemisphere=="South")
rSMrg.ETN_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Elevation"& rSdata$Ecosystem=="Terrestrial" & rSdata$Hemisphere=="North")

# variable pre-selection
chosen_varlat=c("ShiftR", "Position", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif","Source","Class","Family","Genus","Species")
chosen_varele=c("ShiftR", "Position", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif","Source","Class","Family","Genus","Species")

#####################################################
#### MODELING: random effect structure selection ####
#####################################################
# OPTIMUM x LATITUDINAL RANGE SHITS x TERRESTRIAL ECOSYSTEM x NORTH HEMISPHERE
# oLTN hereafter
# Disparities between Classes
Class <- as.data.frame(table(rSdata$Class[rSO.LTN_posit]))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(rSdata$ShiftR[rSO.LTN_posit], rSdata$Class[rSO.LTN_posit], mean))
Class=Class[which(Class$Freq>30), ]
Class
#          Class Freq       Shift (km/yr)
#           Aves 3155  0.98827642
#        Insecta  481  1.29305325
#     Liliopsida  845 -0.06177678
#  Magnoliopsida 2422 -0.01411711
#      Pinopsida   49  0.04723436
# Polypodiopsida   95 -0.03881998

# Selecting observation and variables to analyse range shifts
# Criteria: Class > 30 obs
data <- rSdata[rSO.LTN_posit, chosen_varlat]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 7047 obs x 15 vars

t1=data.frame(table(data$Class))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$Class=relevel(data$Class,ref=ref1) #Intercept value of the model will be the range shift estimation of the class having the greatest number of observation

class_ref=data.frame(model="oLTN",ref=levels(data$Class)[1])

testSingAl=function(x,f,data,n=1:length(x)){
  #function testing for random effect structure
    #x=set of variables to test as random effect
    #f=formula of the fixed structure of the model
    #data= data base with all variables
    #n= levels of variable combination to test (start from 1 when signle variable are tested to the size of the set of variables to test in combination)
  b=1
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

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb)
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif)
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)


#selection of the best random effect structure
x=c("(1|AreaF)","(1|StartF)","(1|NtaxaF)","(1|PrAb)","(1|Sampling)","(1|Grain)","(1|Quality)","(1|Signif)")
f="ShiftR ~ Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_oLTN.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:lenth(x))
  setwd(rep_out)
  write.table(tx,"sing_oLTN.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_oLTN.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Class, data)
        }
      }
    }
  }
}

final.oLTN=mod1

# OPTIMUM x LATITUDINAL RANGE SHITS x TERRESTRIAL ECOSYSTEM x SOUTH HEMISPHERE
# oLTS hereafter
# Disparities between Classes
Class <- as.data.frame(table(rSdata$Class[rSO.LTS_posit]))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(rSdata$ShiftR[rSO.LTS_posit], rSdata$Class[rSO.LTS_posit], mean))
Class=Class[which(Class$Freq>30), ]
Class
#    Class Freq        Shift (km/yr)
# Amphibia  201 -0.008455518
#     Aves  465  0.686393765

# Selecting observation and variables to analyse range shifts
# Criteria: Class > 30 obs
data <- rSdata[rSO.LTS_posit, chosen_varlat]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 666 obs x 15 vars

t1=data.frame(table(data$Class))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$Class=relevel(data$Class,ref=ref1) #Intercept value of the model will be the range shift estimation of the class having the greatest number of observation

class_ref=rbind(class_ref,data.frame(model="oLTS",ref=levels(data$Class)[1]))

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb) #not selected
table(data$Sampling)
table(data$Grain) #not selected
table(data$Signif) #not selected
table(data$Quality)
table(data$AreaF) #not selected
table(data$StartF) #not selected
table(data$NtaxaF) #not selected


#selection of the best random effect structure
x=c("(1|Sampling)","(1|Quality)")
f="ShiftR ~ Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_oLTS.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_oLTS.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_oLTS.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Class, data)
        }
      }
    }
  }
}

final.oLTS=mod1


# OPTIMUM x LATITUDINAL RANGE SHITS x Marine ECOSYSTEM x NORTH HEMISPHERE
# oLMN hereafter
# Disparities between Classes
Class <- as.data.frame(table(rSdata$Class[rSO.LMN_posit]))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(rSdata$ShiftR[rSO.LMN_posit], rSdata$Class[rSO.LMN_posit], mean))
Class=Class[which(Class$Freq>30), ]
Class
#    Class       Freq   Shift (km/yr)
# Actinopterygii  511   1.3277541
#       Bivalvia   31   3.0135520
# Chondrichthyes   51   0.8507834
#     Gastropoda   45   1.9653404
#   Malacostraca   83   0.1875997
#     Polychaeta   46   5.0198741

# Selecting observation and variables to analyse range shifts
# Criteria: Class > 30 obs
data <- rSdata[rSO.LMN_posit, chosen_varlat]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 767 obs x 15 vars

t1=data.frame(table(data$Class))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$Class=relevel(data$Class,ref=ref1) #Intercept value of the model will be the range shift estimation of the class having the greatest number of observation

class_ref=rbind(class_ref,data.frame(model="oLMN",ref=levels(data$Class)[1]))

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$PrAb)
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif)
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure
x=c("(1|PrAb)","(1|Sampling)","(1|Grain)","(1|Quality)","(1|Signif)","(1|AreaF)","(1|StartF)","(1|NtaxaF)")
f="ShiftR ~ Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_oLMN.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_oLMN.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_oLMN.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Class, data)
        }
      }
    }
  }
}

final.oLMN=mod1


# OPTIMUM x LATITUDINAL RANGE SHITS x Marine ECOSYSTEM x SOUTH HEMISPHERE
# oLMS hereafter
# Disparities between Classes
Class <- as.data.frame(table(rSdata$Class[rSO.LMS_posit]))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(rSdata$ShiftR[rSO.LMS_posit], rSdata$Class[rSO.LMS_posit], mean))
Class=Class[which(Class$Freq>30), ]
Class
#          Class  Freq  Shift (km/yr)
# Actinopterygii    71  4.364287

# Selecting observation and variables to analyse range shifts
# Criteria: Class > 30 obs
data <- rSdata[rSO.LMS_posit, chosen_varlat]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 71 obs x 15 vars

t1=data.frame(table(data$Class))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$Class=relevel(data$Class,ref=ref1) #Intercept value of the model will be the range shift estimation of the class having the greatest number of observation

class_ref=rbind(class_ref,data.frame(model="oLMS",ref=levels(data$Class)[1]))

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)table(data$PrAb)
table(data$PrAb) #not selected
table(data$Sampling) #not selected
table(data$Grain) #not selected
table(data$Quality) #not selected
table(data$Signif) #not selected
table(data$AreaF) #not selected
table(data$StartF) #not selected
table(data$NtaxaF) #not selected
table(data$Source)

#selection of the best random effect structure
x=c("(1|Source)")
f="ShiftR ~ 1 + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1)
setwd(rep_out)
write.table(tx,"sing_oLMS.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ 1 + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_oLMS.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ 1"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_oLMS.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ 1 + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ 1 + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ 1 + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~1, data)
        }
      }
    }
  }
}

final.oLMS=mod1

# OPTIMUM x ELEVATIONAL RANGE SHITS x TERRESTRIAL ECOSYSTEM x NORTH HEMISPHERE
# oETN hereafter
# Disparities between Classes
Class <- as.data.frame(table(rSdata$Class[rSO.ETN_posit]))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(rSdata$ShiftR[rSO.ETN_posit], rSdata$Class[rSO.ETN_posit], mean))
Class=Class[which(Class$Freq>30), ]
Class
#      Class      Freq      Shift (m/yr)
#  Actinopterygii   32 1.24505952
#            Aves 1061 0.69903773
#       Bryopsida   62 0.82023064
#         Insecta  507 1.81817683
# Lecanoromycetes   32 0.57754630
#      Liliopsida  895 0.87418069
#   Magnoliopsida 4513 0.57693785
#        Mammalia   33 0.88353323
#       Pinopsida  101 0.03245349
#  Polypodiopsida   97 0.42529087

# Selecting observation and variables to analyse range shifts
# Criteria: Class > 30 obs
data <- rSdata[rSO.ETN_posit, chosen_varele]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 7333 obs x 15 vars

t1=data.frame(table(data$Class))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$Class=relevel(data$Class,ref=ref1)

class_ref=rbind(class_ref,data.frame(model="oETN",ref=levels(data$Class)[1])) #Intercept value of the model will be the range shift estimation of the class having the greatest number of observation

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)table(data$PrAb)
table(data$PrAb) 
table(data$Sampling) 
table(data$Grain) 
table(data$Quality) 
table(data$Signif) 
table(data$AreaF) 
table(data$StartF) 
table(data$NtaxaF) 


#selection of the best random effect structure
x=c("(1|PrAb)","(1|Sampling)","(1|Grain)","(1|Quality)","(1|Signif)","(1|AreaF)","(1|StartF)","(1|NtaxaF)")
f="ShiftR ~ Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,1:length(x))
setwd(rep_out)
write.table(tx,"sing_oETN.csv",sep=";",dec=".",row=F)
tx$Issue=grepl("failed to converge",as.character(tx$warning))
tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
tx1=tx1[order(tx1$aicc,decreasing=F),] 

if(nrow(tx1)==0){ #in case of no model selected considering our criteria (convergence and non signularity), a new random model structure is tested: (1|Family/Genus) + (1|Source)
  f1=paste(f,"+(1|Source)",sep="")
  mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
}else{
  z=2
  iS=T
  while(iS==T & z<=nrow(tx1)){
    f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mod1)[[1]]
    z=z+1
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted substituing (1|Family/Genus) by (1|Genus)
  f="ShiftR ~ Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_oETN.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_oETN.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Class, data)
        }
      }
    }
  }
}

final.oETN=mod1

# OPTIMUM x ELEVATIONAL RANGE SHITS x TERRESTRIAL ECOSYSTEM x SOUTH HEMISPHERE
# oETS hereafter
# Disparities between Classes
Class <- as.data.frame(table(rSdata$Class[rSO.ETS_posit]))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(rSdata$ShiftR[rSO.ETS_posit], rSdata$Class[rSO.ETS_posit], mean))
Class=Class[which(Class$Freq>30), ]
Class
#         Class Freq    Shift (m/yr)
#      Amphibia  119 1.900178
# Magnoliopsida   35 2.460000

# Selecting observation and variables to analyse range shifts
# Criteria: Class > 30 obs
data <- rSdata[rSO.ETS_posit, chosen_varele]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 154 obs x 15 vars

t1=data.frame(table(data$Class))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$Class=relevel(data$Class,ref=ref1) #Intercept value of the model will be the range shift estimation of the class having the greatest number of observation

class_ref=rbind(class_ref,data.frame(model="oETS",ref=levels(data$Class)[1]))

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)table(data$PrAb)
table(data$PrAb) 
table(data$Sampling) #not selected
table(data$Grain) #not selected
table(data$Signif) #not selected 
table(data$Quality) 
table(data$AreaF) #not selected
table(data$StartF) #not selected
table(data$NtaxaF) #not selected

#selection of the best random effect structure
x=c("(1|PrAb)","(1|Quality)")
f="ShiftR ~ Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_oETS.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_oETS.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_oETS.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Class, data)
        }
      }
    }
  }
}

final.oETS=mod1

# MARGIN x LATITUDINAL RANGE SHITS x TERRESTRIAL ECOSYSTEM x NORTH HEMISPHERE
# MrgLTN hereafter
# Disparities between Classes
data <- rSdata[rSMrg.LTN_posit, chosen_varlat]
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))
Class=Class[which(Class$Freq>30), ]
Class
#          Class Freq       Shift
#      Arachnida  444  6.37105308
#           Aves 1262  1.01506535
#      Bryopsida  288 -0.46128360
#        Insecta 3748  1.78748756
#     Liliopsida  264 -0.01241958
#  Magnoliopsida 1195 -0.24298664
#   Malacostraca   37  2.11678530
#      Pinopsida   49 -0.12796589
#       Reptilia   72  0.79513408


# Selecting observation and variables to analyse range shifts
# Criteria: Class > 30 obs
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 7359 obs x 15 vars

t1=data.frame(table(data$Class))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$Class=relevel(data$Class,ref=ref1) #Intercept value of the model will be the range shift estimation of the class having the greatest number of observation

class_ref=rbind(class_ref,data.frame(model="MrgLTN",ref=levels(data$Class)[1]))

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)table(data$PrAb)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|PrAb)","(1|Sampling)", "(1|Grain)","(1|Quality)","(1|Signif)","(1|AreaF)","(1|StartF)","(1|NtaxaF)")
f="ShiftR ~ Position*Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Position*Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Position*Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLTN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Position*Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Position*Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Position*Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Position*Class, data)
        }
      }
    }
  }
}

final.MrgLTN_int=mod1

#selection of the best random effect structure (model with additive effect between Class and Position)
x=c("(1|PrAb)","(1|Sampling)", "(1|Grain)","(1|Quality)","(1|Signif)","(1|AreaF)","(1|StartF)","(1|NtaxaF)")
f="ShiftR ~ Position + Class + (1|Family/Genus)"
a1=Sys.time()
tx=testSingAl(x,f,data,n=1:length(x))
Sys.time()-a1
setwd(rep_out)
write.table(tx,"sing_MrgLTN_add.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Position + Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLTN_add.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  tx$Issue=grepl("failed to converge",as.character(tx$warning))
  tx1=subset(tx,singular==F & is.na(source)==T & Issue==F) #selection of random effect structure for which the model converged and had not any singularity issue
  tx1=tx1[order(tx1$aicc,decreasing=F),] 
  if(nrow(tx1)==0){
    z=1
    iS=T
    while(iS==T & z<=nrow(tx1)){
      f1=paste(f,"+",tx1$test[z],sep="") #best model selection based on AICc
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      iS=isSingular(mod1)[[1]]
      z=z+1
    }
  }else{
    f1=paste(f,"+",tx1$test[1],sep="") #best model selection based on AICc
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  }
}

if(isSingular(mod1)==T){ #in case of the selected model (with REML=T) is singular, a new random model structure selection is conducted without any phylogenetic effect
  f="ShiftR ~ Position + Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLTN_add.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  if(isSingular(mod1)==T){ #in case of no model selected considering our criteria (convergence and non signularity), a more simple random model structure is tested: (1|Family/Genus)
    f1="ShiftR ~ Position + Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Position + Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Position + Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Position + Class, data)
        }
      }
    }
  }
}

final.MrgLTN_add=mod1

#selection of best model: Class*Param vs Class+Param
mInt=lmer(formula(final.MrgLTN_int),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
mAdd=lmer(formula(final.MrgLTN_add),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
bmod=AICc(mInt,mAdd)

if(bmod$AICc[1]-bmod$AICc[2]>2){
  print("Class+Param is the best model")
  final.MrgLTN=final.MrgLTN_add
}else{
  print("Class*Param is the best model") #best model
  final.MrgLTN=final.MrgLTN_int
}

# MARGIN x LATITUDINAL RANGE SHITS x TERRESTRIAL ECOSYSTEM x SOUTH HEMISPHERE
data <- rSdata[rSMrg.LTS_posit, chosen_varlat] # 8 observations
#unsufficient number of observation: 8 obs
#no model fitted

# MARGIN x LATITUDINAL RANGE SHITS x Marine ECOSYSTEM x SOUTH HEMISPHERE
# MrgLMN hereafter
# Disparities between Classes
data <- rSdata[rSMrg.LMN_posit, chosen_varlat]

Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))
Class=Class[which(Class$Freq>30), ]
Class
#           Class Freq    Shift
#  Actinopterygii   89 3.961302
#        Bivalvia   37 1.543026
#      Gastropoda   37 4.817614
#    Malacostraca   37 2.251287
#      Polychaeta   95 3.171517

# Selecting observation and variables to analyse range shifts
# Criteria: Class > 30 obs
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 295 obs x 15 variables

t1=data.frame(table(data$Class))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$Class=relevel(data$Class,ref=ref1) #Intercept value of the model will be the range shift estimation of the class having the greatest number of observation

class_ref=rbind(class_ref,data.frame(model="MrgLMN",ref=levels(data$Class)[1]))

#selection of the set of qualitative variables to test as random effect (more than one modality with 10 obsersation at least and correlation among variables)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF) #not selected (only one modality)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|PrAb)","(1|Sampling)", "(1|Grain)","(1|Quality)","(1|Signif)","(1|AreaF)","(1|StartF)")
f="ShiftR ~ Position*Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_MrgLMN_int.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Position*Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLMN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Position*Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLMN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Position*Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Position*Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Position*Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Position*Class, data)
        }
      }
    }
  }
}

final.MrgLMN_int=mod1

#selection of the best random effect structure (model with additive effect between Class and Position)
x=c("(1|PrAb)","(1|Sampling)", "(1|Grain)","(1|Quality)","(1|Signif)","(1|AreaF)","(1|StartF)")
f="ShiftR ~ Position + Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_MrgLMN_add.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Position + Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLMN_add.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Position + Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLMN_add.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Position + Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Position + Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Position + Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Position + Class, data)
        }
      }
    }
  }
}

final.MrgLMN_add=mod1

#selection of best model: Class*Param vs Class+Param
mInt=lmer(formula(final.MrgLMN_int),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
mAdd=lmer(formula(final.MrgLMN_add),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
bmod=AICc(mInt,mAdd)

if(bmod$AICc[1]-bmod$AICc[2]>2){
  print("Class+Param is the best model")
  final.MrgLMN=final.MrgLMN_add
}else{
  print("Class*Param is the best model") #best model
  final.MrgLMN=final.MrgLMN_int
}


# MARGIN x LATITUDINAL RANGE SHITS x MARINE ECOSYSTEM x SOUTH HEMISPHERE
# MrgLMS hereafter
# Disparities between Classes
data <- rSdata[rSMrg.LMS_posit, chosen_varlat]

Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))
Class=Class[which(Class$Freq>30), ]
Class
#             Class Freq    Shift
#   Florideophyceae   35 1.630896

# Selecting observation and variables to analyse range shifts
# Criteria: Class > 30 obs
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 194 obs x 15 variables

t1=data.frame(table(data$Class))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$Class=relevel(data$Class,ref=ref1) #Intercept value of the model will be the range shift estimation of the class having the greatest number of observation

class_ref=rbind(class_ref,data.frame(model="MrgLMS",ref=levels(data$Class)[1]))

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)table(data$PrAb)
table(data$PrAb) #not selected
table(data$Sampling) #not selected
table(data$Grain) #not selected
table(data$Signif) #not selected
table(data$Quality)
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF) #not selected

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|Quality)","(1|AreaF)","(1|StartF)")
f="ShiftR ~ Position*Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_MrgLMS_int.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Position*Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLMS_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Position*Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLMS_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Position*Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Position*Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Position*Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Position*Class, data)
        }
      }
    }
  }
}

final.MrgLMS_int=mod1

#selection of the best random effect structure (model with additive effect between Class and Position)
x=c("(1|Quality)","(1|AreaF)","(1|StartF)")
f="ShiftR ~ Position + Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_MrgLMS_add.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Position + Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLMS_add.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Position + Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgLMS_add.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Position + Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Position + Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Position + Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Position + Class, data)
        }
      }
    }
  }
}

final.MrgLMS_add=mod1

#selection of best model: Class*Param vs Class+Param
mInt=lmer(formula(final.MrgLMS_int),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
mAdd=lmer(formula(final.MrgLMS_add),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
bmod=AICc(mInt,mAdd)

if(bmod$AICc[1]-bmod$AICc[2]>2){
  print("Class+Param is the best model")
  final.MrgLMS=final.MrgLMS_add
}else{
  print("Class*Param is the best model") #best model
  final.MrgLMS=final.MrgLMS_int
}


# MARGIN x ELEVATIONAL RANGE SHITS x TERRESTRIAL ECOSYSTEM x NORTH HEMISPHERE
# MrgETN hereafter
# Disparities between Classes
data <- rSdata[rSMrg.ETN_posit, chosen_varele]

Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))
Class=Class[which(Class$Freq>30), ]
Class
#             Class Freq      Shift
#    Actinopterygii   64  2.8987946
#          Amphibia   40 -1.1527778
#              Aves  945  1.2695998
#           Insecta  891  3.0808897
#   Lecanoromycetes   33 -0.1043771
#        Liliopsida  675  1.8809083
#     Magnoliopsida 2418  2.1410703
#          Mammalia  322  0.8955100
#         Pinopsida   80  1.0601670
#    Polypodiopsida   79  1.6188206

# Selecting observation and variables to analyse range shifts
# Criteria: Class > 30 obs
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 5547 obs x 15 variables

t1=data.frame(table(data$Class))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$Class=relevel(data$Class,ref=ref1) #Intercept value of the model will be the range shift estimation of the class having the greatest number of observation

class_ref=rbind(class_ref,data.frame(model="MrgETN",ref=levels(data$Class)[1]))

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)table(data$PrAb)
table(data$PrAb)
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif)
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|PrAb)","(1|Sampling)","(1|Grain)","(1|Quality)","(1|Signif)","(1|AreaF)","(1|StartF)","(1|AreaF)")
f="ShiftR ~ Position*Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_MrgETN_int.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Position*Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgETN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Position*Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgETN_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Position*Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Position*Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Position*Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Position*Class, data)
        }
      }
    }
  }
}

final.MrgETN_int=mod1

#selection of the best random effect structure (model with additive effect between Class and Position)
x=c("(1|PrAb)","(1|Sampling)","(1|Grain)","(1|Quality)","(1|Signif)","(1|AreaF)","(1|StartF)","(1|AreaF)")
f="ShiftR ~ Position + Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_MrgETN_add.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Position + Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgETN_add.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Position + Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgETN_add.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Position + Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Position + Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Position + Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Position + Class, data)
        }
      }
    }
  }
}

final.MrgETN_add=mod1

#selection of best model: Class*Param vs Class+Param
mInt=lmer(formula(final.MrgETN_int),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
mAdd=lmer(formula(final.MrgETN_add),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
bmod=AICc(mInt,mAdd)

if(bmod$AICc[1]-bmod$AICc[2]>2){
  print("Class+Param is the best model")
  final.MrgETN=final.MrgETN_add
}else{
  print("Class*Param is the best model") #best model
  final.MrgETN=final.MrgETN_int
}


# MARGIN x ELEVATIONAL RANGE SHITS x TERRESTRIAL ECOSYSTEM x SOUTH HEMISPHERE
# MrgETS hereafter
# Disparities between Classes
data <- rSdata[rSMrg.ETS_posit, chosen_varele]

Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))
Class=Class[which(Class$Freq>30), ]
Class
#         Class Freq    Shift
#      Amphibia   38 8.433014
#          Aves  221 2.320779
# Magnoliopsida   99 3.522383

# Selecting observation and variables to analyse range shifts
# Criteria: Class > 30 obs
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 358 obs x 15 variables

t1=data.frame(table(data$Class))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$Class=relevel(data$Class,ref=ref1) #Intercept value of the model will be the range shift estimation of the class having the greatest number of observation

class_ref=rbind(class_ref,data.frame(model="MrgETS",ref=levels(data$Class)[1]))

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)table(data$PrAb)
table(data$PrAb) #not selected
table(data$Sampling) #not selected
table(data$Grain) #not selected
table(data$Signif) #not selected
table(data$Quality)
table(data$AreaF) #not selected
table(data$StartF)
table(data$NtaxaF)

#selection of the best random effect structure (model with interaction between Class and Position)
x=c("(1|Quality)","(1|StartF)","(1|NtaxaF)")
f="ShiftR ~ Position*Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_MrgETS_int.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Position*Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgETS_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Position*Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgETS_int.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Position*Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Position*Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Position*Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Position*Class, data)
        }
      }
    }
  }
}

final.MrgETS_int=mod1

#selection of the best random effect structure (model with additive effect between Class and Position)
x=c("(1|Quality)","(1|StartF)","(1|NtaxaF)")
f="ShiftR ~ Position + Class + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x))
setwd(rep_out)
write.table(tx,"sing_MrgETS_add.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ Position + Class + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgETS_add.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ Position + Class"
  tx=testSingAl(x,f,data,n=1:length(x))
  setwd(rep_out)
  write.table(tx,"sing_MrgETS_add.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ Position + Class + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ Position + Class + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ Position + Class + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR~Position + Class, data)
        }
      }
    }
  }
}

final.MrgETS_add=mod1

#selection of best model: Class*Param vs Class+Param
mInt=lm(formula(final.MrgETS_int),data,na.action="na.fail")
mAdd=lm(formula(final.MrgETS_add),data,na.action="na.fail")
bmod=AICc(mInt,mAdd)

if(bmod$AICc[1]-bmod$AICc[2]>2){
  print("Class+Param is the best model")
  final.MrgETS=final.MrgETS_add
}else{
  print("Class*Param is the best model") #best model
  final.MrgETS=final.MrgETS_int
}

#write.table(class_ref,"class_ref.csv",sep=";",dec=".",row=F)
save.image("analysis_Class.RData")

#############################
#### MODELING: bootstrap ####
#############################
rm(list=ls())
gc(reset=T)

#required packages
#library(nlme) #version 3.1.140
library(MuMIn) #version 1.43.6
library(lme4) #version 1.1.21
library(optimx) #version 2018-7.10
library(data.table) #version 1.12.2

#repository/file location
rep_data="....../Script" #unzip Scrip.zip available at figshare
rep_out="......." #inform the location where all results have to be saved

setwd(rep_out)
load("analysis_Class.RData")

dir.create("model_cl6000")
setwd("model_cl6000")

#Function
  # exporting coefficents
fixEff<- function(.) {
  summary(.)$coeff
}

  #bootstrap model function adapted from the "bootMer" function (merTools package).
    #see help ?bootMer in R for a detailed description of the function
bootMer2=function (x, FUN, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA, 
                   type = c("parametric", "semiparametric"), verbose = FALSE, 
                   .progress = "none", PBargs = list(), parallel = c("no", "multicore", 
                                                                     "snow"), ncpus = getOption("boot.ncpus", 1L), cl = NULL) 
{
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  
  if (!is.null(seed)) 
    set.seed(seed)
  else if (!exists(".Random.seed", envir = .GlobalEnv)) 
    runif(1)
  mc <- match.call()
  t0 <- FUN(x)
  mle <- list(beta = getME(x, "beta"), theta = getME(x, "theta"))
  if (isLMM(x)) 
    mle <- c(mle, list(sigma = sigma(x)))
  if (type == "parametric") {
    argList <- list(x, nsim = nsim, na.action = na.exclude)
    if (!missing(re.form)) {
      argList <- c(argList, list(re.form = re.form)) #equivalent to use.u=F when re.form=NA
    }
    else {
      argList <- c(argList, list(use.u = use.u))
    }
    ss <- do.call(simulate, argList) 
  }
  else {
    if (!missing(re.form)) 
      stop(paste(sQuote("re.form")), "cannot be used with semiparametric bootstrapping")
    if (use.u) {
      if (isGLMM(x)) 
        warning("semiparametric bootstrapping is questionable for GLMMs")
      ss <- replicate(nsim, fitted(x) + sample(residuals(x, 
                                                         "response"), replace = TRUE), simplify = FALSE)
    }
    else {
      stop("semiparametric bootstrapping with use.u=FALSE not yet implemented")
    }
  }
  ffun <- local({
    FUN
    refit
    x
    ss
    length.t0 <- length(t0)
    f1 <- function(i) FUN(refit(x, ss[[i]])) #refit
    function(i) {
      ret <- f1(i)
      ret
    }
  })
  simvec <- seq_len(nsim)
  res <- parallel::mclapply(simvec, ffun, mc.cores = ncpus)
  res
}

  #Computing R2 and coefficient values from bootrap model
statF=function(.){
  a=summary(.)$optinfo$conv$lme4$messages
  b=isSingular(.)[[1]]
  r2=r.squaredGLMM(.)
  if(is.null(a)==F){
    a=paste(a,sep="",collapse="_")
  }else{
    a=NA
  }
  return(list(r2=data.frame(r2m=r2[1,1],r2c=r2[1,2],warnM=as.character(a[[1]]),sing=b),coeff=data.frame(summary(.)$coeff,var=row.names(summary(.)$coeff)),randEff=as.data.frame(ranef(.,condVar=T))))
}


#Computation
r2=NULL
names=c("oLTN","oLTS","oLMN","oLMS","oETN","oETS","MrgLTN","MrgLMN","MrgLMS","MrgETN")
models=list(final.oLTN,final.oLTS,final.oLMN,final.oLMS,final.oETN,final.oETS,final.MrgLTN,final.MrgLMN,final.MrgLMS,final.MrgETN)

nB=6000 #set the bootstrap number
s1=10 #the seed
ncpu=9 #the number of cpu used for parallel computation

for (m in 1:length(models)){
  gc(reset=T)
  print(models[[m]])
  final=models[[m]]
  
  boot2.time <- system.time(
    boot2 <- bootMer2(final, statF, nsim=nB, use.u=F, type="parametric",ncpus=ncpu,seed=s1,parallel="multicore",.progress = "none")
  )
  
  #exporting results
  ess=simplify2array(boot2)
  r2=rbindlist(ess[1,])
  coeff=rbindlist(ess[2,])
  rand=rbindlist(ess[3,])
  
  write.table(r2,paste("r2_",names[m],".csv",sep=""),sep=";",dec=".",row=F)
  #output:
    #r2m= marginal R2 of every bootstrap model (each line is a different model)
    #r2c= conditional R2 of every bootstrap model
    #warnM= inform for warning during the model fit (such as convergence issue)
    #sing= output of the model singularity test (FALSE= no signularity; TRUE= singularity issue)

  write.table(coeff,paste("coeff_",names[m],".csv",sep=""),sep=";",dec=".",row=F)
  #output is the result of model$coefficient for everey bootstrap model:
    #Estimate= coefficient value computed from every boostrap model
    #Std..Error= Standard error computed from every boostrap model
    #t.value= t-value computed from every boostrap model
    #var= names of the term
}

set.seed(s1)
#special case of MrgETS
for (m in 1:nB){
  gc(reset=T)
  data1=data[sample(1:nrow(data),size=nrow(data),rep=T),]
  lm1=lm(formula(final.MrgETS),data=data1)
  
  #exporting results
  x1=summary(lm1)
  coeff=data.frame(Estimate=x1$coefficients[,1],Std..Error=x1$coefficients[,2],t.value=x1$coefficients[,3],var=rownames(x1$coefficients))
  r2=data.frame(r2m=x1$r.squared,r2c=NA,warnM=NA,sing=NA)
  if(m==1){
    r2_ok=r2
    coeff_ok=coeff
  }else{
    coeff_ok=rbind(coeff_ok,coeff)
    r2_ok=rbind(r2_ok,r2)
  }
}
write.table(r2_ok,"r2_MrgETS.csv",sep=";",dec=".",row=F)
  #output:
    #r2m= marginal R2 of every bootstrap model (each line is a different model)
    #r2c= conditional R2 of every bootstrap model
    #warnM= inform for warning during the model fit (such as convergence issue)
    #sing= output of the model singularity test (FALSE= no signularity; TRUE= singularity issue)
  
write.table(coeff_ok,"coeff_MrgETS.csv",sep=";",dec=".",row=F)
  #output is the result of model$coefficient for everey bootstrap model:
    #Estimate= coefficient value computed from every boostrap model
    #Std..Error= Standard error computed from every boostrap model
    #t.value= t-value computed from every boostrap model
    #var= names of the term


#################################
#### Bootstrap model analysis####
#################################
rm(list=ls())
gc(reset=T)

#repository/file location
rep_data="....../Script" #unzip Scrip.zip available at figshare
rep_out="......." #inform the location where all results have to be saved

setwd(rep_out)
load("analysis_Class.RData")
setwd("model_cl6000")

mod=c("oLTN","oLTS","oLMN","oLMS","oETN","oETS","MrgLTN","MrgLMN","MrgLMS","MrgETN")

#bootstrap test: H0= coefficient equal to 0; H1= coefficient higher than 0
test1=function(x,mu=0){
  return(1-(length(x[x>mu])/length(x)))
}

#bootstrap test: H0= coefficient equal to 0; H1= coefficient lower than 0
test2=function(x,mu=0){
  return(1-(length(x[x<mu])/length(x)))
}

a=1
s1=10
for(i in 1:length(mod)){
  print(a)
  r2=read.csv2(paste("r2_",mod[i],".csv",sep=""),sep=";",dec=".")
  r2$nB=1:nrow(r2)
  r2a=subset(r2,warnM=="boundary (singular) fit: see ?isSingular" | is.na(warnM)==T)
  print(paste(mod[i],': ',nrow(r2a),sep=""))
  if(nrow(r2a)>=5000){ #random selection of 5000 bootstrap model (criteria: model convergence)
    set.seed(s1)
    r2a=r2a[sample(1:nrow(r2a),size=5000,replace=F),]
  }
  if(nrow(r2a)<5000){
    print("sampling size is less than 5000 bootstraps")
  }else{
    resR2=data.frame(type="r2m",moy=mean(r2a$r2m),sd=sd(r2a$r2m),median=median(r2a$r2m),q025=quantile(r2a$r2m,probs=0.025),p975=quantile(r2a$r2m,probs=0.975))
    resR2=rbind(resR2,data.frame(type="r2c",moy=mean(r2a$r2c),sd=sd(r2a$r2c),median=median(r2a$r2c),q025=quantile(r2a$r2c,probs=0.025),p975=quantile(r2a$r2c,probs=0.975)))
    resR2$n_sing=nrow(subset(r2a,sing==T))
    resR2$model=mod[i]
    
    c1=read.csv2(paste("coeff_",mod[i],".csv",sep=""),sep=";",dec=".",h=T)
    c1$nB=rep(1:6000,each=nrow(c1)/6000)
    c1=merge(c1,data.frame(nB=r2a$nB),by.x="nB",by.y="nB")
    cX=c1
    c1$int=grepl("PositionTrailing edge:",as.character(c1$var))
    if(grepl("Mrg",mod[i])==T){
      if(nrow(subset(c1,int==T))>0){
        c1a=subset(c1,var=="(Intercept)") #intercept général
        c1b=subset(c1,var=="PositionTrailing edge") #intercept général
        rr3=subset(c1,var!="PositionTrailing edge" & var!="(Intercept)" & int==F)
        rr3$Class=gsub("Class","",as.character(rr3$var))
        rr3=merge(rr3,data.frame(nB=c1a$nB,Intercept=c1a$Estimate),by.x="nB",by.y="nB")
        rr3$val=rr3$Estimate+rr3$Intercept #computing class range shifts at the leading edge
        rr3=rr3[,-8]
        c1a$val=c1a$Estimate
        c1a$Class=NA
        rr3=rbind(c1a,rr3)
        rr3$pos="LE"
        rr3=rr3[order(rr3$nB),]
        rr3$Class[is.na(rr3$Class)==T]=as.character(class_ref$ref[class_ref==mod[i]][[1]])
        resLE=data.frame(Class=names(tapply(rr3$val,rr3$Class,mean)),moy=tapply(rr3$val,rr3$Class,mean),sd=tapply(rr3$val,rr3$Class,sd),median=tapply(rr3$val,rr3$Class,median),q025=tapply(rr3$val,rr3$Class,quantile,probs=0.025),
                         p975=tapply(rr3$val,rr3$Class,quantile,probs=0.975), pv.inf0=tapply(rr3$val,rr3$Class,test2,mu=0),pv.sup0=tapply(rr3$val,rr3$Class,test1,mu=0))
        resLE$pos="LE"
        
        rr4=subset(c1,int==T)
        rr4$Class=gsub("PositionTrailing edge:Class","",as.character(rr4$var))
        rr4$v2=paste(rr4$Class,"_",rr4$nB,sep="")
        rr3$v2=paste(rr3$Class,"_",rr3$nB,sep="")
        rr4=merge(rr4,data.frame(v2=rr3$v2,val=rr3$val),by.x="v2",by.y="v2",all.x=T)
        rr4$val[is.na(rr4$val)]=0
        rr4$val2=rr4$Estimate+rr4$val #computing class range shifts at the trailing edge
        rr4=merge(rr4,data.frame(nB=c1b$nB,val3=c1b$Estimate),by.x="nB",by.y="nB",all.x=T)
        rr4$val2=rr4$val2+rr4$val3
        rr4=rr4[,c(-2,-9,-11)]
        names(rr4)[ncol(rr4)]="val"
        rr4$pos="TE"
        rr3=rr3[,-1*ncol(rr3)]
        
        rr4a=subset(rr3,var=="(Intercept)")
        rr4a=merge(rr4a,data.frame(nB=c1b$nB,val3=c1b$Estimate),by.x="nB",by.y="nB",all.x=T)
        rr4a$val=rr4a$val+rr4a$val3
        rr4a=rr4a[,-1*ncol(rr4a)]
        rr4a$pos="TE"
        rr4=rbind(rr4,rr4a)
        rr4=rr4[order(rr4$nB),]
        rr3=rr3[order(rr3$nB),]
        rr4$Class[is.na(rr4$Class)==T]=as.character(class_ref$ref[class_ref==mod[i]][[1]])
        resTE=data.frame(Class=names(tapply(rr4$val,rr4$Class,mean)),moy=tapply(rr4$val,rr4$Class,mean),sd=tapply(rr4$val,rr4$Class,sd),median=tapply(rr4$val,rr4$Class,median),q025=tapply(rr4$val,rr4$Class,quantile,probs=0.025),
                         p975=tapply(rr4$val,rr4$Class,quantile,probs=0.975), pv.inf0=tapply(rr4$val,rr4$Class,test2,mu=0),pv.sup0=tapply(rr4$val,rr4$Class,test1,mu=0))
        resTE$pos="TE"
        rr3=rbind(rr3,rr4)
        rr3$model=mod[i]
        res=rbind(resTE,resLE)
        res$model=mod[i]
        cX$pos="LETE"
        
      }else{
        c1a=subset(c1,var=="(Intercept)") #intercept général
        c1b=subset(c1,var=="PositionTrailing edge") #intercept général
        rr3=subset(c1,var!="PositionTrailing edge" & var!="(Intercept)")
        rr3$Class=gsub("Class","",as.character(rr3$var))
        rr3=merge(rr3,data.frame(nB=c1a$nB,Intercept=c1a$Estimate),by.x="nB",by.y="nB")
        rr3$val=rr3$Estimate+rr3$Intercept #computing class range shifts at the leading edge
        rr3=rr3[,-8]
        c1a$val=c1a$Estimate
        c1a$Class=NA
        rr3=rbind(c1a,rr3)
        rr3$pos="LE"
        rr3=rr3[order(rr3$nB),]
        rr3$Class[is.na(rr3$Class)==T]=as.character(class_ref$ref[class_ref==mod[i]][[1]])
        resLE=data.frame(Class=names(tapply(rr3$val,rr3$Class,mean)),moy=tapply(rr3$val,rr3$Class,mean),sd=tapply(rr3$val,rr3$Class,sd),median=tapply(rr3$val,rr3$Class,median),q025=tapply(rr3$val,rr3$Class,quantile,probs=0.025),
                         p975=tapply(rr3$val,rr3$Class,quantile,probs=0.975), pv.inf0=tapply(rr3$val,rr3$Class,test2,mu=0),pv.sup0=tapply(rr3$val,rr3$Class,test1,mu=0))
        resLE$pos="LE"
        
        rr4=rr3
        rr4=merge(rr4,data.frame(nB=c1b$nB,TE=c1b$Estimate),by.x="nB",by.y="nB")
        rr4$val=rr4$val+rr4$TE #computing class range shifts at the trailing edge
        rr4$pos="TE"
        rr4=rr4[,-1*ncol(rr4)]
        rr4=rr4[order(rr4$nB),]
        resTE=data.frame(Class=names(tapply(rr4$val,rr4$Class,mean)),moy=tapply(rr4$val,rr4$Class,mean),sd=tapply(rr4$val,rr4$Class,sd),median=tapply(rr4$val,rr4$Class,median),q025=tapply(rr4$val,rr4$Class,quantile,probs=0.025),
                         p975=tapply(rr4$val,rr4$Class,quantile,probs=0.975), pv.inf0=tapply(rr4$val,rr4$Class,test2,mu=0),pv.sup0=tapply(rr4$val,rr4$Class,test1,mu=0))
        resTE$pos="TE"
        rr3=rbind(rr3,rr4)
        rr3$model=mod[i]
        res=rbind(resTE,resLE)
        res$model=mod[i]
        cX$pos="LETE"
      }
    }else{
      c1a=subset(c1,var=="(Intercept)")
      rr3=subset(c1,var!="(Intercept)")
      if(nrow(rr3)>0){
        rr3$Class=gsub("Class","",as.character(rr3$var))
        rr3=merge(rr3,data.frame(nB=c1a$nB,Intercept=c1a$Estimate),by.x="nB",by.y="nB")
        rr3$val=rr3$Estimate+rr3$Intercept
        rr3=rr3[,-8]
        c1a$val=c1a$Estimate
        c1a$Class=NA
        rr3=rbind(c1a,rr3)
      }
      else{
        rr3=c1a
        rr3$val=rr3$Estimate
        rr3$Class=NA
      }
      rr3$Class[is.na(rr3$Class)==T]=as.character(class_ref$ref[class_ref==mod[i]][[1]])
      rr3$pos="CE"
      rr3=rr3[order(rr3$nB),]
      res=data.frame(Class=names(tapply(rr3$val,rr3$Class,mean)),moy=tapply(rr3$val,rr3$Class,mean),sd=tapply(rr3$val,rr3$Class,sd),median=tapply(rr3$val,rr3$Class,median),q025=tapply(rr3$val,rr3$Class,quantile,probs=0.025),
                     p975=tapply(rr3$val,rr3$Class,quantile,probs=0.975), pv.inf0=tapply(rr3$val,rr3$Class,test2,mu=0),pv.sup0=tapply(rr3$val,rr3$Class,test1,mu=0))
      res$pos="CE"
      cX$pos="CE"
      rr3$model=mod[i]
      res$model=mod[i]
    }
    cX$model=mod[i]
    if(a==1){
      write.table(res,"datashift_class.csv",sep=";",dec=".",append=F,row.names=F,col.names=T)  #used to make Tables S4
      #range shift statistics for each taxonomic class
        #Class =  Taxonomic class
        #moy = mean range shift (in km/yr and m/yr along the latitudinal and elevational gradient) computed from the 5000 bootstrapped mode
        #sd = standard deviation of range shift computed from the 5000 bootstrapped model
        #med = median coefficient computed from the 5000 bootstrapped model
        #q025 = 2.5% quantile of the bootstrap distribution of range shifts
        #q975 = 97.5% quantile of the bootstrap distribution of range shifts
        #pv.sup0 = p-value of the bootstrap test testing the significance of the range shifts (H0= coefficient equal to 0; H1= coefficient greater than 0)
        #pv.inf0 = p-value of the bootstrap test testing the significance of the range shifts (H0= coefficient equal to 0; H1= coefficient lower than 0)
        #pos= range shift position (CE= centroid; LE=Leading edge; TE=Trailing edge)
        #model= model name abbreviation
      
      write.table(rr3,"bootshift_class.csv",sep=";",dec=".",append=F,row.names=F,col.names=T)  #used to draw Fig2
      #range shift estimation for each taxonomic class and each bootstrap model
        #nB= ID of the bootstrap model
        #val= range shift estimation
        #Class= Taxonomic class
        #pos= range shift position (CE= centroid; LE=Leading edge; TE=Trailing edge)
        #model= model name abbreviation
      
      
      write.table(resR2,"R2summary.csv",sep=";",dec=".",append=F,row.names=F,col.names=T) #used to make Table S3
      #R2 statistic computed from 5000 boostrap model:
        #type= R2 type (r2m= marginal R2; r2c= conditional R2)
        #moy= mean R2 value
        #sd= R2 standard deviation value
        #median= median R2 value
        #q025= 2.5% quantile of the bootstrap distribution
        #q975= 97.5% quantile of the bootstrap distribution
        #n_sing= number of singular model
        #model= model name abbreviation
      
      #write.table(cX,"bootcoeff.csv",sep=";",dec=".",append=F,row.names=F,col.names=T)
        
      
    }else{
      write.table(res,"datashift_class.csv",sep=";",dec=".",append=T,row.names=F,col.names=F) #used to make Tables S4
      write.table(rr3,"bootshift_class.csv",sep=";",dec=".",append=T,row.names=F,col.names=F) # used to draw Fig2
      write.table(resR2,"R2summary.csv",sep=";",dec=".",append=T,row.names=F,col.names=F) #used to make Table S3
     # write.table(cX,"bootcoeff.csv",sep=";",dec=".",append=T,row.names=F,col.names=F)
    }
    a=a+1
  }
}

#adding output of the simple linear model MrgETS
s1=10
mod="MrgETS"
for(i in 1:length(mod)){
  print(a)
  r2=read.csv2(paste("r2_",mod[i],".csv",sep=""),sep=";",dec=".")
  r2$nB=1:nrow(r2)
  r2a=subset(r2,warnM=="boundary (singular) fit: see ?isSingular" | is.na(warnM)==T)
  print(paste(mod[i],': ',nrow(r2a),sep=""))
  if(nrow(r2a)>=5000){ #random selection of 5000 bootstrap model (criteria: model convergence)
    set.seed(s1)
    r2a=r2a[sample(1:nrow(r2a),size=5000,replace=F),]
  }
  if(nrow(r2a)<5000){
    print("sampling size is less than 5000 bootstraps")
  }else{
    resR2=data.frame(type="r2m",moy=mean(r2a$r2m),sd=sd(r2a$r2m),median=median(r2a$r2m),q025=quantile(r2a$r2m,probs=0.025),p975=quantile(r2a$r2m,probs=0.975))
    resR2=rbind(resR2,data.frame(type="r2c",moy=NA,sd=NA,median=NA,q025=NA,p975=NA))
    resR2$n_sing=NA
    resR2$model=mod[i]
    
    c1=read.csv2(paste("coeff_",mod[i],".csv",sep=""),sep=";",dec=".",h=T)
    c1$nB=rep(1:6000,each=nrow(c1)/6000)
    c1=merge(c1,data.frame(nB=r2a$nB),by.x="nB",by.y="nB")
    cX=c1
    c1$int=grepl("PositionTrailing edge:",as.character(c1$var))
    if(grepl("Mrg",mod[i])==T){
      if(nrow(subset(c1,int==T))>0){
        c1a=subset(c1,var=="(Intercept)") #intercept général
        c1b=subset(c1,var=="PositionTrailing edge") #intercept général
        rr3=subset(c1,var!="PositionTrailing edge" & var!="(Intercept)" & int==F)
        rr3$Class=gsub("Class","",as.character(rr3$var))
        rr3=merge(rr3,data.frame(nB=c1a$nB,Intercept=c1a$Estimate),by.x="nB",by.y="nB")
        rr3$val=rr3$Estimate+rr3$Intercept #computing class range shifts at the leading edge
        rr3=rr3[,-8]
        c1a$val=c1a$Estimate
        c1a$Class=NA
        rr3=rbind(c1a,rr3)
        rr3$pos="LE"
        rr3=rr3[order(rr3$nB),]
        rr3$Class[is.na(rr3$Class)==T]=as.character(class_ref$ref[class_ref==mod[i]][[1]])
        resLE=data.frame(Class=names(tapply(rr3$val,rr3$Class,mean)),moy=tapply(rr3$val,rr3$Class,mean),sd=tapply(rr3$val,rr3$Class,sd),median=tapply(rr3$val,rr3$Class,median),q025=tapply(rr3$val,rr3$Class,quantile,probs=0.025),
                         p975=tapply(rr3$val,rr3$Class,quantile,probs=0.975), pv.inf0=tapply(rr3$val,rr3$Class,test2,mu=0),pv.sup0=tapply(rr3$val,rr3$Class,test1,mu=0))
        resLE$pos="LE"
        
        rr4=subset(c1,int==T)
        rr4$Class=gsub("PositionTrailing edge:Class","",as.character(rr4$var))
        rr4$v2=paste(rr4$Class,"_",rr4$nB,sep="")
        rr3$v2=paste(rr3$Class,"_",rr3$nB,sep="")
        rr4=merge(rr4,data.frame(v2=rr3$v2,val=rr3$val),by.x="v2",by.y="v2",all.x=T)
        rr4$val[is.na(rr4$val)]=0
        rr4$val2=rr4$Estimate+rr4$val #computing class range shifts at the trailing edge
        rr4=merge(rr4,data.frame(nB=c1b$nB,val3=c1b$Estimate),by.x="nB",by.y="nB",all.x=T)
        rr4$val2=rr4$val2+rr4$val3
        rr4=rr4[,c(-2,-9,-11)]
        names(rr4)[ncol(rr4)]="val"
        rr4$pos="TE"
        rr3=rr3[,-1*ncol(rr3)]
        
        rr4a=subset(rr3,var=="(Intercept)")
        rr4a=merge(rr4a,data.frame(nB=c1b$nB,val3=c1b$Estimate),by.x="nB",by.y="nB",all.x=T)
        rr4a$val=rr4a$val+rr4a$val3
        rr4a=rr4a[,-1*ncol(rr4a)]
        rr4a$pos="TE"
        rr4=rbind(rr4,rr4a)
        rr4=rr4[order(rr4$nB),]
        rr3=rr3[order(rr3$nB),]
        rr4$Class[is.na(rr4$Class)==T]=as.character(class_ref$ref[class_ref==mod[i]][[1]])
        resTE=data.frame(Class=names(tapply(rr4$val,rr4$Class,mean)),moy=tapply(rr4$val,rr4$Class,mean),sd=tapply(rr4$val,rr4$Class,sd),median=tapply(rr4$val,rr4$Class,median),q025=tapply(rr4$val,rr4$Class,quantile,probs=0.025),
                         p975=tapply(rr4$val,rr4$Class,quantile,probs=0.975), pv.inf0=tapply(rr4$val,rr4$Class,test2,mu=0),pv.sup0=tapply(rr4$val,rr4$Class,test1,mu=0))
        resTE$pos="TE"
        rr3=rbind(rr3,rr4)
        rr3$model=mod[i]
        res=rbind(resTE,resLE)
        res$model=mod[i]
        
        cX$pos="LETE"
        cX$model=mod[i]
        write.table(res,"datashift_class.csv",sep=";",dec=".",append=T,row.names=F,col.names=F) #used to make Tables S4
        write.table(rr3,"bootshift_class.csv",sep=";",dec=".",append=T,row.names=F,col.names=F) # used to draw Fig2
        write.table(resR2,"R2summary.csv",sep=";",dec=".",append=T,row.names=F,col.names=F) #used to make Table S3
        #write.table(cX,"bootcoeff.csv",sep=";",dec=".",append=T,row.names=F,col.names=F)
        
      }
    }
  }
}


