#work conducted on R version 3.6.1 (2019-07-05); Platform: x86_64-pc-linux-gnu (64-bit); Running under: Ubuntu 16.04.5 LTS
#authors: Bertrand Romain, Lise Comte, Bourgaud Luana and Jonathan Lenoir
#AIM= Exploring the relationship between species range shifts and velocity of climate warming, baseline temperature and human pressure.

#required packages
#library(nlme) #version 3.1.140
library(MuMIn) #version 1.43.6
library(lme4) #version 1.1.21
library(optimx) #version 2018-7.10
library(jtools) #version 2.0.1

#repository/file location
rep_data="....../Script" #unzip Scrip.zip available at figshare
rep_out="......." #inform the location where all results have to be saved

#############################
#### DATASET preparation ####
#############################
setwd(rep_data)
rSdata = read.table("Table_S1.csv",sep=";",h=T,dec=".",stringsAsFactors = FALSE) #exported from Table_S1.xlsx with separator=";" and decimal="."
rSdata = rSdata[order(rSdata$n_sp),]

#transforming continuous method variables in qualitative variables
rSdata$NtaxaF = as.numeric(cut(rSdata$Ntaxa,quantile(rSdata$Ntaxa,c(0,0.25,0.5,0.75,1)),include.lowest = T))  #75th quantile is also the max
rSdata$StartF = as.numeric(cut(rSdata$Start,quantile(rSdata$Start,c(0,0.25,0.5,0.75,1)),include.lowest = T))  #75th quantile is also the max
rSdata$AreaF = as.numeric(cut(rSdata$Area,quantile(rSdata$Area,c(0,0.25,0.5,0.75,1)),include.lowest = T))  #75th quantile is also the max

#relevel variable to limit singularity issue
rSdata$Sampling = ifelse(rSdata$Sampling == "TWO","TWO","MULT")
rSdata$Quality = ifelse(rSdata$Quality == "BALANCED","RESURVEYED",rSdata$Quality)

# variable pre-selection
chosen_varlatT=c("ShiftR", "LatVeloT", "BaselineT", "HFI", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif", "Source", "Class", "Family", "Genus", "Species", "Gradient", "Ecosystem", "LifeForm")
chosen_varlatM=c("ShiftR", "LatVeloT", "BaselineT", "HFI", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif", "Source", "Class", "Family", "Genus", "Species", "Gradient", "Ecosystem", "LifeForm")
chosen_varele=c("ShiftR", "EleVeloT", "BaselineT", "HFI", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif", "Source", "Class", "Family", "Genus", "Species", "Gradient", "Ecosystem", "LifeForm")

#####################################################
#### MODELING: random effect structure selection ####
#####################################################
#Elevational range shift
Class <- as.data.frame(table(rSdata$Class[rSdata$Gradient == "Elevation"]))
names(Class) <- c("Class", "Freq")
Class=Class[which(Class$Freq>30), ]       # Criteria: Class > 30 obs
Class
#           Class Freq
#  Actinopterygii   96
#        Amphibia  197
#            Aves 2227
#       Bryopsida   86
#         Insecta 1418
# Lecanoromycetes   65
#      Liliopsida 1574
#  Lycopodiopsida   41
#   Magnoliopsida 7065
#        Mammalia  355
#       Pinopsida  181
#  Polypodiopsida  182
#        Reptilia   33

# Selecting observation and variables to analyse range shifts
data <- rSdata[rSdata$Gradient == "Elevation", chosen_varele]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data$LifeForm <- as.factor(as.character(data$LifeForm))
data$StartF=as.factor(data$StartF)
data$NtaxaF=as.factor(data$NtaxaF)
data$AreaF=as.factor(data$AreaF)
data <- na.omit(data)
data=droplevels(data)
dim(data) # 13459 obs

#scaling explaining continuous variables
data$BaselineT2=gscale(data$BaselineT)^2
data$g_velT=gscale(data$EleVeloT)
data$g_basT=gscale(data$BaselineT)
data$g_HFI=gscale(data$HFI)

t1=data.frame(table(data$LifeForm))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$LifeForm=relevel(data$LifeForm,ref=ref1) #Intercept value of the model will be the range shift estimation of the life form having the greatest number of observation

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
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)
table(data$PrAb)
table(data$Sampling)  
table(data$Quality)
table(data$Signif)
table(data$Grain)

#selection of the best random effect structure
x=c("(1|AreaF)","(1|StartF)","(1|NtaxaF)","(1|PrAb)","(1|Sampling)","(1|Grain)","(1|Quality)","(1|Signif)")
f="ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x)) #take some time
setwd(rep_out)
write.table(tx,"sing_ELE.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x)) #take some time
  setwd(rep_out)
  write.table(tx,"sing_ELE.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm"
  tx=testSingAl(x,f,data,n=1:length(x)) #take some time
  setwd(rep_out)
  write.table(tx,"sing_ELE.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm, data)
        }
      }
    }
  }
}

if(isSingular(mod1)==F){
  mod2=lmer(formula(mod1),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  ms1=dredge(mod2) #take some time: ms1 is the table reporting results of the model selection through AICc
  write.table(ms1,"elevation_selAIC.csv",sep=";",dec=".",row=F)
  ms1E=ms1
  
  #best model selection with no singularity
  a=1
  iS=T
  while(iS==T){
    print(a)
    m1=get.models(ms1E,subset=a)[[1]]
    mbe=lmer(formula(m1), data=data,REML=T,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mbe)[[1]]
    a=a+1
  }
  summary(mbe)
  r.squaredGLMM(mbe)
}

#Latitudinal and marine range shifts
Class <- as.data.frame(table(rSdata$Class[rSdata$Gradient == "Latitudinal" & rSdata$Ecosystem=="Marine"]))
names(Class) <- c("Class", "Freq")
Class=Class[which(Class$Freq>30), ]       # Criteria: Class > 30 obs
Class
#           Class Freq
#  Actinopterygii  720
#        Bivalvia   71
#  Chondrichthyes   83
# Florideophyceae   70
#      Gastropoda  134
#    Malacostraca  123
#     Maxillopoda   55
#    Phaeophyceae   59
#      Polychaeta  143

# Selecting observation and variables to analyse range shifts
data <- rSdata[rSdata$Gradient == "Latitudinal" & rSdata$Ecosystem=="Marine", chosen_varlatM]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data$LifeForm <- as.factor(as.character(data$LifeForm))
data$StartF=as.factor(data$StartF)
data$NtaxaF=as.factor(data$NtaxaF)
data$AreaF=as.factor(data$AreaF)
data <- na.omit(data)
data=droplevels(data)
dim(data) # 1403 obs

#scaling explaining continuous variables
data$BaselineT2=gscale(data$BaselineT)^2
data$g_velT=gscale(data$LatVeloT)
data$g_basT=gscale(data$BaselineT)
data$g_HFI=gscale(data$HFI)

t1=data.frame(table(data$LifeForm))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$LifeForm=relevel(data$LifeForm,ref=ref1) #Intercept value of the model will be the range shift estimation of the life form having the greatest number of observation

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)
table(data$PrAb)
table(data$Sampling)  
table(data$Quality)
table(data$Signif)
table(data$Grain)

#selection of the best random effect structure
x=c("(1|AreaF)","(1|StartF)","(1|NtaxaF)","(1|PrAb)","(1|Sampling)","(1|Grain)","(1|Quality)","(1|Signif)") #set of random effects
f="ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x)) #take some time
setwd(rep_out)
write.table(tx,"sing_LATM.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x)) #take some time
  setwd(rep_out)
  write.table(tx,"sing_LATM.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm"
  tx=testSingAl(x,f,data,n=1:length(x)) #take some time
  setwd(rep_out)
  write.table(tx,"sing_LATM.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm, data)
        }
      }
    }
  }
}

if(isSingular(mod1)==F){
  mod2=lmer(formula(mod1),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  ms1=dredge(mod2) #take some time: ms1 is the table reporting results of the model selection through AICc
  setwd(rep_out)
  write.table(ms1,"latM_selAIC.csv",sep=";",dec=".",row=F)
  ms1Lm=ms1
  
  #best model selection with no singularity
  a=1
  iS=T
  while(iS==T){
    print(a)
    m1=get.models(ms1Lm,subset=a)[[1]]
    mblm=lmer(formula(m1), data=data,REML=T,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mblm)[[1]]
    a=a+1
  }
  summary(mblm)
  r.squaredGLMM(mblm)
}


#Latitudinal and terrestrial range shifts
Class <- as.data.frame(table(rSdata$Class[rSdata$Gradient == "Latitudinal" & rSdata$Ecosystem=="Terrestrial"]))
names(Class) <- c("Class", "Freq")
Class=Class[which(Class$Freq>30), ]       # Criteria: Class > 30 obs
Class
#          Class Freq
#       Amphibia  211
#      Arachnida  444
#           Aves 4882
#      Bryopsida  288
#  Equisetopsida   34
#        Insecta 4237
#     Liliopsida 1109
#  Magnoliopsida 3617
#   Malacostraca   37
#      Pinopsida   98
# Polypodiopsida  105
#       Reptilia   90

# Selecting observation and variables to analyse range shifts
data <- rSdata[rSdata$Gradient == "Latitudinal" & rSdata$Ecosystem=="Terrestrial", chosen_varlatT]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data$LifeForm <- as.factor(as.character(data$LifeForm))
data$StartF=as.factor(data$StartF)
data$NtaxaF=as.factor(data$NtaxaF)
data$AreaF=as.factor(data$AreaF)
data <- na.omit(data)
data=droplevels(data)
dim(data) # 15118 obs

#scaling explaining continuous variables
data$BaselineT2=gscale(data$BaselineT)^2
data$g_velT=gscale(data$LatVeloT)
data$g_basT=gscale(data$BaselineT)
data$g_HFI=gscale(data$HFI)

t1=data.frame(table(data$LifeForm))
t1=t1[order(t1$Freq,decreasing=T),]
ref1=as.character(t1[1,1])
data$LifeForm=relevel(data$LifeForm,ref=ref1) #Intercept value of the model will be the range shift estimation of the life form having the greatest number of observation

#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)
table(data$PrAb)
table(data$Sampling)  
table(data$Quality)
table(data$Signif)
table(data$Grain)

#selection of the best random effect structure
x=c("(1|AreaF)","(1|StartF)","(1|NtaxaF)","(1|PrAb)","(1|Sampling)","(1|Grain)","(1|Quality)","(1|Signif)")
f="ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Family/Genus)"
tx=testSingAl(x,f,data,n=1:length(x)) #take some times
setwd(rep_out)
write.table(tx,"sing_LATT.csv",sep=";",dec=".",row=F)
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
  f="ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Genus)" #
  tx=testSingAl(x,f,data,n=1:length(x)) #take some times
  setwd(rep_out)
  write.table(tx,"sing_LATT.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
  f="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm"
  tx=testSingAl(x,f,data,n=1:length(x)) #take some times
  setwd(rep_out)
  write.table(tx,"sing_LATT.csv",sep=";",dec=".",row=F,col.names=F,append=T)
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
    f1="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Family/Genus)"
    mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
      f1="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Genus)"
      mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Source)
        f1="ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm + (1|Source)"
        mod1=lmer(as.formula(f1),data,REML=TRUE,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
        if(isSingular(mod1)==T){ #in case of singularity, we finally fit a simple linear model
          print("Impossible to fit the linear mixed-effect model due to singularity issue")
          mod1=lm(ShiftR ~ ShiftR ~ BaselineT2 + g_velT * g_basT + g_velT * g_HFI + g_velT * LifeForm, data)
        }
      }
    }
  }
}

if(isSingular(mod1)==F){
  mod2=lmer(formula(mod1),data,REML=F,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  ms1=dredge(mod2) #take some time: ms1 is the table reporting results of the model selection through AICc
  setwd(rep_out)
  write.table(ms1,"latT_selAIC.csv",sep=";",dec=".",row=F)
  ms1Lt=ms1
  
  #best model selection with no singularity
  a=1
  iS=T
  while(iS==T){
    print(a)
    m1=get.models(ms1Lt,subset=a)[[1]]
    mblt=lmer(formula(m1), data=data,REML=T,na.action="na.fail",control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    iS=isSingular(mblt)[[1]]
    a=a+1
  }
  summary(mblt)
  r.squaredGLMM(mblt)
}

setwd(rep_out)
save.image("analysis_SpRangeShift.RData")

#############################
#### MODELING: bootstrap ####
#############################
rm(list=ls())
gc(reset=T)

#required packages
library(MuMIn) #version 1.43.6
library(lme4) #version 1.1.21
library(optimx) #version 2018-7.10
library(data.table) #version 1.12.2
library(jtools) #version 2.0.1

#repository/file location
rep_data="....../Script" #unzip Scrip.zip available at figshare
rep_out="......." #inform the location where all results have to be saved

setwd(rep_out)
load("analysis_SpRangeShift.RData")

dir.create("model_sp6000")
setwd("model_sp6000")

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
names=c("elev","latM","latT")
models=list(mbe,mblm,mblt)

nB=6000 #set the bootstrap number
s1=10 #the seed
ncpu=10 #the number of cpu used for parallel computation

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

#################################
#### Bootstrap model analysis####
#################################
rm(list=ls())
gc(reset=T)

#required packages
library(MuMIn) #version 1.43.6
library(lme4) #version 1.1.21
library(optimx) #version 2018-7.10
library(data.table) #version 1.12.2
library(jtools) #version 2.0.1

#repository/file location
rep_data="....../Script" #unzip Scrip.zip available at figshare
rep_out="......." #inform the location where all results have to be saved

setwd(rep_out)
load("analysis_SpRangeShift.RData")
setwd("model_sp6000")

mod=c("elev","latM","latT")

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
    resR2=data.frame(type="r2m",moy=mean(r2a$r2m),sd=sd(r2a$r2m),median=median(r2a$r2m),q025=quantile(r2a$r2m,probs=0.025),q975=quantile(r2a$r2m,probs=0.975))
    resR2=rbind(resR2,data.frame(type="r2c",moy=mean(r2a$r2c),sd=sd(r2a$r2c),median=median(r2a$r2c),q025=quantile(r2a$r2c,probs=0.025),p975=quantile(r2a$r2c,probs=0.975)))
    resR2$n_sing=nrow(subset(r2a,sing==T))
    resR2$model=mod[i]
    
    c1=read.csv2(paste("coeff_",mod[i],".csv",sep=""),sep=";",dec=".",h=T)
    c1$nB=rep(1:6000,each=nrow(c1)/6000)
    c1=merge(c1,data.frame(nB=r2a$nB),by.x="nB",by.y="nB")
    rr3=c1
    
    res=data.frame(var=names(tapply(rr3$Estimate,rr3$var,mean)),moy=tapply(rr3$Estimate,rr3$var,mean),sd=tapply(rr3$Estimate,rr3$var,sd),median=tapply(rr3$Estimate,rr3$var,median),q025=tapply(rr3$Estimate,rr3$var,quantile,probs=0.025),
                   p975=tapply(rr3$Estimate,rr3$var,quantile,probs=0.975), pv.inf0=tapply(rr3$Estimate,rr3$var,test2,mu=0),pv.sup0=tapply(rr3$Estimate,rr3$var,test1,mu=0))
    res$model=mod[i]
    if(a==1){
      resR2_ok=resR2
      res_ok=res
    }
    else{
      resR2_ok=rbind(resR2_ok,resR2)
      res_ok=rbind(res_ok,res)
    }
    write.table(rr3,paste("coeff_",mod[i],"_5000.csv",sep=""),sep=";",dec=".",row=F)
    #Coefficient value for each of the 5000 selected bootstrap models
      #nB= ID of the bootstrap model
      #Estimate= coefficient value computed from every boostrap model
      #Std..Error= Standard error computed from every boostrap model
      #t.value= t-value computed from every boostrap model
      #var= names of the term
    
    write.table(r2a,paste("r2_",mod[i],"_5000.csv",sep=""),sep=";",dec=".",row=F)
    #R2 values for each of the 5000 selected bootstrap models :
      #r2m= marginal R2 of every bootstrap model (each line is a different model)
      #r2c= conditional R2 of every bootstrap model
      #warnM= inform for warning during the model fit (such as convergence issue)
      #sing= output of the model singularity test (FALSE= no signularity; TRUE= singularity issue) 
      #nB= ID of the bootstrap model
    
    a=a+1
  }
}

write.table(resR2_ok,"summary_R2.csv",sep=";",dec=".",row=F) #used to make table S5
#R2 statistic computed from 5000 boostrap model:
  #type= R2 type (r2m= marginal R2; r2c= conditional R2)
  #moy= mean R2 value
  #sd= R2 standard deviation value
  #median= median R2 value
  #q025= 2.5% quantile of the bootstrap distribution
  #q975= 97.5% quantile of the bootstrap distribution
  #n_sing= number of singular model
  #model= model name abbreviation

write.table(res_ok,"summary_coeff.csv",sep=";",dec=".",row=F) #used to make table S6
#coefficient statistics
  #Var =  Term of the model
  #moy = mean effect computed from the 5000 bootstrap model
  #sd = standard deviation computed from the 5000 bootstrap model
  #med = median effect computed from the 5000 bootstrap model
  #q025 = 2.5% quantile of the bootstrap coefficient distribution
  #q975 = 97.5% quantile of the bootstrap coefficient distribution
  #pv.sup0 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 0; H1= coefficient greater than 0)
  #pv.inf0 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 0; H1= coefficient lower than 0)
  #model= model name abbreviation

####################################################################################################################################
#### Test of the effect of climate warming velocity, baseline temperature conditions and human pressure on species range shifts ####
####################################################################################################################################
rm(list=ls())
gc(reset=T)

library(data.table) #version 1.12

#repository/file location
rep_data="....../Script" #unzip Scrip.zip available at figshare
rep_out="......." #inform the location where all results have to be saved

setwd(rep_data)
rSdata = read.table("Table_S1.csv",sep=";",h=T,dec=".",stringsAsFactors = FALSE) #exported from Table_S1.xlsx with separator=";" and decimal="."
rSdata = rSdata[order(rSdata$n_sp),]

#transforming continuous method variables in qualitative variables
rSdata$NtaxaF = as.numeric(cut(rSdata$Ntaxa,quantile(rSdata$Ntaxa,c(0,0.25,0.5,0.75,1)),include.lowest = T))  #75th quantile is also the max
rSdata$StartF = as.numeric(cut(rSdata$Start,quantile(rSdata$Start,c(0,0.25,0.5,0.75,1)),include.lowest = T))  #75th quantile is also the max
rSdata$AreaF = as.numeric(cut(rSdata$Area,quantile(rSdata$Area,c(0,0.25,0.5,0.75,1)),include.lowest = T))  #75th quantile is also the max

#relevel variable to limit singularity issue
rSdata$Sampling = ifelse(rSdata$Sampling == "TWO","TWO","MULT")
rSdata$Quality = ifelse(rSdata$Quality == "BALANCED","RESURVEYED",rSdata$Quality)

# variable pre-selection
chosen_varlatT=c("ShiftR", "LatVeloT", "BaselineT", "HFI", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif", "Source", "Class", "Family", "Genus", "Species", "Gradient", "Ecosystem", "LifeForm")
chosen_varlatM=c("ShiftR", "LatVeloT", "BaselineT", "HFI", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif", "Source", "Class", "Family", "Genus", "Species", "Gradient", "Ecosystem", "LifeForm")
chosen_varele=c("ShiftR", "EleVeloT", "BaselineT", "HFI", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif", "Source", "Class", "Family", "Genus", "Species", "Gradient", "Ecosystem", "LifeForm")

#Elevation range shifts
Class <- as.data.frame(table(rSdata$Class[rSdata$Gradient == "Elevation"]))
names(Class) <- c("Class", "Freq")
Class=Class[which(Class$Freq>30), ]       # Criteria: Class > 30 obs
Class
#           Class Freq
#  Actinopterygii   96
#        Amphibia  197
#            Aves 2227
#       Bryopsida   86
#         Insecta 1418
# Lecanoromycetes   65
#      Liliopsida 1574
#  Lycopodiopsida   41
#   Magnoliopsida 7065
#        Mammalia  355
#       Pinopsida  181
#  Polypodiopsida  182
#        Reptilia   33

# Selecting observation and variables to analyse range shifts
data <- rSdata[rSdata$Gradient == "Elevation", chosen_varele]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data <- na.omit(data)
dim(data) # 13459 obs

setwd(rep_out)
setwd("model_sp6000")
t1=read.csv2("summary_coeff.csv",sep=";",dec=".",h=T)
t1a=subset(t1,model=="elev")
t2=read.csv2("coeff_elev_5000.csv",sep=";",dec=".",h=T)

v2=seq(-20,50,by=0.1) #range of baseline temperature tested
v1=seq(-10,20,by=0.1) #range of climate warming velocity tested

#Estimation of the effect of climate warming velocity (unscaled)
#Climate warming velocity parameters
#The effect of climate warming velocity depends to baseline temperatures
a2=t1a$moy[t1a$var=="g_velT"] #direct effect of climate warming velocity
a5=t1a$moy[t1a$var=="g_basT:g_velT"] #interaction between climate warming velocity and baseline temperatures

res=matrix(NA,ncol=4,nrow=length(v2))
#ectotherms
i0=t1a$moy[t1a$var=="(Intercept)"]+t1a$moy[t1a$var=="LifeFormecto"]
a3=t1a$moy[t1a$var=="g_velT:LifeFormecto"]
for(z in 1:length(v2)){
  res[z,1]=v2[z]
  res[z,2]=a3+a2+(a5*(v2[z]-mean(data$baselineT))/(2*sd(data$baselineT)))
  p1=res[z,2]*(v1-mean(data$EleVeloT))/(2*sd(data$EleVeloT)) + i0
  lm1=lm(p1~v1)
  res[z,3]=lm1$coefficients[[2]]
  res[z,4]=lm1$coefficients[[1]]
}
resM=data.frame(res)
resM=resM[,-2]
names(resM)=c("basT","slop_velT","int")
resM$LifeForm="ecto"
write.table(resM,"slopVelT_ELEV.csv",sep=";",dec=".",row=F)

#endotherms
i0=t1a$moy[t1a$var=="(Intercept)"]+t1a$moy[t1a$var=="LifeFormendo"]
a3=t1a$moy[t1a$var=="g_velT:LifeFormendo"]
for(z in 1:length(v2)){
  res[z,1]=v2[z]
  res[z,2]=a3+a2+(a5*(v2[z]-mean(data$baselineT))/(2*sd(data$baselineT)))
  p1=res[z,2]*(v1-mean(data$EleVeloT))/(2*sd(data$EleVeloT)) + i0
  lm1=lm(p1~v1)
  res[z,3]=lm1$coefficients[[2]]
  res[z,4]=lm1$coefficients[[1]]
}
resM=data.frame(res)
resM=resM[,-2]
names(resM)=c("basT","slop_velT","int")
resM$LifeForm="endo"
write.table(resM,"slopVelT_ELEV.csv",sep=";",dec=".",row=F,col.names=F,append=T)

#phanerogams
i0=t1a$moy[t1a$var=="(Intercept)"]
for(z in 1:length(v2)){
  res[z,1]=v2[z]
  res[z,2]=a2+(a5*(v2[z]-mean(data$baselineT))/(2*sd(data$baselineT)))
  p1=res[z,2]*(v1-mean(data$EleVeloT))/(2*sd(data$EleVeloT)) + i0
  lm1=lm(p1~v1)
  res[z,3]=lm1$coefficients[[2]]
  res[z,4]=lm1$coefficients[[1]]
}
resM=data.frame(res)
resM=resM[,-2]
names(resM)=c("basT","slop_velT","int")
resM$LifeForm="phanero"
write.table(resM,"slopVelT_ELEV.csv",sep=";",dec=".",row=F,col.names=F,append=T)

#cryptogams
i0=t1a$moy[t1a$var=="(Intercept)"]+t1a$moy[t1a$var=="LifeFormcrypto"]
a3=t1a$moy[t1a$var=="g_velT:LifeFormcrypto"]
for(z in 1:length(v2)){
  res[z,1]=v2[z]
  res[z,2]=a3+a2+(a5*(v2[z]-mean(data$baselineT))/(2*sd(data$baselineT)))
  p1=res[z,2]*(v1-mean(data$EleVeloT))/(2*sd(data$EleVeloT)) + i0
  lm1=lm(p1~v1)
  res[z,3]=lm1$coefficients[[2]]
  res[z,4]=lm1$coefficients[[1]]
}
resM=data.frame(res)
resM=resM[,-2]
names(resM)=c("basT","slop_velT","int")
resM$LifeForm="crypto"
write.table(resM,"slopVelT_ELEV.csv",sep=";",dec=".",row=F,col.names=F,append=T)

#Test of the effect of climate warming velocity in each boostrap model
nB=unique(t2$nB)
a=1
for(j in nB){ #take some time
  print(j)
  t2a=subset(t2,nB==j)
  t2a$velT=grepl("g_velT",as.character(t2a$var))
  t2a=subset(t2a,velT==T)
  #Climate warming velocity parameters computed from the focal model
    #The effect of climate warming velocity depends to the life forms and baseline temperature conditions
  a2=t2a$Estimate[t2a$var=="g_velT"] #direct effect of climate warming velocity for phanerogams
  ecto=t2a$Estimate[t2a$var=="g_velT:LifeFormecto"] #effect of climate warming velocity for ectotherms
  endo=t2a$Estimate[t2a$var=="g_velT:LifeFormendo"] #effect of climate warming velocity for endotherms
  crypto=t2a$Estimate[t2a$var=="g_velT:LifeFormcrypto"] #effect of climate warming velocity for cryptogams
  a4=t2a$Estimate[t2a$var=="g_basT:g_velT"] #interaction between climate warming velocity and baseline temperatures
  
  #Intercept
  t2a=subset(t2,nB==j)
  i_phanero=t2a$Estimate[t2a$var=="(Intercept)"] #for phanerogams
  i_crypto=t2a$Estimate[t2a$var=="(Intercept)"] + t2a$Estimate[t2a$var=="LifeFormcrypto"] #for cryptogams
  i_ecto=t2a$Estimate[t2a$var=="(Intercept)"] + t2a$Estimate[t2a$var=="LifeFormecto"] #for ectotherms
  i_endo=t2a$Estimate[t2a$var=="(Intercept)"] + t2a$Estimate[t2a$var=="LifeFormendo"] #for endotherms
  
  res=matrix(NA,ncol=4,nrow=length(v2))
  
  #phanerogams
  for(i in 1:length(v2)){
    res[i,1]=v2[i]
    res[i,2]=a2+(a4*(v2[i]-mean(data$BaselineT))/(2*sd(data$BaselineT)))
    p1=res[i,2]*(v1-mean(data$EleVeloT))/(2*sd(data$EleVeloT))+i_phanero
    lm1=lm(p1~v1)
    res[i,3]=lm1$coefficients[[2]]
    res[i,4]=lm1$coefficients[[1]]
  }
  resP=data.frame(res)
  resP=resP[,-2]
  names(resP)=c("basT","slop_velT","int")
  
  #ectotherms
  res=matrix(NA,ncol=4,nrow=length(v2))
  for(i in 1:length(v2)){
    res[i,1]=v2[i]
    res[i,2]=ecto+a2+(a4*(v2[i]-mean(data$BaselineT))/(2*sd(data$BaselineT)))
    p1=res[i,2]*(v1-mean(data$EleVeloT))/(2*sd(data$EleVeloT))+i_phanero+i_ecto
    lm1=lm(p1~v1)
    res[i,3]=lm1$coefficients[2]
    res[i,4]=lm1$coefficients[1]
  }
  resEc=data.frame(res)
  resEc=resEc[,-2]
  names(resEc)=c("basT","slop_velT","int")
  
  #endotherms
  res=matrix(NA,ncol=4,nrow=length(v2))
  for(i in 1:length(v2)){
    res[i,1]=v2[i]
    res[i,2]=endo+a2+(a4*(v2[i]-mean(data$BaselineT))/(2*sd(data$BaselineT)))
    p1=res[i,2]*(v1-mean(data$EleVeloT))/(2*sd(data$EleVeloT))+i_phanero+i_endo
    lm1=lm(p1~v1)
    res[i,3]=lm1$coefficients[2]
    res[i,4]=lm1$coefficients[1]
  }
  resEn=data.frame(res)
  resEn=resEn[,-2]
  names(resEn)=c("basT","slop_velT","int")
  
  #cryptogams
  res=matrix(NA,ncol=4,nrow=length(v2))
  for(i in 1:length(v2)){
    res[i,1]=v2[i]
    res[i,2]=crypto+a2+(a4*(v2[i]-mean(data$BaselineT))/(2*sd(data$BaselineT)))
    p1=res[i,2]*(v1-mean(data$EleVeloT))/(2*sd(data$EleVeloT))+i_phanero+i_crypto
    lm1=lm(p1~v1)
    res[i,3]=lm1$coefficients[2]
    res[i,4]=lm1$coefficients[1]
  }
  resC=data.frame(res)
  resC=resC[,-2]
  names(resC)=c("basT","slop_velT","int")
  
  resC$nB=j 
  resP$nB=j
  resEn$nB=j
  resEc$nB=j
  
  if(a==1){
    write.table(resC,"boostrap_slopTvel_crypto_ELEV.csv",sep=";",dec=".",row=F)
    write.table(resEc,"boostrap_slopTvel_ecto_ELEV.csv",sep=";",dec=".",row=F)
    write.table(resEn,"boostrap_slopTvel_endo_ELEV.csv",sep=";",dec=".",row=F)
    write.table(resP,"boostrap_slopTvel_phanero_ELEV.csv",sep=";",dec=".",row=F)
  }else{
    write.table(resC,"boostrap_slopTvel_crypto_ELEV.csv",sep=";",dec=".",row=F,col.names=F,append=T)
    write.table(resEc,"boostrap_slopTvel_ecto_ELEV.csv",sep=";",dec=".",row=F,col.names=F,append=T)
    write.table(resEn,"boostrap_slopTvel_endo_ELEV.csv",sep=";",dec=".",row=F,col.names=F,append=T)
    write.table(resP,"boostrap_slopTvel_phanero_ELEV.csv",sep=";",dec=".",row=F,col.names=F,append=T)
  }
  a=a+1
}
rm(resC,resEc,resEn,resP)
gc(reset=T)

#Statistics of the relationships for a range of combinations of baseline temperature conditions and climate warming velocity for each species life form
#bootstrap test: H0= coefficient equal to 0; H1= coefficient greater than 0
test1=function(x,mu=0){
  return(1-(length(x[x>mu])/length(x)))
}

#bootstrap test: H0= coefficient equal to 0; H1= coefficient lower than 0
test2=function(x,mu=0){
  return(1-(length(x[x<mu])/length(x)))
}

  #Cryptogams
t1=fread("boostrap_slopTvel_crypto_ELEV.csv",sep=";",dec=".",h=T)
gc(reset=T)
tx=data.frame(Tbas_moy=as.numeric(names(tapply(t1$slop_velT,t1$basT,mean))),moy=tapply(t1$slop_velT,t1$basT,mean),sd=tapply(t1$slop_velT,t1$basT,sd),q025=tapply(t1$slop_velT,t1$basT,quantile,probs=0.025),
              q975=tapply(t1$slop_velT,t1$basT,quantile,probs=0.975), pv.inf0=tapply(t1$slop_velT,t1$basT,test2,mu=0),pv.sup0=tapply(t1$slop_velT,t1$basT,test1,mu=0),pv.inf1=tapply(t1$slop_velT,t1$basT,test2,mu=1),pv.sup1=tapply(t1$slop_velT,t1$basT,test1,mu=1))
gc(reset=T)
tx$LifeForm="crypto"
t1a=subset(tx,pv.inf1>0.05 & pv.sup1>0.05)
range(t1a$Tbas_moy) #velocity slope never equal to 1
resOK=tx
rm(t1)
gc(reset=T)

  #Phanerogams
t1=fread("boostrap_slopTvel_phanero_ELEV.csv",sep=";",dec=".",h=T)
gc(reset=T)
tx=data.frame(Tbas_moy=as.numeric(names(tapply(t1$slop_velT,t1$basT,mean))),moy=tapply(t1$slop_velT,t1$basT,mean),sd=tapply(t1$slop_velT,t1$basT,sd),q025=tapply(t1$slop_velT,t1$basT,quantile,probs=0.025),
              q975=tapply(t1$slop_velT,t1$basT,quantile,probs=0.975), pv.inf0=tapply(t1$slop_velT,t1$basT,test2,mu=0),pv.sup0=tapply(t1$slop_velT,t1$basT,test1,mu=0),pv.inf1=tapply(t1$slop_velT,t1$basT,test2,mu=1),pv.sup1=tapply(t1$slop_velT,t1$basT,test1,mu=1))
gc(reset=T)
tx$LifeForm="phanero"
t1a=subset(tx,pv.inf1>0.05 & pv.sup1>0.05)
range(t1a$Tbas_moy) #velocity slope never equal to 1
resOK=rbind(resOK,tx)
rm(t1)
gc(reset=T)

  #Ectotherms
t1=fread("boostrap_slopTvel_ecto_ELEV.csv",sep=";",dec=".",h=T)
gc(reset=T)
tx=data.frame(Tbas_moy=as.numeric(names(tapply(t1$slop_velT,t1$basT,mean))),moy=tapply(t1$slop_velT,t1$basT,mean),sd=tapply(t1$slop_velT,t1$basT,sd),q025=tapply(t1$slop_velT,t1$basT,quantile,probs=0.025),
              q975=tapply(t1$slop_velT,t1$basT,quantile,probs=0.975), pv.inf0=tapply(t1$slop_velT,t1$basT,test2,mu=0),pv.sup0=tapply(t1$slop_velT,t1$basT,test1,mu=0),pv.inf1=tapply(t1$slop_velT,t1$basT,test2,mu=1),pv.sup1=tapply(t1$slop_velT,t1$basT,test1,mu=1))
gc(reset=T)
tx$LifeForm="ecto"
t1a=subset(tx,pv.inf1>0.05 & pv.sup1>0.05)
range(t1a$Tbas_moy) #velocity slope never equal to 1
resOK=rbind(resOK,tx)
rm(t1)
gc(reset=T)

  #Endotherms
t1=fread("boostrap_slopTvel_endo_ELEV.csv",sep=";",dec=".",h=T)
gc(reset=T)
tx=data.frame(Tbas_moy=as.numeric(names(tapply(t1$slop_velT,t1$basT,mean))),moy=tapply(t1$slop_velT,t1$basT,mean),sd=tapply(t1$slop_velT,t1$basT,sd),q025=tapply(t1$slop_velT,t1$basT,quantile,probs=0.025),
              q975=tapply(t1$slop_velT,t1$basT,quantile,probs=0.975), pv.inf0=tapply(t1$slop_velT,t1$basT,test2,mu=0),pv.sup0=tapply(t1$slop_velT,t1$basT,test1,mu=0),pv.inf1=tapply(t1$slop_velT,t1$basT,test2,mu=1),pv.sup1=tapply(t1$slop_velT,t1$basT,test1,mu=1))
gc(reset=T)
tx$LifeForm="endo"
t1a=subset(tx,pv.inf1>0.05 & pv.sup1>0.05)
range(t1a$Tbas_moy) #velocity slope never equal to 1
resOK=rbind(resOK,tx)
rm(t1)
gc(reset=T)

write.table(resOK,"resTestSlopTvel_ELEV.csv",sep=";",dec=".",row=F)
#Significance, magniture and variation of the effect of climate warming velocity on species range shifts depending to baseline temperature conditions and species life form
  #Tbas_moy =  Baseline temperature
  #moy = mean effect of climate warming velocity computed from the 5000 bootstrap models
  #sd = standard deviation of the effect of effect of climate warming velocity computed from the 5000 bootstrap models
  #med = median effect of climate warming velocity computed from the 5000 bootstrap models
  #q025 = 2.5% quantile of the bootstrap coefficient distribution
  #q975 = 97.5% quantile of the bootstrap coefficient distribution
  #pv.sup0 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 0; H1= coefficient greater than 0)
  #pv.inf0 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 0; H1= coefficient lower than 0)
  #pv.sup1 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 1; H1= coefficient greater than 1)
  #pv.inf1 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 1; H1= coefficient lower than 1)
  #LifeForm= Species life forms (ecto= ectotherms; endo=endotherms; phanero=phanerogams; crypto=cryptogams)

#Latitudinal and marine range shifts
Class <- as.data.frame(table(rSdata$Class[rSdata$Gradient == "Latitudinal" & rSdata$Ecosystem=="Marine"]))
names(Class) <- c("Class", "Freq")
Class=Class[which(Class$Freq>30), ]       # Criteria: Class > 30 obs
Class
#           Class Freq
#  Actinopterygii  720
#        Bivalvia   71
#  Chondrichthyes   83
# Florideophyceae   70
#      Gastropoda  134
#    Malacostraca  123
#     Maxillopoda   55
#    Phaeophyceae   59
#      Polychaeta  143

# Selecting observation and variables to analyse range shifts
data <- rSdata[rSdata$Gradient == "Latitudinal" & rSdata$Ecosystem=="Marine", chosen_varlatM]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data <- na.omit(data)
dim(data) # 1403 obs

setwd(rep_out)
setwd("model_sp6000")
t1=read.csv2("summary_coeff.csv",sep=";",dec=".",h=T)
t1a=subset(t1,model=="latM")
t2=read.csv2("coeff_latM_5000.csv",sep=";",dec=".",h=T)

v2=seq(-20,50,by=0.1) #range of baseline temperature tested
v1=seq(-10,20,by=0.1) #range of climate warming velocity tested
v3=seq(0,1,by=0.01) #range of standardized HFI tested

#Estimation of the effect of climate warming velocity (unscaled)
  #Climate warming velocity parameters
    #The effect of climate warming velocity depends to baseline temperatures and standardized human footprint
a2=t1a$moy[t1a$var=="g_velT"] #direct effect of climate warming velocity
a4=t1a$moy[t1a$var=="g_basT:g_velT"] #interaction between climate warming velocity and baseline temperatures
a5=t1a$moy[t1a$var=="g_HFI:g_velT"] #interaction between climate warming velocity and standardized HFI

    #ecotherms
i0=t1a$moy[t1a$var=="(Intercept)"] 
res=matrix(NA,ncol=4,nrow=length(v3))
for(i in 1:length(v2)){
  print(i)
  for(z in 1:length(v3)){
    res[z,1]=v3[z]
    res[z,2]=a2+(a4*(v2[i]-mean(data$BaselineT))/(2*sd(data$BaselineT)))+(a5*(v3[z]-mean(data$HFI))/(2*sd(data$HFI)))
    p1=res[z,2]*(v1-mean(data$LatVeloT))/(2*sd(data$LatVeloT)) + i0
    lm1=lm(p1~v1)
    res[z,3]=lm1$coefficients[[2]]
    res[z,4]=lm1$coefficients[[1]]
  }
  if(i==1){
    resM=data.frame(res)
    resM=resM[,-2]
    names(resM)=c("HFI","slop_velT","int")
    resM$basT=v2[i]      
  }else{
    res1=data.frame(res)
    res1=res1[,-2]
    names(res1)=c("HFI","slop_velT","int")
    res1$basT=v2[i] 
    resM=rbind(resM,res1)
  }
}
resM$LifeForm="ecto"
write.table(resM,"slopVelT_LATM.csv",sep=";",dec=".",row=F)

  #cryptogams
i0=t1a$moy[t1a$var=="(Intercept)"]+t1a$moy[t1a$var=="LifeFormcrypto"]
res=matrix(NA,ncol=4,nrow=length(v3))
for(i in 1:length(v2)){
  print(i)
  for(z in 1:length(v3)){
    res[z,1]=v3[z]
    res[z,2]=a2+(a4*(v2[i]-mean(data$BaselineT))/(2*sd(data$BaselineT)))+(a5*(v3[z]-mean(data$HFI))/(2*sd(data$HFI)))
    p1=res[z,2]*(v1-mean(data$LatVeloT))/(2*sd(data$LatVeloT)) + i0
    lm1=lm(p1~v1)
    res[z,3]=lm1$coefficients[[2]]
    res[z,4]=lm1$coefficients[[1]]
  }
  if(i==1){
    resM=data.frame(res)
    resM=resM[,-2]
    names(resM)=c("HFI","slop_velT","int")
    resM$basT=v2[i]      
  }else{
    res1=data.frame(res)
    res1=res1[,-2]
    names(res1)=c("HFI","slop_velT","int")
    res1$basT=v2[i] 
    resM=rbind(resM,res1)
  }
}
resM$LifeForm="crypto"
write.table(resM,"slopVelT_LATM.csv",sep=";",dec=".",row=F,append=T,col.names=F)

#Test of the effect of climate warming velocity in each boostrap model
nB=unique(t2$nB)
a=1
v2=seq(-5,35,by=0.1) #range of baseline temperature tested
for(j in nB){ #take a lot of time (more than 24 hrs)
  print(j)
  t2a=subset(t2,nB==j)
  t2a$velT=grepl("g_velT",as.character(t2a$var))
  t2a=subset(t2a,velT==T)
  #Climate warming velocity parameters computed from the focal model
    #The effect of climate warming velocity depends to the baseline temperature conditions and standardized human footprint index
  a2=t2a$Estimate[t2a$var=="g_velT"] #direct effect of climate warming velocity
  a4=t2a$Estimate[t2a$var=="g_basT:g_velT"] #interaction between climate warming velocity and baseline temperatures
  a5=t2a$Estimate[t2a$var=="g_HFI:g_velT"] #interaction between climate warming velocity and standardized HFI
  
  #Intercept
  t2a=subset(t2,nB==j)
  i0=t2a$Estimate[t2a$var=="(Intercept)"]
  res=matrix(NA,ncol=3,nrow=length(v3))
  for(i in 1:length(v2)){
    for(z in 1:length(v3)){
      res[z,1]=v3[z]
      p2=a2+(a4*(v2[i]-mean(data$BaselineT))/(2*sd(data$BaselineT)))+(a5*(v3[z]-mean(data$HFI))/(2*sd(data$HFI)))
      p1=p2*(v1-mean(data$LatVeloT))/(2*sd(data$LatVeloT)) + i0
      lm1=lm(p1~v1)
      res[z,2]=lm1$coefficients[[2]]
      res[z,3]=lm1$coefficients[[1]]
    }
    if(i==1){
      resM=data.frame(res)
      names(resM)=c("HFI","slop_velT","int")
      resM$basT=v2[i]      
    }else{
      res1=data.frame(res)
      names(res1)=c("HFI","slop_velT","int")
      res1$basT=v2[i] 
      resM=rbind(resM,res1)
    }
  }
  resM$nB=j 
  
  if(a==1){
    resOK=resM
  }else{
    resOK=rbind(resOK,resM)
  }
  a=a+1
}
write.table(resOK,"boostrap_slopTvel_LATM.csv",sep=";",dec=".",row=F)
rm(resOK)
gc(reset=T)

#Statistics of the relationships for a range of combinations of baseline temperature conditions, standardized HFI and climate warming velocity
#bootstrap test: H0= coefficient equal to 0; H1= coefficient greater than 0
test1=function(x,mu=0){
  return(1-(length(x[x>mu])/length(x)))
}

#bootstrap test: H0= coefficient equal to 0; H1= coefficient lower than 0
test2=function(x,mu=0){
  return(1-(length(x[x<mu])/length(x)))
}

t1=fread("boostrap_slopTvel_LATM.csv",sep=";",dec=".",h=T)
gc(reset=T)
t1$v1=paste(t1$basT,"_",t1$HFI,sep="")
gc(reset=T)
tx=data.frame(v1=names(tapply(t1$slop_velT,t1$v1,mean)),moy=tapply(t1$slop_velT,t1$v1,mean),sd=tapply(t1$slop_velT,t1$v1,sd),q025=tapply(t1$slop_velT,t1$v1,quantile,probs=0.025),
              q975=tapply(t1$slop_velT,t1$v1,quantile,probs=0.975), pv.inf0=tapply(t1$slop_velT,t1$v1,test2,mu=0),pv.sup0=tapply(t1$slop_velT,t1$v1,test1,mu=0),pv.inf1=tapply(t1$slop_velT,t1$v1,test2,mu=1),pv.sup1=tapply(t1$slop_velT,t1$v1,test1,mu=1))
gc(reset=T)
n1=strsplit(as.character(tx$v1),"_")
n1=unlist(n1)
tx$Tbas_moy=as.numeric(n1[seq(1,length(n1)-1,by=2)])
tx$HFI=as.numeric(n1[seq(2,length(n1),by=2)])

t1a=subset(tx,pv.inf1>0.05 & pv.sup1>0.05)
range(t1a$Tbas_moy) #range of baseline temperatures for which species range shift ~ climate waring velocity
range(t1a$HFI) #range of standardized HFI for which species range shift ~ climate waring velocity

write.table(tx,"resTestSlopTvel_LATM.csv",sep=";",dec=".",row=F)
#Significance, magniture and variation of the effect of climate warming velocity on species range shifts depending to baseline temperature conditions and standardized HFI
  #v1 = combination of baseline temperatures and standardized HFI
  #moy = mean effect of climate warming velocity computed from the 5000 bootstrap models
  #sd = standard deviation of the effect of effect of climate warming velocity computed from the 5000 bootstrap models
  #med = median effect of climate warming velocity computed from the 5000 bootstrap models
  #q025 = 2.5% quantile of the bootstrap coefficient distribution
  #q975 = 97.5% quantile of the bootstrap coefficient distribution
  #pv.sup0 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 0; H1= coefficient greater than 0)
  #pv.inf0 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 0; H1= coefficient lower than 0)
  #pv.sup1 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 1; H1= coefficient greater than 1)
  #pv.inf1 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 1; H1= coefficient lower than 1)
  #Tbas_moy =  Baseline temperature
  #HFI= standardized human footprint index
rm(t1)
gc(reset=T)

#Latitudinal and terrestrial range shifts
Class <- as.data.frame(table(rSdata$Class[rSdata$Gradient == "Latitudinal" & rSdata$Ecosystem=="Terrestrial"]))
names(Class) <- c("Class", "Freq")
Class=Class[which(Class$Freq>30), ]       # Criteria: Class > 30 obs
Class
#          Class Freq
#       Amphibia  211
#      Arachnida  444
#           Aves 4882
#      Bryopsida  288
#  Equisetopsida   34
#        Insecta 4237
#     Liliopsida 1109
#  Magnoliopsida 3617
#   Malacostraca   37
#      Pinopsida   98
# Polypodiopsida  105
#       Reptilia   90

# Selecting observation and variables to analyse range shifts
data <- rSdata[rSdata$Gradient == "Latitudinal" & rSdata$Ecosystem=="Terrestrial", chosen_varlatT]
data <- data[which(data$Class %in% unique(Class$Class)), ]
data <- na.omit(data)
dim(data) # 15118 obs

setwd(rep_out)
setwd("model_sp6000")
t1=read.csv2("summary_coeff.csv",sep=";",dec=".",h=T)
t1a=subset(t1,model=="latT")
t2=read.csv2("coeff_latT_5000.csv",sep=";",dec=".",h=T)

v1=seq(-10,20,by=0.1) #range of climate warming velocity tested
v3=seq(0,1,by=0.01) #range of standardized HFI tested

#Estimation of the effect of climate warming velocity (unscaled)
  #Climate warming velocity parameters
    #The effect of climate warming velocity depends to the standardized human footprint index
a2=t1a$moy[t1a$var=="g_velT"] #direct effect of climate warming velocity
a5=t1a$moy[t1a$var=="g_HFI:g_velT"] #interaction between climate warming velocity and standardized HFI

res=matrix(NA,ncol=4,nrow=length(v3))
#ectotherms
i0=t1a$moy[t1a$var=="(Intercept)"]
for(z in 1:length(v3)){
  res[z,1]=v3[z]
  res[z,2]=a2+(a5*(v3[z]-mean(data$HFI))/(2*sd(data$HFI)))
  p1=res[z,2]*(v1-mean(data$LatVeloT))/(2*sd(data$LatVeloT)) + i0
  lm1=lm(p1~v1)
  res[z,3]=lm1$coefficients[[2]]
  res[z,4]=lm1$coefficients[[1]]
}
resM=data.frame(res)
resM=resM[,-2]
names(resM)=c("HFI","slop_velT","int")
resM$LifeForm="ecto"
write.table(resM,"slopVelT_LATT.csv",sep=";",dec=".",row=F)

#endotherms
i0=t1a$moy[t1a$var=="(Intercept)"]+t1a$moy[t1a$var=="(LifeFormendo)"]
for(z in 1:length(v3)){
  res[z,1]=v3[z]
  res[z,2]=a2+(a5*(v3[z]-mean(data$HFI))/(2*sd(data$HFI)))
  p1=res[z,2]*(v1-mean(data$LatVeloT))/(2*sd(data$LatVeloT)) + i0
  lm1=lm(p1~v1)
  res[z,3]=lm1$coefficients[[2]]
  res[z,4]=lm1$coefficients[[1]]
}
resM=data.frame(res)
resM=resM[,-2]
names(resM)=c("HFI","slop_velT","int")
resM$LifeForm="endo"
write.table(resM,"slopVelT_LATT.csv",sep=";",dec=".",row=F,col.names=F,append=T)

#phanerogams
i0=t1a$moy[t1a$var=="(Intercept)"]+t1a$moy[t1a$var=="(LifeFormphanero)"]
for(z in 1:length(v3)){
  res[z,1]=v3[z]
  res[z,2]=a2+(a5*(v3[z]-mean(data$HFI))/(2*sd(data$HFI)))
  p1=res[z,2]*(v1-mean(data$LatVeloT))/(2*sd(data$LatVeloT)) + i0
  lm1=lm(p1~v1)
  res[z,3]=lm1$coefficients[[2]]
  res[z,4]=lm1$coefficients[[1]]
}
resM=data.frame(res)
resM=resM[,-2]
names(resM)=c("HFI","slop_velT","int")
resM$LifeForm="phanero"
write.table(resM,"slopVelT_LATT.csv",sep=";",dec=".",row=F,col.names=F,append=T)

#cryptogams
i0=t1a$moy[t1a$var=="(Intercept)"]+t1a$moy[t1a$var=="(LifeFormcrypto)"]
for(z in 1:length(v3)){
  res[z,1]=v3[z]
  res[z,2]=a2+(a5*(v3[z]-mean(data$HFI))/(2*sd(data$HFI)))
  p1=res[z,2]*(v1-mean(data$LatVeloT))/(2*sd(data$LatVeloT)) + i0
  lm1=lm(p1~v1)
  res[z,3]=lm1$coefficients[[2]]
  res[z,4]=lm1$coefficients[[1]]
}
resM=data.frame(res)
resM=resM[,-2]
names(resM)=c("HFI","slop_velT","int")
resM$LifeForm="crypto"
write.table(resM,"slopVelT_LATT.csv",sep=";",dec=".",row=F,col.names=F,append=T)

#Test of the effect of climate warming velocity in each boostrap model
nB=unique(t2$nB)
a=1
for(j in nB){ #take some time
  print(j)
  t2a=subset(t2,nB==j)
  t2a$velT=grepl("g_velT",as.character(t2a$var))
  t2a=subset(t2a,velT==T)
  #Climate warming velocity parameters computed from the focal model
    #The effect of climate warming velocity depends to the standardized human footprint index
  a2=t2a$Estimate[t2a$var=="g_velT"] #direct effect of climate warming velocity
  a5=t2a$Estimate[t2a$var=="g_HFI:g_velT"] #interaction between climate warming velocity and standardized HFI
  
  #Intercept
  t2a=subset(t2,nB==j)
  i0=t2a$Estimate[t2a$var=="(Intercept)"]
  
  res=matrix(NA,ncol=4,nrow=length(v3))
  
  for(z in 1:length(v3)){
    res[z,1]=v3[z]
    res[z,2]=a2+(a5*(v3[z]-mean(data$HFI))/(2*sd(data$HFI)))
    p1=res[z,2]*(v1-mean(data$LatVeloT))/(2*sd(data$LatVeloT)) + i0
    lm1=lm(p1~v1)
    res[z,3]=lm1$coefficients[[2]]
    res[z,4]=lm1$coefficients[[1]]
  }
  resM=data.frame(res)
  resM=resM[,-2]
  names(resM)=c("HFI","slop_velT","int")
  resM$nB=j 
  
  if(a==1){
    resOK=resM
  }else{
    resOK=rbind(resOK,resM)
  }
  a=a+1
}
setwd(rep_out)
write.table(resOK,"boostrap_slopTvel_LATT.csv",sep=";",dec=".",row=F)
rm(resOK)
gc(reset=T)

#Statistics of the relationships for a range of combinations of standardized HFI and climate warming velocity
#bootstrap test: H0= coefficient equal to 0; H1= coefficient greater than 0
test1=function(x,mu=0){
  return(1-(length(x[x>mu])/length(x)))
}

#bootstrap test: H0= coefficient equal to 0; H1= coefficient lower than 0
test2=function(x,mu=0){
  return(1-(length(x[x<mu])/length(x)))
}

t1=fread("boostrap_slopTvel_LATT.csv",sep=";",dec=".",h=T)
gc(reset=T)
tx=data.frame(HFI=as.numeric(names(tapply(t1$slop_velT,t1$HFI,mean))),moy=tapply(t1$slop_velT,t1$HFI,mean),sd=tapply(t1$slop_velT,t1$HFI,sd),q025=tapply(t1$slop_velT,t1$HFI,quantile,probs=0.025),
              q975=tapply(t1$slop_velT,t1$HFI,quantile,probs=0.975), pv.inf0=tapply(t1$slop_velT,t1$HFI,test2,mu=0),pv.sup0=tapply(t1$slop_velT,t1$HFI,test1,mu=0),pv.inf1=tapply(t1$slop_velT,t1$HFI,test2,mu=1),pv.sup1=tapply(t1$slop_velT,t1$HFI,test1,mu=1))
gc(reset=T)
t1a=subset(tx,pv.inf1>0.05 & pv.sup1>0.05)
range(t1a$HFI) #velocity slope never equals to 1

write.table(tx,"resTestSlopTvel_LATT.csv",sep=";",dec=".",row=F)
#Significance, magniture and variation of the effect of climate warming velocity on species range shifts depending to baseline temperature conditions and standardized HFI
  #HFI= standardized human footprint index
  #moy = mean effect of climate warming velocity computed from the 5000 bootstrap models
  #sd = standard deviation of the effect of effect of climate warming velocity computed from the 5000 bootstrap models
  #med = median effect of climate warming velocity computed from the 5000 bootstrap models
  #q025 = 2.5% quantile of the bootstrap coefficient distribution
  #q975 = 97.5% quantile of the bootstrap coefficient distribution
  #pv.sup0 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 0; H1= coefficient greater than 0)
  #pv.inf0 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 0; H1= coefficient lower than 0)
  #pv.sup1 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 1; H1= coefficient greater than 1)
  #pv.inf1 = p-value of the bootstrap test testing the significance of the effect (H0= coefficient equal to 1; H1= coefficient lower than 1)
rm(t1)