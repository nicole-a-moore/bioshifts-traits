#AIM= Inferring range shift at the taxonimic class level
#based on script used for the analysis published in Lenoir et al. 2020 and available here: https://figshare.com/articles/dataset/BioShifts_a_global_geodatabase_of_climate-induced_species_redistribution_over_land_and_sea/7413365?file=22057815

#required packages
#library(nlme) #version 3.1.140
library(MuMIn) #version 1.43.6
library(lme4) #version 1.1.21
library(optimx) #version 2018-7.10

#repository/file location
rep_data="XXXX" #put here the path towards the repository where Table_S1.csv is located.
#rep_out="......." #inform the location where all results have to be saved


#############################
#### DATASET preparation ####
#############################
setwd(rep_data)
rSdata = read.table("Table_S1.csv",sep=";",h=T,dec=".",stringsAsFactors = FALSE) #it's the data used in Lenoir et al. 2020; so it's a subset of Bioshift v1 available here: https://figshare.com/articles/dataset/BioShifts_a_global_geodatabase_of_climate-induced_species_redistribution_over_land_and_sea/7413365?file=22057815
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

rSdata = rSdata[order(rSdata$n_cl),]

#transforming continuous method variables in qualitative variables => two benefits: (1) all variables area qualitative (2) better than quantitative in case the effect is not linear
q1=quantile(rSdata$Start,probs=c(0,0.25,0.5,0.75,1))
rSdata$StartF=cut(rSdata$Start,breaks=q1,include.lowest=T)
q1=quantile(rSdata$Area,probs=c(0,0.25,0.5,0.75,1))
rSdata$AreaF=cut(rSdata$Area,breaks=q1,include.lowest=T)
q1=quantile(rSdata$Ntaxa,probs=c(0,0.25,0.5,0.75,1))
rSdata$NtaxaF=cut(rSdata$Ntaxa,breaks=q1,include.lowest=T)

rSdata$Sampling = ifelse(rSdata$Sampling %in% c("IRR","MULT"),"MULT", rSdata$Sampling)

rSdata$IDn=1:nrow(rSdata)

#data subsetting: Gradient x Position x Ecosystem x Hemisphere
#Core of the distribution (O)
rSO.LT_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Terrestrial")# & rSdata$Hemisphere=="South")
#rSO.LTN_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Terrestrial")# & rSdata$Hemisphere=="North")
rSO.LM_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Marine")# & rSdata$Hemisphere=="South")
#rSO.LMN_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Marine")# & rSdata$Hemisphere=="North")
rSO.ET_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Elevation"& rSdata$Ecosystem=="Terrestrial")# & rSdata$Hemisphere=="South")
#rSO.ETN_posit <- which(rSdata$Position=="Centroid"& rSdata$Gradient=="Elevation"& rSdata$Ecosystem=="Terrestrial")# & rSdata$Hemisphere=="North")

#Margins (LE & TE)
rSMrg.LT_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Terrestrial")# & rSdata$Hemisphere=="South")
#rSMrg.LTN_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Terrestrial")# & rSdata$Hemisphere=="North")
rSMrg.LM_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Marine")# & rSdata$Hemisphere=="South")
#rSMrg.LMN_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Latitudinal"& rSdata$Ecosystem=="Marine")# & rSdata$Hemisphere=="North")
rSMrg.ET_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Elevation"& rSdata$Ecosystem=="Terrestrial")# & rSdata$Hemisphere=="South")
#rSMrg.ETN_posit <- which((rSdata$Position=="Trailing edge"|rSdata$Position=="Leading edge")& rSdata$Gradient=="Elevation"& rSdata$Ecosystem=="Terrestrial")# & rSdata$Hemisphere=="North")

# variable pre-selection (one by gradient)
chosen_varlat=c("ShiftR", "Position", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif","Source","Class","Family","Genus","Species","Start","Area","Ntaxa","IDn")
chosen_varele=c("ShiftR", "Position", "NtaxaF", "StartF", "AreaF", "PrAb", "Sampling", "Grain", "Quality", "Signif","Source","Class","Family","Genus","Species","Start","Area","Ntaxa","IDn")

#####################################################
#### MODELING: random effect structure selection ####
#####################################################
# I decided to fit one model by combination of position (Trailing edge or Leading edge or center), gradient (latitudinal or elevation) and realms (terrestrial or marine)
# An example for: Trailing edge x LATITUDINAL RANGE SHITS x TERRESTRIAL ECOSYSTEM
# teLTN hereafter
# Disparities between Classes and data selection
data <- rSdata[rSMrg.LT_posit,chosen_varlat] #create a data.frame with the data of realm, position and gradient you want to analyse
data=subset(data,Position=="Trailing edge") #here select the position you want
Class <- as.data.frame(table(data$Class))
names(Class) <- c("Class", "Freq")
Class$Shift <- c(tapply(data$ShiftR, data$Class, mean))
Class=Class[which(Class$Freq>10), ] #here I use a criteria of number of observation by class in order to have sufficient replicate among classes (I set it to 10 but you can use another threshold)
Class
Class
#         Class  Freq      Shift
#1          Aves  295 -0.7262751
#4       Insecta  246  2.3350293
#5    Liliopsida  132 -0.1617514
#7 Magnoliopsida  584 -0.2034522
#9      Pinopsida   18 -0.2047738
#12      Reptilia   18 -3.5604597


# Selecting observation and variables to analyse range shifts
# Criteria: Class > 10 obs
data <- data[which(data$Class %in% unique(Class$Class)), ]
data$Class <- as.factor(as.character(data$Class))
data$Family <- as.factor(as.character(data$Family))
data$Genus <- as.factor(as.character(data$Genus))
data <- na.omit(data)
data=droplevels(data)
dim(data) # 7359 obs x 15 vars


#selection of the set of qualitative variables to test as random effect (based on number of observations, and correlation among variables)table(data$PrAb)
table(data$PrAb) 
table(data$Sampling)
table(data$Grain)
table(data$Quality)
table(data$Signif) 
table(data$AreaF)
table(data$StartF)
table(data$NtaxaF)

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
      if(isSingular(mod1)==T){ #in case of singularity, a more simple random model structure is tested: (1|Genus)
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
#statistics of the selected model
R2m=r.squaredGLMM(mod1)[[1]] #marginal R2, ie part of variation explain by fixed effect (so here it's methodological factor)
R2c=r.squaredGLMM(mod1)[[2]] #conditional R2, ie part of variation explain by fixed effect and random effect (so here it's methodological factor)
aic=AIC(mod1)[[1]]
aicc=AICc(mod1)[[1]] #correction AIC usefull for low sampling size

#comparing linear model and LMM coeff estimation
mod0=lm(ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF, data)
summary(mod0)
plot(summary(mod1)$coeff[,1]~mod0$coeff, xlab="coeff. in LM", ylab="coeff. in LMM")
abline(a=0,b=1,col=2)

#here I consider range shift estimation corrected by methodological effect as the residuals of the selected LMM (ie mod0)
#extraction of the resisuals
data$shiftRcor=summary(mod1)$residuals
plot(shiftRcor~ShiftR,data)
abline(a=0,b=1,col=2)
lm0=lm(shiftRcor~ShiftR,data)
summary(lm0)
hist(data$ShiftR,col=rgb(1,0,0,0.5))
hist(data$shiftRcor,col=rgb(0,0,1,0.5),add=T)
#so by using residuals we decrease the variation in the raw range shift obs, it's like we work on a more homogenous variables as all the variations due to methods variation has been substracted to the shift.

x1=tapply(data$ShiftR,data$Class,mean)
x2=tapply(data$shiftRcor,data$Class,mean)#lower mean value than true obs 
plot(x1~x2)
lm0=lm(x1~x2,data)
summary(lm0) #it changes many things

x1=tapply(data$ShiftR,data$Class,var)
x2=tapply(data$shiftRcor,data$Class,var)#lower variance value than in true obs 
plot(x1~x2)
lm0=lm(x1~x2,data)
summary(lm0) #but high positive correlation meaning that relative variation is conserved

x1=tapply(data$ShiftR,data$Species,mean)
x2=tapply(data$shiftRcor,data$Species,mean)#lower mean value than in true obs 
plot(x1~x2)
lm0=lm(x1~x2,data)
summary(lm0) #but at the species level, we observe a high positive relationship so the corrected shifts did not change the pattern of range shift=> good

#looking at the relationship with climate velocity
#WARNING:here you need to adapt the code: if you look at the elevational range shift so you have to use EleVeloT in order to use the corresponding climate velocity
data2=rSdata[data$IDn,]
data2$shiftRcor=data$shiftRcor
par(mfrow=c(1,2))
plot(ShiftR~LatVeloT,data2)
lm0=lm(ShiftR~LatVeloT+I(LatVeloT^2),data2)
summary(lm0)
x1=seq(min(data2$LatVeloT),max(data2$LatVeloT),le=100)
p1=predict(lm0,newdata=data.frame(LatVeloT=x1),type="response")
points(p1~x1,type="l",col=2,lwd=2)
plot(shiftRcor~LatVeloT,data2)#the pattern is changing a lot
lm0=lm(shiftRcor~LatVeloT+I(LatVeloT^2),data2)
summary(lm0)
x1=seq(min(data2$LatVeloT),max(data2$LatVeloT),le=100)
p1=predict(lm0,newdata=data.frame(LatVeloT=x1),type="response")
points(p1~x1,type="l",col=2,lwd=2)
#we observed a relationship for the raw range shift (an unimodal relationship with higher shifts when absolute climate velocity increase).
#for the range shift corrected by method, no relationship with climatic velocity is observed

#looking at the relationship with lag
data2$lag=data2$LatVeloT-data2$ShiftR
par(mfrow=c(1,2))
plot(lag~LatVeloT,data2,main="raw lag estimate")
data2$si=sign(abs(data2$LatVeloT)-abs(data2$ShiftR))
data2$si[sign(data2$LatVeloT)==1 & sign(data2$ShiftR)==(-1)]=1
data2$si[sign(data2$LatVeloT)==(-1) & sign(data2$ShiftR)==1]=1
data2$lag2=data2$si*abs(data2$lag)#positive values are a lag (range shift lower smalled than expected or in opposite directions) and negative values are range shift larger than expected
plot(lag2~LatVeloT,data2,main="raw lag estimate with correction for special case")
#change the shape in extreme negative climatic shifts

par(mfrow=c(1,2))
data2$lagC=data2$LatVeloT-data2$shiftRcor
plot(lagC~LatVeloT,data2,main="corrected lag estimate")
data2$si=sign(abs(data2$LatVeloT)-abs(data2$shiftRcor))
data2$si[sign(data2$LatVeloT)==1 & sign(data2$shiftRcor)==(-1)]=1
data2$si[sign(data2$LatVeloT)==(-1) & sign(data2$shiftRcor)==1]=1
data2$lagC2=data2$si*abs(data2$lagC)#positive values are a lag and negative values are range shift larger than expected
plot(lagC2~LatVeloT,data2,main="corrected lag estimate with correction for special case")
#change more things than the above lag metrics, still the highly change are observed at extreme negative values


plot(lag2~lagC2,data2)
lm0=lm(lag~lagC,data2)
summary(lm0) #we observe a high positive relationship so the corrected lag did not change the pattern of the observed lag=> good

par(mfrow=c(1,2))
plot(lag2~LatVeloT,data2)
lm0=lm(lag2~LatVeloT+I(LatVeloT^2),data2)
summary(lm0) #significant; R2=26%
x1=seq(min(data2$LatVeloT),max(data2$LatVeloT),le=100)
p1=predict(lm0,newdata=data.frame(LatVeloT=x1),type="response")
points(p1~x1,type="l",col=2,lwd=2) #bimodal relationship

plot(lagC2~LatVeloT,data2) 
lm0=lm(lagC2~LatVeloT+I(LatVeloT^2),data2)
summary(lm0)#significant; R2=51.8
x1=seq(min(data2$LatVeloT),max(data2$LatVeloT),le=100)
p1=predict(lm0,newdata=data.frame(LatVeloT=x1),type="response")
points(p1~x1,type="l",col=2,lwd=2) #bimodal relationship

#final comments
##two future direction for the analysis
#1#conducting genetic diversity analysis on the lag correcting for special case only (variable called lag2 above) => could change output of the analysis as special case in this example are observed for 34 observations
#2#conducting genetic diversity analysis on the lag both correcting for special case only and accounting for the method of range shift estimation in studies (varaible called lagC2 above) => could change considerably the outputs of the analysis

#The script can be used for every data you want infer a corrected range shift. You simply need to change the criteria of data selection at lines 63-64 of this script.
