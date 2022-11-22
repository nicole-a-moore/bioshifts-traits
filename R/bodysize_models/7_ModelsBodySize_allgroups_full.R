bod<-read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/CompilationBodySizeBioshiftsv1harmonized.csv")
trait <- read.csv('G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Bioshiftv1_checked.csv')
trait$spnameoriginal <- trait$new_name
trait$SpeciesChecked[is.na(trait$SpeciesChecked)] <- trait$new_name[is.na(trait$SpeciesChecked)] #a few ID problem


#----------------------
#merge different classes of body sizes
bod$UnitB <- bod$Code
bod$Code[bod$Code=="VegetativeHeight"] <- 'BodyLengthMean' #use vegetative length for plants
bod$Code[bod$Code=="WingspanMean"] <- 'BodyLengthMean' #use wingspan for butterflies

#--------------------------
#use insects orders to avoid mixing butterflies with others
bod$class[bod$class=="Insecta"] <- bod$order[bod$class=="Insecta"]
trait$class[trait$class=="Insecta"] <- trait$order[trait$class=="Insecta"]

#----------------------
#correct and reorganize methodological variables
trait$data[trait$data=="occurence-based"] <- "occurrence-based"
trait$sampling <- ifelse(trait$sampling == "TWO","TWO","MULTIPLE")
trait$grain_size <- ifelse(trait$grain_size %in% c("large","very_large"),"large",trait$grain_size)
trait$uncertainty_distribution <- ifelse(trait$uncertainty_distribution %in% c("RESAMPLING","RESAMPLING(same)"),"RESAMPLING",
                                         ifelse(trait$uncertainty_distribution %in% c("MODEL","MODEL+RESAMPLING(same)","RESAMPLING+MODEL"),"MODEL",
                                                ifelse(trait$uncertainty_distribution %in% c("DETECTABILITY","RESAMPLING(same)+DETECTABILITY"),"DETECTABILITY",
                                                       trait$uncertainty_distribution)))
trait$uncertainty_distribution <- ifelse(trait$uncertainty_distribution == "RAW","OPPORTUNISTIC","PROCESSED")  

#transform study area
trait$id.area <- log(trait$id.area)

#----------------------
#remove freshwater fishes
trait <- trait[-which(trait$class == "Actinopterygii" & trait$eco=="T"),]

#remove marine birds
trait <- trait[-which(trait$class == "Aves" & trait$eco=="M"),]

#characteristics dataset
nlevels(factor(trait$id[trait$type=="LAT" & trait$param != "O"]))#113
nlevels(factor(trait$id[trait$type=="LAT" & trait$param != "O" & trait$eco=="M"]))#36
nlevels(factor(trait$SpeciesChecked[trait$type=="LAT" & trait$param != "O"]))#4094
length((trait$SpeciesChecked[trait$type=="LAT" & trait$param != "O"]))#7678

#----------------------
#transform body size
BodyLengthMean <- bod[bod$Code == "BodyLengthMean",]
BodyLengthMean$BodyLength <- as.numeric(as.character(BodyLengthMean$BodyLength))
BodyLengthMean <- BodyLengthMean[which(is.na(BodyLengthMean$BodyLength)==F),]
hist(BodyLengthMean$BodyLength)
hist(BodyLengthMean$BodyLength[BodyLengthMean$kingdom=='Plantae'])
hist(BodyLengthMean$BodyLength[BodyLengthMean$kingdom=='Animalia'])

BodyLengthMean$BodyLength <- log(BodyLengthMean$BodyLength)
hist(BodyLengthMean$BodyLength)

table(BodyLengthMean$Database)

#----------------------
#prepare datasets for each gradient

#elevation
ele <- trait[trait$type == "ELE",]
ele <- ele[which(is.na(ele$v.ele.mean)==F),]
ele$lags <- ele$shift - ele$v.ele.mean
hist(ele$lags)


#latitude
lat <- trait[trait$type == "LAT",]
lat$lags <- lat$shift - lat$v.lat.mean
hist(lat$lags)

#add body size
ele$BodyLengthMean <- BodyLengthMean$BodyLength[match(ele$SpeciesChecked,BodyLengthMean$SpeciesChecked)]
lat$BodyLengthMean <- BodyLengthMean$BodyLength[match(lat$SpeciesChecked,BodyLengthMean$SpeciesChecked)]

ele$UnitB <- BodyLengthMean$UnitB[match(ele$SpeciesChecked,BodyLengthMean$SpeciesChecked)]
lat$UnitB <- BodyLengthMean$UnitB[match(lat$SpeciesChecked,BodyLengthMean$SpeciesChecked)]

# #remove optimum
# lat <- lat[lat$param %in% c("TE","LE"),]
# ele <- ele[ele$param %in% c("TE","LE"),]
# 
# lat$compiltrait <- ifelse(is.na(lat$BodyLengthMean)==F,"found","no")
# table(lat$compiltrait)[1]/sum(table(lat$compiltrait))
# aa <- table(lat$compiltrait,lat$eco)
# barplot(t(round(t(aa)*100/apply(aa,2,sum))),legend=T)
# aa <- table(lat$compiltrait,lat$kingdom)
# barplot(t(round(t(aa)*100/apply(aa,2,sum))),legend=F,axes=F)
# axis(2,las=2,tck=-0.02)
# 
# aa <- table(lat$compiltrait,lat$class)
# aa <- aa[,-which(aa[1,] < 15)]
# par(mar=c(8,2,1,1))
# barplot(t(round(t(aa)*100/apply(aa,2,sum))),legend=F,axes=F,las=2,col=c(grey(0.2),"white"))
# axis(2,las=2,tck=-0.02)
# 
# t(round(t(aa)*100/apply(aa,2,sum)))
#----------------------------------
#prepare datasets
#----------------------------------
##select only classes with >10 species
dat.MeanBodySize.ele <- data.frame(ele[which(is.na(ele$BodyLengthMean)==F),])
dat.MeanBodySize.ele <- dat.MeanBodySize.ele[which(is.na(dat.MeanBodySize.ele$v.ele.mean)==F),]

#remove optimum
dat.MeanBodySize.ele <- dat.MeanBodySize.ele[dat.MeanBodySize.ele$param %in% c("TE","LE"),]

dsp <- unique(dat.MeanBodySize.ele[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl > 15))
dat.MeanBodySize.ele <- dat.MeanBodySize.ele[which(dat.MeanBodySize.ele$class %in% nn),]
dat.MeanBodySize.ele <-droplevels(dat.MeanBodySize.ele)

#scale variables
dat.MeanBodySize.ele$BodyLengthMean <- scale(dat.MeanBodySize.ele$BodyLengthMean)
dat.MeanBodySize.ele$n <- scale(dat.MeanBodySize.ele$n)
dat.MeanBodySize.ele$id.area <- scale(dat.MeanBodySize.ele$id.area)
dat.MeanBodySize.ele$start <- scale(dat.MeanBodySize.ele$start)
dat.MeanBodySize.ele$dur <- scale(dat.MeanBodySize.ele$dur)
dat.MeanBodySize.ele$shift <- scale(dat.MeanBodySize.ele$shift)
dat.MeanBodySize.ele$v.ele.mean <- scale(dat.MeanBodySize.ele$v.ele.mean)

#----------------------------------
##select only classes with >10 species
dat.MeanBodySize.lat <- data.frame(lat[which(is.na(lat$BodyLengthMean)==F),])

#characteristics dataset
nlevels(factor(dat.MeanBodySize.lat$id))#68
nlevels(factor(dat.MeanBodySize.lat$SpeciesChecked))#1697
length((dat.MeanBodySize.lat$SpeciesChecked))#3690

dsp <- unique(dat.MeanBodySize.lat[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >15))

#remove optimum
dat.MeanBodySize.lat <- dat.MeanBodySize.lat[dat.MeanBodySize.lat$param %in% c("TE","LE"),]
dat.MeanBodySize.lat <- dat.MeanBodySize.lat[which(dat.MeanBodySize.lat$class %in% nn),]
dat.MeanBodySize.lat <-droplevels(dat.MeanBodySize.lat)

#scale variables
dat.MeanBodySize.lat$BodyLengthMean <- scale(dat.MeanBodySize.lat$BodyLengthMean)
dat.MeanBodySize.lat$n <- scale(dat.MeanBodySize.lat$n)
dat.MeanBodySize.lat$id.area <- scale(dat.MeanBodySize.lat$id.area)
dat.MeanBodySize.lat$start <- scale(dat.MeanBodySize.lat$start)
dat.MeanBodySize.lat$dur <- scale(dat.MeanBodySize.lat$dur)
dat.MeanBodySize.lat$shift <- scale(dat.MeanBodySize.lat$shift)
dat.MeanBodySize.lat$v.lat.mean <- scale(dat.MeanBodySize.lat$v.lat.mean)

#characteristics dataset
table(dat.MeanBodySize.lat$param)
length(dat.MeanBodySize.lat$param)

#----------------------------------
#per parameter
#-----------------------------------
#----------------------------------
ele_te <- ele[ele$param=="TE",]
ele_le <- ele[ele$param=="LE",]
ele_o <- ele[ele$param=="O",]
lat_te <- lat[lat$param=="TE",]
lat_le <- lat[lat$param=="LE",]
lat_o <- lat[lat$param=="O",]

#----------------------------------
#check variables
hist(lat$v.lat.mean)#negative velocities!!
hist(lat_le$v.lat.mean,add=T,col="red")
hist(lat_te$v.lat.mean,add=T,col="blue")

# #remove negative velocities?
# lat_le <- lat_le[-which(lat_le$v.lat.mean <0),]
# lat_te <- lat_te[-which(lat_te$v.lat.mean <0),]
# lat_o <- lat_o[-which(lat_o$v.lat.mean <0),]
# 
# ele_le <- ele_le[-which(ele_le$v.lat.mean <0),]
# ele_te <- ele_te[-which(ele_te$v.lat.mean <0),]
# ele_o <- ele_o[-which(ele_o$v.lat.mean <0),]

#-----------------------------------
lat_le <- lat_le[which(is.na(lat_le$BodyLengthMean)==F),]
dsp <- unique(lat_le[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
lat_le <- lat_le[which(lat_le$class %in% nn),]
lat_le <-droplevels(lat_le)

#scale variables
lat_le$BodyLengthMean <- scale(lat_le$BodyLengthMean)
lat_le$n <- scale(lat_le$n)
lat_le$id.area <- scale(lat_le$id.area)
lat_le$start <- scale(lat_le$start)
lat_le$dur <- scale(lat_le$dur)
lat_le$shift <- scale(lat_le$shift)
lat_le$v.lat.mean <- scale(lat_le$v.lat.mean)

#----------------------------------
lat_te <- lat_te[which(is.na(lat_te$BodyLengthMean)==F),]
dsp <- unique(lat_te[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
lat_te <- lat_te[which(lat_te$class %in% nn),]
lat_te <-droplevels(lat_te)

#scale variables
lat_te$BodyLengthMean <- scale(lat_te$BodyLengthMean)
lat_te$n <- scale(lat_te$n)
lat_te$id.area <- scale(lat_te$id.area)
lat_te$start <- scale(lat_te$start)
lat_te$dur <- scale(lat_te$dur)
lat_te$shift <- scale(lat_te$shift)
lat_te$v.lat.mean <- scale(lat_te$v.lat.mean)

#----------------------------------
lat_o <- lat_o[which(is.na(lat_o$BodyLengthMean)==F),]
dsp <- unique(lat_o[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
lat_o <- lat_o[which(lat_o$class %in% nn),]
lat_o <-droplevels(lat_o)

#scale variables
lat_o$BodyLengthMean <- scale(lat_o$BodyLengthMean)
lat_o$n <- scale(lat_o$n)
lat_o$id.area <- scale(lat_o$id.area)
lat_o$start <- scale(lat_o$start)
lat_o$dur <- scale(lat_o$dur)
lat_o$shift <- scale(lat_o$shift)
lat_o$v.lat.mean <- scale(lat_o$v.lat.mean)

#----------------------------------
ele_te <- ele_te[which(is.na(ele_te$BodyLengthMean)==F),]
dsp <- unique(ele_te[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
ele_te <- ele_te[which(ele_te$class %in% nn),]
ele_te <-droplevels(ele_te)

#scale variables
ele_te$BodyLengthMean <- scale(ele_te$BodyLengthMean)
ele_te$n <- scale(ele_te$n)
ele_te$id.area <- scale(ele_te$id.area)
ele_te$start <- scale(ele_te$start)
ele_te$dur <- scale(ele_te$dur)
ele_te$shift <- scale(ele_te$shift)
ele_te$v.ele.mean <- scale(ele_te$v.ele.mean)

#----------------------------------
ele_le <- ele_le[which(is.na(ele_le$BodyLengthMean)==F),]
dsp <- unique(ele_le[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
ele_le <- ele_le[which(ele_le$class %in% nn),]
ele_le <-droplevels(ele_le)

#scale variables
ele_le$BodyLengthMean <- scale(ele_le$BodyLengthMean)
ele_le$n <- scale(ele_le$n)
ele_le$id.area <- scale(ele_le$id.area)
ele_le$start <- scale(ele_le$start)
ele_le$dur <- scale(ele_le$dur)
ele_le$shift <- scale(ele_le$shift)
ele_le$v.ele.mean <- scale(ele_le$v.ele.mean)

#----------------------------------
ele_o <- ele_o[which(is.na(ele_o$BodyLengthMean)==F),]
dsp <- unique(ele_o[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
ele_o <- ele_o[which(ele_o$class %in% nn),]
ele_o <-droplevels(ele_o)

#scale variables
ele_o$BodyLengthMean <- scale(ele_o$BodyLengthMean)
ele_o$n <- scale(ele_o$n)
ele_o$id.area <- scale(ele_o$id.area)
ele_o$start <- scale(ele_o$start)
ele_o$dur <- scale(ele_o$dur)
ele_o$shift <- scale(ele_o$shift)
ele_o$v.ele.mean <- scale(ele_o$v.ele.mean)
# 
# #-------------------------------------------
#models
#------------------------------------------
library(lme4);library(lmerTest); library(corrplot); library(MuMIn)
#--------------------------

#---------------------------------------
#MODELS with methods only
#--------------------------------------
meth <- trait[,names(trait) %in% c("param","data","sampling","grain_size","uncertainty_distribution","uncertainty_parameter","start","duration","id.area","dur")] #remove "n" as can capture diversity gradients
meth <- na.omit(meth)
meth[,names(meth) %in% c("param","data","sampling","grain_size","uncertainty_distribution","uncertainty_parameter")] <- apply(meth[,names(meth) %in% c("param","data","sampling","grain_size","uncertainty_distribution","uncertainty_parameter")],2,function(x) as.numeric(factor(x)))
cormat <- cor(meth)
corrplot(cormat)
#id.area and grain size are correlated so pick only area
#start and duration are correlated so pick only dur

##elevation
table(dat.MeanBodySize.ele$param)
table(dat.MeanBodySize.ele$data)
table(dat.MeanBodySize.ele$sampling)
table(dat.MeanBodySize.ele$grain_size)
table(dat.MeanBodySize.ele$uncertainty_distribution)
table(dat.MeanBodySize.ele$uncertainty_parameter)
table(dat.MeanBodySize.ele$article_id)

par(mfrow=c(2,2))
hist(dat.MeanBodySize.ele$n)
hist(dat.MeanBodySize.ele$start)
hist(dat.MeanBodySize.ele$dur)
hist(dat.MeanBodySize.ele$id.area)

mod.meth.ele <- lmer(shift ~  param + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + 
                       (1|article_id) + (1|class),data=dat.MeanBodySize.ele,REML=TRUE,na.action="na.fail",
                     control=lmerControl(optimizer="Nelder_Mead",
                                         optCtrl=list(maxfun=1e4)))
summary(mod.meth.ele)

R2.MeanBodySize.ele<-r.squaredGLMM(mod.meth.ele)
R2.MeanBodySize.ele


##latitude
table(dat.MeanBodySize.lat$param)
table(dat.MeanBodySize.lat$data)
table(dat.MeanBodySize.lat$sampling)
table(dat.MeanBodySize.lat$grain_size)
table(dat.MeanBodySize.lat$uncertainty_distribution)
table(dat.MeanBodySize.lat$uncertainty_parameter)
table(dat.MeanBodySize.lat$article_id)
table(dat.MeanBodySize.lat$SpeciesChecked)#not enough repetitions

par(mfrow=c(2,2))
hist(dat.MeanBodySize.lat$n)
hist(dat.MeanBodySize.lat$start)
hist(dat.MeanBodySize.lat$dur)
hist(dat.MeanBodySize.lat$id.area)
mod.meth.lat <- lmer(shift ~ param + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (1|article_id) + (1|class),data=dat.MeanBodySize.lat,REML=TRUE,na.action="na.fail")
summary(mod.meth.lat)

R2.MeanBodySize.lat<-r.squaredGLMM(mod.meth.lat)
R2.MeanBodySize.lat

#--------------------------------------------
##full models
#-------------------------------------------
# #elevation
par(mfrow=c(1,3))
hist(dat.MeanBodySize.ele$shift)
hist(dat.MeanBodySize.ele$v.ele.mean)
hist(dat.MeanBodySize.ele$BodyLengthMean)

table(dat.MeanBodySize.ele$param,dat.MeanBodySize.ele$class)
mod.full.ele <- lmer(shift ~ v.ele.mean * BodyLengthMean +  dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (1|article_id)+(v.ele.mean * BodyLengthMean|param/class),data=dat.MeanBodySize.ele,REML=TRUE,na.action="na.fail")
summary(mod.full.ele)
rr <- ranef(mod.full.ele, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.ele <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.ele <- confint(mod.full.ele, method="Wald")
cimod.full.ele <- cbind(fixef(mod.full.ele),na.omit(cimod.full.ele))
R2.mod.full.ele<-r.squaredGLMM(mod.full.ele)
R2.mod.full.ele
# 
# mod.full.ele2 <- lmer(shift ~ v.ele.mean * BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter+uncertainty_distribution+(1|article_id)+(v.ele.mean * BodyLengthMean|param/class),data=dat.MeanBodySize.ele,REML=TRUE,na.action="na.fail")
# summary(mod.full.ele2)
# rr <- ranef(mod.full.ele2, condVar=TRUE)
# dd <- as.data.frame(rr)
# ddmod.full.ele2 <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
# lattice::dotplot(rr)
#  cimod.full.ele2 <- confint(mod.full.ele2, method="Wald")
#  cimod.full.ele2 <- cbind(fixef(mod.full.ele2),na.omit(cimod.full.ele2))
# R2.mod.full.ele2<-r.squaredGLMM(mod.full.ele2)
# R2.mod.full.ele2


mod.full.ele.TE <- lmer(shift ~ v.ele.mean * BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (v.ele.mean * BodyLengthMean|class)+(1|article_id),data=ele_te,REML=TRUE,na.action="na.fail")
summary(mod.full.ele.TE)
rr <- ranef(mod.full.ele.TE, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.ele.TE <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.ele.TE <- confint(mod.full.ele.TE, method="Wald")
cimod.full.ele.TE <- cbind(fixef(mod.full.ele.TE),na.omit(cimod.full.ele.TE))
lattice::dotplot(cimod.full.ele.TE)
R2.mod.full.ele.TE<-r.squaredGLMM(mod.full.ele.TE)
R2.mod.full.ele.TE

mod.full.ele.LE <- lmer(shift ~ v.ele.mean * BodyLengthMean + (v.ele.mean * BodyLengthMean|class)+(1|article_id),data=ele_le,REML=TRUE,na.action="na.fail")
summary(mod.full.ele.LE)
rr <- ranef(mod.full.ele.LE, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.ele.LE <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.ele.LE <- confint(mod.full.ele.LE, method="Wald")
cimod.full.ele.LE <- cbind(fixef(mod.full.ele.LE),na.omit(cimod.full.ele.LE))
lattice::dotplot(cimod.full.ele.LE)
R2.mod.full.ele.LE<-r.squaredGLMM(mod.full.ele.LE)
R2.mod.full.ele.LE

mod.full.ele.O <- lmer(shift ~ v.ele.mean * BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (v.ele.mean * BodyLengthMean|class)+(1|article_id),data=ele_o,REML=TRUE,na.action="na.fail")
summary(mod.full.ele.O)
rr <- ranef(mod.full.ele.O, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.ele.O <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.ele.O <- confint(mod.full.ele.O, method="Wald")
cimod.full.ele.O <- cbind(fixef(mod.full.ele.O),na.omit(cimod.full.ele.O))
lattice::dotplot(cimod.full.ele.O)
R2.mod.full.ele.O<-r.squaredGLMM(mod.full.ele.O)
mod.full.ele.O
# 

#with lags
dat.MeanBodySize.ele$lags <- scale(dat.MeanBodySize.ele$lags)
mod.full.ele.lag <- lmer(lags ~ BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (1|article_id)+(BodyLengthMean|param/class),data=dat.MeanBodySize.ele,REML=TRUE,na.action="na.fail")
summary(mod.full.ele.lag)
rr <- ranef(mod.full.ele.lag, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.ele.lag <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.ele.lag <- confint(mod.full.ele.lag, method="Wald")
cimod.full.ele.lag <- cbind(fixef(mod.full.ele.lag),na.omit(cimod.full.ele.lag))
R2.mod.full.ele.lag<-r.squaredGLMM(mod.full.ele.lag)
R2.mod.full.ele.lag

# #latitude
par(mfrow=c(1,3))
hist(dat.MeanBodySize.lat$shift)
hist(dat.MeanBodySize.lat$v.lat.mean)
hist(dat.MeanBodySize.lat$BodyLengthMean)
# 
table(dat.MeanBodySize.lat$param,dat.MeanBodySize.lat$class)
mod.full.lat <- lmer(shift ~ v.lat.mean * BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (1|article_id)+(v.lat.mean * BodyLengthMean|param/class),data=dat.MeanBodySize.lat,REML=TRUE,na.action="na.fail")
summary(mod.full.lat)
rr <- ranef(mod.full.lat, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.lat <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.lat <- confint(mod.full.lat, method="Wald")
cimod.full.lat <- cbind(fixef(mod.full.lat),na.omit(cimod.full.lat))
lattice::dotplot(cimod.full.lat)
R2.mod.full.lat<-r.squaredGLMM(mod.full.lat)
R2.mod.full.lat

# mod.full.lat2 <- lmer(shift ~ v.lat.mean * BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter+uncertainty_distribution+(1|article_id)+(v.lat.mean * BodyLengthMean|param/class),data=dat.MeanBodySize.lat,REML=TRUE,na.action="na.fail")
# summary(mod.full.lat2)
# rr <- ranef(mod.full.lat2, condVar=TRUE)
# dd <- as.data.frame(rr)
# ddmod.full.lat2 <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
# lattice::dotplot(rr)
# cimod.full.lat2 <- confint(mod.full.lat2, method="Wald")
# cimod.full.lat2 <- cbind(fixef(mod.full.lat2),na.omit(cimod.full.lat2))

# lattice::dotplot(cimod.full.lat2)
# R2.mod.full.lat2<-r.squaredGLMM(mod.full.lat2)
# R2.mod.full.lat2

#with lags
dat.MeanBodySize.lat$lags <- scale(dat.MeanBodySize.lat$lags)
mod.full.lat.lag <- lmer(lags ~ BodyLengthMean +  dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (1|article_id)+(BodyLengthMean|param/class),data=dat.MeanBodySize.lat,REML=TRUE,na.action="na.fail")
summary(mod.full.lat.lag)
rr <- ranef(mod.full.lat.lag, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.lat.lag <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.lat.lag <- confint(mod.full.lat.lag, method="Wald")
cimod.full.lat.lag <- cbind(fixef(mod.full.lat.lag),na.omit(cimod.full.lat.lag))
R2.mod.full.lat.lag<-r.squaredGLMM(mod.full.lat.lag)
R2.mod.full.lat.lag

#
mod.full.lat.TE <- lmer(shift ~ v.lat.mean * BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (v.lat.mean * BodyLengthMean|class)+(1|article_id),data=lat_te,REML=TRUE,na.action="na.fail")
summary(mod.full.lat.TE)
rr <- ranef(mod.full.lat.TE, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.lat.TE <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.lat.TE <- confint(mod.full.lat.TE, method="Wald")
cimod.full.lat.TE <- cbind(fixef(mod.full.lat.TE),na.omit(cimod.full.lat.TE))
lattice::dotplot(cimod.full.lat.TE)
R2.mod.full.lat.TE<-r.squaredGLMM(mod.full.lat.TE)
mod.full.lat.TE

#library(jtools)
#interact_plot(mod.full.lat.TE,pred = v.lat.mean, modx = BodyLengthMean)
curve(rr$class[1,1]+rr$class[1,2]*x) #mean effect
yy <- range(lat_te$v.lat.mean)
curve(rr$class[1,1]+rr$class[1,2]*x+rr$class[1,3]*yy[1]+rr$class[1,3]*yy[1]*x,add=T,col="blue")
curve(rr$class[1,1]+rr$class[1,2]*x+rr$class[1,3]*yy[2]+rr$class[1,3]*yy[2]*x,add=T,col="red")

mod.full.lat.LE <- lmer(shift ~ v.lat.mean * BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (v.lat.mean * BodyLengthMean|class)+(1|article_id),data=lat_le,REML=TRUE,na.action="na.fail")
summary(mod.full.lat.LE)
rr <- ranef(mod.full.lat.LE, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.lat.LE <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.lat.LE <- confint(mod.full.lat.LE, method="Wald")
cimod.full.lat.LE <- cbind(fixef(mod.full.lat.LE),na.omit(cimod.full.lat.LE))
lattice::dotplot(cimod.full.lat.LE)
R2.mod.full.lat.LE<-r.squaredGLMM(mod.full.lat.LE)
R2.mod.full.lat.LE

mod.full.lat.O <- lmer(shift ~ v.lat.mean * BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (v.lat.mean * BodyLengthMean|class)+(1|article_id),data=lat_o,REML=TRUE,na.action="na.fail")
summary(mod.full.lat.O)
rr <- ranef(mod.full.lat.O, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.lat.O <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.lat.O <- confint(mod.full.lat.O, method="Wald")
cimod.full.lat.O <- cbind(fixef(mod.full.lat.O),na.omit(cimod.full.lat.O))
lattice::dotplot(cimod.full.lat.O)
R2.mod.full.lat.O<-r.squaredGLMM(mod.full.lat.O)
R2.mod.full.lat.O

#with lags
#
hist(lat_te$lags)
lat_te$lags <- scale(lat_te$lags)
mod.full.lat.TE.lag <- lmer(lags ~ BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (BodyLengthMean|class)+(1|article_id),data=lat_te,REML=TRUE,na.action="na.fail")
mod.full.lat.TE.lag <- lmer(lags ~ BodyLengthMean + (BodyLengthMean|class)+(1|article_id),data=lat_te,REML=TRUE,na.action="na.fail")
summary(mod.full.lat.TE.lag)
rr <- ranef(mod.full.lat.TE.lag, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.lat.TE.lag <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.lat.TE.lag <- confint(mod.full.lat.TE.lag, method="Wald")
cimod.full.lat.TE.lag <- cbind(fixef(mod.full.lat.TE.lag),na.omit(cimod.full.lat.TE.lag))
lattice::dotplot(cimod.full.lat.TE.lag)
R2.mod.full.lat.TE.lag<-r.squaredGLMM(mod.full.lat.TE.lag)


hist(lat_le$lags)
lat_le$lags <- scale(lat_le$lags)
mod.full.lat.LE.lag <- lmer(lags ~ BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (BodyLengthMean|class)+(1|article_id),data=lat_le,REML=TRUE,na.action="na.fail")
summary(mod.full.lat.LE.lag)
rr <- ranef(mod.full.lat.LE.lag, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.lat.LE.lag <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.lat.LE.lag <- confint(mod.full.lat.LE.lag, method="Wald")
cimod.full.lat.LE.lag <- cbind(fixef(mod.full.lat.LE.lag),na.omit(cimod.full.lat.LE.lag))
lattice::dotplot(cimod.full.lat.LE.lag)
R2.mod.full.lat.LE.lag<-r.squaredGLMM(mod.full.lat.LE.lag)

hist(lat_o$lags)
lat_o$lags <- scale(lat_o$lags)
mod.full.lat.O.lag <- lmer(lags ~ BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (BodyLengthMean|class)+(1|article_id),data=lat_o,REML=TRUE,na.action="na.fail")
summary(mod.full.lat.O.lag)
rr <- ranef(mod.full.lat.O.lag, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.lat.O.lag <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.lat.O.lag <- confint(mod.full.lat.O.lag, method="Wald")
cimod.full.lat.O.lag <- cbind(fixef(mod.full.lat.O.lag),na.omit(cimod.full.lat.O.lag))
lattice::dotplot(cimod.full.lat.O.lag)
R2.mod.full.lat.O.lag<-r.squaredGLMM(mod.full.lat.O.lag)
R2.mod.full.lat.O.lag

hist(ele_te$lags)
ele_te$lags <- scale(ele_te$lags)
mod.full.ele.TE.lag <- lmer(lags ~ BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (BodyLengthMean|class)+(1|article_id),data=ele_te,REML=TRUE,na.action="na.fail")
summary(mod.full.ele.TE.lag)
rr <- ranef(mod.full.ele.TE.lag, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.ele.TE.lag <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.ele.TE.lag <- confint(mod.full.ele.TE.lag, method="Wald")
cimod.full.ele.TE.lag <- cbind(fixef(mod.full.ele.TE.lag),na.omit(cimod.full.ele.TE.lag))
lattice::dotplot(cimod.full.ele.TE.lag)
R2.mod.full.ele.TE.lag<-r.squaredGLMM(mod.full.ele.TE.lag)

hist(ele_le$lags)
ele_le$lags <- scale(ele_le$lags)
mod.full.ele.LE.lag <- lmer(lags ~ BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (BodyLengthMean|class)+(1|article_id),data=ele_le,REML=TRUE,na.action="na.fail")
summary(mod.full.ele.LE.lag)
rr <- ranef(mod.full.ele.LE.lag, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.ele.LE.lag <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.ele.LE.lag <- confint(mod.full.ele.LE.lag, method="Wald")
cimod.full.ele.LE.lag <- cbind(fixef(mod.full.ele.LE.lag),na.omit(cimod.full.ele.LE.lag))
lattice::dotplot(cimod.full.ele.LE.lag)
R2.mod.full.ele.LE.lag<-r.squaredGLMM(mod.full.ele.LE.lag)

hist(ele_o$lags)
ele_o$lags <- scale(ele_o$lags)
mod.full.ele.O.lag <- lmer(lags ~ BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter + uncertainty_distribution + (BodyLengthMean|class)+(1|article_id),data=ele_o,REML=TRUE,na.action="na.fail")
summary(mod.full.ele.O.lag)
rr <- ranef(mod.full.ele.O.lag, condVar=TRUE)
dd <- as.data.frame(rr)
ddmod.full.ele.O.lag <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
lattice::dotplot(rr)
cimod.full.ele.O.lag <- confint(mod.full.ele.O.lag, method="Wald")
cimod.full.ele.O.lag <- cbind(fixef(mod.full.ele.O.lag),na.omit(cimod.full.ele.O.lag))
lattice::dotplot(cimod.full.ele.O.lag)
R2.mod.full.ele.O.lag<-r.squaredGLMM(mod.full.ele.O.lag)
R2.mod.full.ele.O.lag

#--------------------------------------------
#check convergence
#--------------------------------------------
library(AICcmodavg)
checkConv(mod.full.ele.LE)
checkConv(mod.full.ele.TE)
checkConv(mod.full.lat.LE)
checkConv(mod.full.lat.TE)
checkConv(mod.full.ele)
checkConv(mod.full.ele)
checkConv(mod.full.lat)
checkConv(mod.full.lat)
checkConv(mod.full.ele.lag)
checkConv(mod.full.ele.lag)
checkConv(mod.full.lat.lag)
checkConv(mod.full.lat.lag)

#---------------------------------------------
#concatenate datasets
#--------------------------------------------
R2 <- rbind(R2.mod.full.ele.LE,R2.mod.full.ele.TE,R2.mod.full.ele.O,R2.mod.full.lat.LE,R2.mod.full.lat.TE,R2.mod.full.lat.O)
R2 <- data.frame(R2,Gradient = c(rep("ELE",3),rep("LAT",3)),Type = rep(c("LE","TE","O"),2))


dd <- rbind(ddmod.full.ele.LE,ddmod.full.ele.TE,ddmod.full.ele.O,ddmod.full.lat.LE,ddmod.full.lat.TE,ddmod.full.lat.O)
dd <- data.frame(dd,Gradient = c(rep("ELE",nrow(ddmod.full.ele.LE)+nrow(ddmod.full.ele.TE)+nrow(ddmod.full.ele.O)),
                                 rep("LAT",nrow(ddmod.full.lat.LE)+nrow(ddmod.full.lat.TE)+nrow(ddmod.full.lat.O))),
                 Type = c(rep("LE",nrow(ddmod.full.ele.LE)),rep("TE",nrow(ddmod.full.ele.TE)),rep("O",nrow(ddmod.full.ele.O)),rep("LE",nrow(ddmod.full.lat.LE)),rep("TE",nrow(ddmod.full.lat.TE)),rep("O",nrow(ddmod.full.lat.O)))
                 
)
ci <- rbind(cimod.full.ele.LE,cimod.full.ele.TE,cimod.full.ele.O,cimod.full.lat.LE,cimod.full.lat.TE,cimod.full.lat.O)

ci <- data.frame(ci,Variable = rownames(ci),Gradient = c(rep("ELE",nrow(cimod.full.ele.LE)+nrow(cimod.full.ele.TE)+nrow(cimod.full.ele.O)),
                                                         rep("LAT",nrow(cimod.full.lat.LE)+nrow(cimod.full.lat.TE)+nrow(cimod.full.lat.O))),
                 Type = c(rep("LE",nrow(cimod.full.ele.LE)),rep("TE",nrow(cimod.full.ele.TE)),rep("O",nrow(cimod.full.ele.O)),rep("LE",nrow(cimod.full.lat.LE)),rep("TE",nrow(cimod.full.lat.TE)),rep("O",nrow(cimod.full.lat.O)))
)

ci <- na.omit(ci)

write.csv(ci,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps_full.csv",row.names=F)
write.csv(dd,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps_full.csv",row.names=F)
write.csv(R2,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps_full.csv",row.names=F)

#full models (one model)
R2 <- rbind(R2.mod.full.ele,R2.mod.full.lat)
R2 <- data.frame(R2,Gradient = c("ELE","LAT"))

dd <- rbind(ddmod.full.ele,ddmod.full.lat)
dd <- data.frame(dd,Gradient = c(rep("ELE",nrow(ddmod.full.ele)),
                                 rep("LAT",nrow(ddmod.full.lat))),
                 Type = sapply(strsplit(as.character(dd$grp),":"),'[',2),
                 Gp = sapply(strsplit(as.character(dd$grp),":"),'[',1))
                  
                 
dd$Type[is.na(dd$Type)] <- as.character(dd$grp[is.na(dd$Type)])

ci <- rbind(cimod.full.ele,cimod.full.lat)
ci <- data.frame(ci,Variable = rownames(ci),Gradient = c(rep("ELE",nrow(cimod.full.ele)),
                                                         rep("LAT",nrow(cimod.full.lat))))

write.csv(ci,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps_singlemodel.csv",row.names=F)
write.csv(dd,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps_singlemodel.csv",row.names=F)
write.csv(R2,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps_singlemodel.csv",row.names=F)

#lag model (one model)
R2 <- rbind(R2.mod.full.ele.lag,R2.mod.full.lat.lag)
R2 <- data.frame(R2,Gradient = c("ELE","LAT"))

dd <- rbind(ddmod.full.ele.lag,ddmod.full.lat.lag)
dd <- data.frame(dd,Gradient = c(rep("ELE",nrow(ddmod.full.ele.lag)),
                                 rep("LAT",nrow(ddmod.full.lat.lag))),
                 Type = sapply(strsplit(as.character(dd$grp),":"),'[',2),
                 Gp = sapply(strsplit(as.character(dd$grp),":"),'[',1))


dd$Type[is.na(dd$Type)] <- as.character(dd$grp[is.na(dd$Type)])

ci <- rbind(cimod.full.ele.lag,cimod.full.lat.lag)
ci <- data.frame(ci,Variable = rownames(ci),Gradient = c(rep("ELE",nrow(cimod.full.ele.lag)),
                                                         rep("LAT",nrow(cimod.full.lat.lag))))

write.csv(ci,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps_lag.csv",row.names=F)
write.csv(dd,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps_lag.csv",row.names=F)
write.csv(R2,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps_lag.csv",row.names=F)

#lag models
R2 <- rbind(R2.mod.full.ele.LE.lag,R2.mod.full.ele.TE.lag,R2.mod.full.ele.O.lag,R2.mod.full.lat.LE.lag,R2.mod.full.lat.TE.lag,R2.mod.full.lat.O.lag)
R2 <- data.frame(R2,Gradient = c(rep("ELE",3),rep("LAT",3)),Type = rep(c("LE","TE","O"),2))


dd <- rbind(ddmod.full.ele.LE.lag,ddmod.full.ele.TE.lag,ddmod.full.ele.O.lag,ddmod.full.lat.LE.lag,ddmod.full.lat.TE.lag,ddmod.full.lat.O.lag)
dd <- data.frame(dd,Gradient = c(rep("ELE",nrow(ddmod.full.ele.LE.lag)+nrow(ddmod.full.ele.TE.lag)+nrow(ddmod.full.ele.O.lag)),
                                 rep("LAT",nrow(ddmod.full.lat.LE.lag)+nrow(ddmod.full.lat.TE.lag)+nrow(ddmod.full.lat.O.lag))),
                 Type = c(rep("LE",nrow(ddmod.full.ele.LE.lag)),rep("TE",nrow(ddmod.full.ele.TE.lag)),rep("O",nrow(ddmod.full.ele.O.lag)),rep("LE",nrow(ddmod.full.lat.LE.lag)),rep("TE",nrow(ddmod.full.lat.TE.lag)),rep("O",nrow(ddmod.full.lat.O.lag)))
                 
)
ci <- rbind(cimod.full.ele.LE.lag,cimod.full.ele.TE.lag,cimod.full.ele.O.lag,cimod.full.lat.LE.lag,cimod.full.lat.TE.lag,cimod.full.lat.O.lag)

ci <- data.frame(ci,Variable = rownames(ci),Gradient = c(rep("ELE",nrow(cimod.full.ele.LE.lag)+nrow(cimod.full.ele.TE.lag)+nrow(cimod.full.ele.O.lag)),
                                                         rep("LAT",nrow(cimod.full.lat.LE.lag)+nrow(cimod.full.lat.TE.lag)+nrow(cimod.full.lat.O.lag))),
                 Type = c(rep("LE",nrow(cimod.full.ele.LE.lag)),rep("TE",nrow(cimod.full.ele.TE.lag)),rep("O",nrow(cimod.full.ele.O.lag)),rep("LE",nrow(cimod.full.lat.LE.lag)),rep("TE",nrow(cimod.full.lat.TE.lag)),rep("O",nrow(cimod.full.lat.O.lag)))
)

ci <- na.omit(ci)

write.csv(ci,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps_lags.csv",row.names=F)
write.csv(dd,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps_lags.csv",row.names=F)
write.csv(R2,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps_lags.csv",row.names=F)

#---------------------------------------------------------------------------------------------
library(plotrix); library(randomcoloR)

dd <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps_full.csv")
ci <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps_full.csv")
R2 <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps_full.csv")

n <- length(factor(c("AAA",sort(unique(as.character(dd$grp[dd$grpvar=="class"]))))))
palette <- data.frame(col=distinctColorPalette(n))
#palette$grp <- factor(c("Global",sort(unique(as.character(dd$grp[dd$grpvar=="class"])))))
dd$kingdom <- trait$kingdom[match(dd$grp,trait$class)]
palette$grp <- factor(c("Global",sort(unique(as.character(dd$grp[dd$grpvar=="class" & dd$kingdom=="Animalia"]))),sort(unique(as.character(dd$grp[dd$grpvar=="class" & dd$kingdom=="Plantae"])))))

pgp <- table(dd$grp[dd$grpvar=="class"],dd$Gradient[dd$grpvar=="class"])
palette$ELE <- ifelse(palette$grp %in% rownames(pgp)[which(pgp[,1]>0)],1,0)
palette$LAT <- ifelse(palette$grp %in% rownames(pgp)[which(pgp[,2]>0)],1,0)
palette$ELE[1] <- 1
palette$LAT[1] <- 1
palette$grp[1] <- "Global"

#----------------------------------
#latitude
#----------------------------------
#velocity
fixed <- ci[ci$Variable == 'v.lat.mean',c(1:3,6)]
fixed$class <- "Global"

names(fixed) <- c("condval","lwr","upr","Type","grp")
random <- dd[dd$term == 'v.lat.mean',c(4,6:7,9,3)]

datplot <- rbind(fixed[fixed$Type=="O",],random[random$Type=="O",],
                 fixed[fixed$Type=="TE",],random[random$Type=="TE",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",]
)


#                  
# x <- c(1:length(datplot$Type[datplot$Type=="O"]),
#        (length(datplot$Type[datplot$Type=="O"])+2):(length(datplot$Type[datplot$Type=="O"])+length(datplot$Type[datplot$Type=="LE"])+1),
#        ((length(datplot$Type[datplot$Type=="O"])+length(datplot$Type[datplot$Type=="LE"]))+2+1):((length(datplot$Type[datplot$Type=="O"])+2+length(datplot$Type[datplot$Type=="LE"]))+length(datplot$Type[datplot$Type=="TE"]))
# )
# 
# plotCI(x = x,y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=palette[match(datplot$grp,palette$grp),1],xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
# axis(2,tck=-0.02)
# box()
# abline(h=0,lty=2)
# abline(v=10,lty=1)
# abline(v=20,lty=1)
# legend("topleft",col=palette$col,legend=palette$grp,bty="n",pch=16,cex=1.5)
# 
#alternative figure
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelVelocity_allgps_full.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
#plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))
legend("topright",pch=c(16,17),legend=c("Leading","Trailing"),pt.cex=2,bty="n",cex=1.2)

mtext(outer=T,side=3,"Velocity",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()



#body size
fixed <- ci[ci$Variable == 'BodyLengthMean' & ci$Gradient=="LAT",c(1:3,6)]
fixed$class <- "Global"
names(fixed) <- c("condval","lwr","upr","Type","grp")
random <- dd[dd$term == 'BodyLengthMean'  & dd$Gradient=="LAT",c(4,6:7,9,3)]

datplot <- rbind(fixed[fixed$Type=="O" ,],random[random$Type=="O",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",],
                 fixed[fixed$Type=="TE" ,],random[random$Type=="TE",]
)


# 
# x <- c(1:length(datplot$Type[datplot$Type=="O"]),
#        (length(datplot$Type[datplot$Type=="O"])+2):(length(datplot$Type[datplot$Type=="O"])+length(datplot$Type[datplot$Type=="LE"])+1),
#        ((length(datplot$Type[datplot$Type=="O"])+length(datplot$Type[datplot$Type=="LE"]))+2+1):((length(datplot$Type[datplot$Type=="O"])+2+length(datplot$Type[datplot$Type=="LE"]))+length(datplot$Type[datplot$Type=="TE"]))
# )
# 
# plotCI(x = x,y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=palette[match(datplot$grp,palette$grp),1],xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
# axis(2,tck=-0.02)
# box()
# abline(h=0,lty=2)
# abline(v=10,lty=1)
# abline(v=20,lty=1)
# legend("bottomleft",col=palette$col,legend=palette$grp,bty="n",pch=16,cex=1.5)

#alternative figure
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelBodySize_allgps_full.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
#plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))
legend("topright",pch=c(16,17),legend=c("Leading","Trailing"),pt.cex=2,bty="n",cex=1.2)

mtext(outer=T,side=3,"BodySize",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()


#interaction body size:velocity
fixed <- ci[ci$Variable == 'v.lat.mean:BodyLengthMean' & ci$Gradient=="LAT",c(1:3,6)]
fixed$class <- "Global"
names(fixed) <- c("condval","lwr","upr","Type","grp")
random <- dd[dd$term == 'v.lat.mean:BodyLengthMean'  & dd$Gradient=="LAT",c(4,6:7,9,3)]

datplot <- rbind(fixed[fixed$Type=="O" ,],random[random$Type=="O",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",],
                 fixed[fixed$Type=="TE" ,],random[random$Type=="TE",]
)



# x <- c(1:length(datplot$Type[datplot$Type=="O"]),
#        (length(datplot$Type[datplot$Type=="O"])+2):(length(datplot$Type[datplot$Type=="O"])+length(datplot$Type[datplot$Type=="LE"])+1),
#        ((length(datplot$Type[datplot$Type=="O"])+length(datplot$Type[datplot$Type=="LE"]))+2+1):((length(datplot$Type[datplot$Type=="O"])+2+length(datplot$Type[datplot$Type=="LE"]))+length(datplot$Type[datplot$Type=="TE"]))
# )
# 
# plotCI(x = x,y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=palette[match(datplot$grp,palette$grp),1],xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
# axis(2,tck=-0.02)
# box()
# abline(h=0,lty=2)
# abline(v=10,lty=1)
# abline(v=20,lty=1)
# legend("bottomleft",col=palette$col,legend=palette$grp,bty="n",pch=16,cex=1.5)

#alternative figure
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelInteractionBodySizeVelocity_allgps_full.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
#plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))
legend("topright",pch=c(16,17),legend=c("Leading","Trailing"),pt.cex=2,bty="n",cex=1.2)

mtext(outer=T,side=3,"Velocity:BodySize",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()

#----------------------------------
#elevation
#----------------------------------
#velocity
fixed <- ci[ci$Variable == 'v.ele.mean',c(1:3,6)]
fixed$class <- "Global"

names(fixed) <- c("condval","lwr","upr","Type","grp")
random <- dd[dd$term == 'v.ele.mean',c(4,6:7,9,3)]

datplot <- rbind(fixed[fixed$Type=="O",],random[random$Type=="O",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",],
                 fixed[fixed$Type=="TE",],random[random$Type=="TE",]
)

#alternative figure
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelVelocity_allgps_full_ele.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$ELE==1],legend=palette$grp[palette$ELE==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))
legend("topright",pch=c(16,15,17),legend=c("Center","Leading","Trailing"),pt.cex=1.5,bty="n")

mtext(outer=T,side=3,"Velocity",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()



#body size
fixed <- ci[ci$Variable == 'BodyLengthMean' & ci$Gradient=="ELE",c(1:3,6)]
fixed$class <- "Global"
names(fixed) <- c("condval","lwr","upr","Type","grp")
random <- dd[dd$term == 'BodyLengthMean'  & dd$Gradient=="ELE",c(4,6:7,9,3)]

datplot <- rbind(fixed[fixed$Type=="O" ,],random[random$Type=="O",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",],
                 fixed[fixed$Type=="TE" ,],random[random$Type=="TE",]
)


#alternative figure
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelBodySize_allgps_full_ele.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$ELE==1],legend=palette$grp[palette$ELE==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))
mtext(outer=T,side=3,"BodySize",cex=1.3,line=-3,font=3,adj=0.1)
legend("topright",pch=c(16,15,17),legend=c("Center","Leading","Trailing"),pt.cex=1.5,bty="n")

dev.off()


#interaction body size:velocity
fixed <- ci[ci$Variable == 'v.ele.mean:BodyLengthMean' & ci$Gradient=="ELE",c(1:3,6)]
fixed$class <- "Global"
names(fixed) <- c("condval","lwr","upr","Type","grp")
random <- dd[dd$term == 'v.ele.mean:BodyLengthMean'  & dd$Gradient=="ELE",c(4,6:7,9,3)]

datplot <- rbind(fixed[fixed$Type=="O" ,],random[random$Type=="O",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",],
                 fixed[fixed$Type=="TE" ,],random[random$Type=="TE",]
)

#alternative figure
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelInteractionBodySizeVelocity_allgps_full_ele.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$ELE==1],legend=palette$grp[palette$ELE==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))
legend("topright",pch=c(16,15,17),legend=c("Center","Leading","Trailing"),pt.cex=1.5,bty="n")

mtext(outer=T,side=3,"Velocity:BodySize",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()

#------------------------------------------------------------------
#FULL MODELS
#-------------------------------------------------------------------
library(plotrix); library(randomcoloR)

dd <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps_singlemodel.csv")
ci <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps_singlemodel.csv")
R2 <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps_singlemodel.csv")

n <- length(factor(c("AAA",sort(unique(as.character(dd$Gp[dd$grpvar=="class:param"]))))))
palette <- data.frame(col=distinctColorPalette(n))
dd$kingdom <- trait$kingdom[match(dd$Gp,trait$class)]
palette$grp <- factor(c("Global",sort(unique(as.character(dd$Gp[dd$grpvar=="class:param" & dd$kingdom=="Animalia"]))),sort(unique(as.character(dd$Gp[dd$grpvar=="class:param" & dd$kingdom=="Plantae"])))))

pgp <- table(dd$Gp[dd$grpvar=="class:param"],dd$Gradient[dd$grpvar=="class:param"])
palette$ELE <- ifelse(palette$grp %in% rownames(pgp)[which(pgp[,1]>0)],1,0)
palette$LAT <- ifelse(palette$grp %in% rownames(pgp)[which(pgp[,2]>0)],1,0)
palette$ELE[1] <- 1
palette$LAT[1] <- 1
palette$grp[1] <- "Global"

#----------------------------------
#latitude
#----------------------------------
#velocity
# fixed <- ci[ci$Variable == 'v.lat.mean',c(1:3)]
# fixed$class <- "Global"
# fixed$grp <- 'all'
# names(fixed) <- c("condval","lwr","upr","Type",'Gp')


random <- dd[dd$term == 'v.lat.mean',c(4,6:7,9,10)]

fixed <- random[random$Gp %in% c("LE","O","TE"),]
random <- random[-which(random$Gp %in% c("LE","O","TE")),]

datplot <- rbind(fixed[fixed$Type=="O",],random[random$Type=="O",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",],
                 fixed[fixed$Type=="TE",],random[random$Type=="TE",]
)

datplot$grp <- datplot$Gp
datplot$grp[datplot$grp %in% c("LE","O","TE") ] <- "Global"
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelVelocity_allgps_singlemodel.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
#plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=7,inset=c(0,-0.15))
legend("topright",pch=c(16,17),legend=c("Leading","Trailing"),pt.cex=1.5,bty="n")

mtext(outer=T,side=3,"Velocity",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()

#body size
# fixed <- ci[ci$Variable == 'BodyLengthMean' & ci$Gradient=="LAT",c(1:3)]
# fixed$class <- "Global"
# fixed$grp <- 'all'
# names(fixed) <- c("condval","lwr","upr","Type",'Gp')

random <- dd[dd$term == 'BodyLengthMean' & dd$Gradient=="LAT",c(4,6:7,9,10)]

fixed <- random[random$Gp %in% c("LE","O","TE"),]
random <- random[-which(random$Gp %in% c("LE","O","TE")),]

datplot <- rbind(fixed[fixed$Type=="O",],random[random$Type=="O",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",],
                 fixed[fixed$Type=="TE",],random[random$Type=="TE",]
)

datplot$grp <- datplot$Gp
datplot$grp[datplot$grp %in% c("LE","O","TE") ] <- "Global"
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelBodySizeMean_allgps_singlemodel.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
#plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=7,inset=c(0,-0.15))
legend("topright",pch=c(16,17),legend=c("Leading","Trailing"),pt.cex=1.5,bty="n")

mtext(outer=T,side=3,"Body size",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()


#body size:velocity
# fixed <- ci[ci$Variable == 'BodyLengthMean' & ci$Gradient=="LAT",c(1:3)]
# fixed$class <- "Global"
# fixed$grp <- 'all'
# names(fixed) <- c("condval","lwr","upr","Type",'Gp')

random <- dd[dd$term == 'v.lat.mean:BodyLengthMean' & dd$Gradient=="LAT",c(4,6:7,9,10)]

fixed <- random[random$Gp %in% c("LE","O","TE"),]
random <- random[-which(random$Gp %in% c("LE","O","TE")),]

datplot <- rbind(fixed[fixed$Type=="O",],random[random$Type=="O",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",],
                 fixed[fixed$Type=="TE",],random[random$Type=="TE",]
)

datplot$grp <- datplot$Gp
datplot$grp[datplot$grp %in% c("LE","O","TE") ] <- "Global"
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelInteractionBodySizeMeanVelocity_allgps_singlemodel.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
#plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=7,inset=c(0,-0.15))
legend("topright",pch=c(16,17),legend=c("Leading","Trailing"),pt.cex=1.5,bty="n")

mtext(outer=T,side=3,"Velocity:Body size",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()

#--------------------------------------------------------------
#FULL LAGS
#---------------------------------------------------------------
library(plotrix); library(randomcoloR)

dd <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps_lag.csv")
ci <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps_lag.csv")
R2 <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps_lag.csv")

n <- length(factor(c("AAA",sort(unique(as.character(dd$Gp[dd$grpvar=="class:param"]))))))
palette <- data.frame(col=distinctColorPalette(n))
dd$kingdom <- trait$kingdom[match(dd$Gp,trait$class)]
palette$grp <- factor(c("Global",sort(unique(as.character(dd$Gp[dd$grpvar=="class:param" & dd$kingdom=="Animalia"]))),sort(unique(as.character(dd$Gp[dd$grpvar=="class:param" & dd$kingdom=="Plantae"])))))

pgp <- table(dd$Gp[dd$grpvar=="class:param"],dd$Gradient[dd$grpvar=="class:param"])
palette$ELE <- ifelse(palette$grp %in% rownames(pgp)[which(pgp[,1]>0)],1,0)
palette$LAT <- ifelse(palette$grp %in% rownames(pgp)[which(pgp[,2]>0)],1,0)
palette$ELE[1] <- 1
palette$LAT[1] <- 1
palette$grp[1] <- "Global"

#----------------------------------
#lags
#----------------------------------
#body size
# fixed <- ci[ci$Variable == dd$term == 'BodyLengthMean' & dd$Gradient=="LAT",c(1:3)]
# fixed$class <- "Global"
# fixed$grp <- 'all'
# names(fixed) <- c("condval","lwr","upr","Type",'Gp')


random <- dd[dd$term == 'BodyLengthMean' & dd$Gradient=="LAT",c(4,6:7,9,10)]

fixed <- random[random$Gp %in% c("LE","O","TE"),]
random <- random[-which(random$Gp %in% c("LE","O","TE")),]

datplot <- rbind(fixed[fixed$Type=="O",],random[random$Type=="O",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",],
                 fixed[fixed$Type=="TE",],random[random$Type=="TE",]
)

datplot$grp <- datplot$Gp
datplot$grp[datplot$grp %in% c("LE","O","TE") ] <- "Global"
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelBodySize_allgps_lag.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
#plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=7,inset=c(0,-0.15))
legend("topright",pch=c(16,17),legend=c("Leading","Trailing"),pt.cex=2,bty="n",cex=1.2)

mtext(outer=T,side=3,"Body Size",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()
# 

#------------------------------------------
#clean figures
#------------------------------------------
#---------------------------------------------------------------------------------------------
library(plotrix); library(randomcoloR)

dd <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps_full.csv")
ci <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps_full.csv")
R2 <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps_full.csv")

n <- length(factor(c("AAA",sort(unique(as.character(dd$grp[dd$grpvar=="class"]))))))
palette <- data.frame(col=distinctColorPalette(n))
dd$kingdom <- trait$kingdom[match(dd$grp,trait$class)]
palette$grp <- factor(c("Global",sort(unique(as.character(dd$grp[dd$grpvar=="class" & dd$kingdom=="Animalia"]))),sort(unique(as.character(dd$grp[dd$grpvar=="class" & dd$kingdom=="Plantae"])))))

pgp <- table(dd$grp[dd$grpvar=="class" & dd$Type %in% c("LE","TE")],dd$Gradient[dd$grpvar=="class"& dd$Type %in% c("LE","TE")])
palette$ELE <- ifelse(palette$grp %in% rownames(pgp)[which(pgp[,1]>0)],1,0)
palette$LAT <- ifelse(palette$grp %in% rownames(pgp)[which(pgp[,2]>0)],1,0)
palette$ELE[1] <- 1
palette$LAT[1] <- 1
palette$grp[1] <- "Global"

#----------------------------------
#latitude
#----------------------------------
#velocity
fixed <- ci[ci$Variable == 'v.lat.mean',c(1:3,6)]
fixed$class <- "Global"

names(fixed) <- c("condval","lwr","upr","Type","grp")
random <- dd[dd$term == 'v.lat.mean',c(4,6:7,9,3)]

datplot <- rbind(fixed[fixed$Type=="O",],random[random$Type=="O",],
                 fixed[fixed$Type=="TE",],random[random$Type=="TE",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",]
)

#clean table
aa <- table(datplot$grp,datplot$Type) 
rem <- names(which(apply(aa,1,sum)==1 & aa[,2]==1))
datplot <- datplot[-which(datplot$grp %in% rem),]
aa <- table(datplot$grp,datplot$Type) 
addTE <- names(which(aa[,3]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="TE",grp=addTE))
addLE <- names(which(aa[,1]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="LE",grp=addLE))
addO <- names(which(aa[,2]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="O",grp=addO))

#alternative figure
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelVelocity_allgps_full.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.5,axes=F)
#plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15),pt.cex=2)
legend("bottomleft",pch=c(17,15),legend=c("Trailing","Leading"),pt.cex=2.2,bty="n",cex=1.4)

mtext(outer=T,side=3,"Velocity",cex=1.5,line=-3,font=3,adj=0.1)
dev.off()



#body size
fixed <- ci[ci$Variable == 'BodyLengthMean' & ci$Gradient=="LAT",c(1:3,6)]
fixed$class <- "Global"
names(fixed) <- c("condval","lwr","upr","Type","grp")
random <- dd[dd$term == 'BodyLengthMean'  & dd$Gradient=="LAT",c(4,6:7,9,3)]

datplot <- rbind(fixed[fixed$Type=="O" ,],random[random$Type=="O",],
                 fixed[fixed$Type=="TE",],random[random$Type=="TE",],
                 fixed[fixed$Type=="LE" ,],random[random$Type=="LE",]
)

#clean table
aa <- table(datplot$grp,datplot$Type) 
rem <- names(which(apply(aa,1,sum)==1 & aa[,2]==1))
datplot <- datplot[-which(datplot$grp %in% rem),]
aa <- table(datplot$grp,datplot$Type) 
addTE <- names(which(aa[,3]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="TE",grp=addTE))
addLE <- names(which(aa[,1]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="LE",grp=addLE))
addO <- names(which(aa[,2]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="O",grp=addO))

#alternative figure
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelBodySize_allgps_full.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.5,axes=F)
#plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15),pt.cex=2)
legend("bottomleft",pch=c(17,15),legend=c("Trailing","Leading"),pt.cex=2.2,bty="n",cex=1.4)

mtext(outer=T,side=3,"BodySize",cex=1.5,line=-3,font=3,adj=0.1)
dev.off()


#interaction body size:velocity
fixed <- ci[ci$Variable == 'v.lat.mean:BodyLengthMean' & ci$Gradient=="LAT",c(1:3,6)]
fixed$class <- "Global"
names(fixed) <- c("condval","lwr","upr","Type","grp")
random <- dd[dd$term == 'v.lat.mean:BodyLengthMean'  & dd$Gradient=="LAT",c(4,6:7,9,3)]

datplot <- rbind(fixed[fixed$Type=="O" ,],random[random$Type=="O",],
                 fixed[fixed$Type=="TE",],random[random$Type=="TE",],
                 fixed[fixed$Type=="LE" ,],random[random$Type=="LE",]
)

#clean table
aa <- table(datplot$grp,datplot$Type) 
rem <- names(which(apply(aa,1,sum)==1 & aa[,2]==1))
datplot <- datplot[-which(datplot$grp %in% rem),]
aa <- table(datplot$grp,datplot$Type) 
addTE <- names(which(aa[,3]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="TE",grp=addTE))
addLE <- names(which(aa[,1]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="LE",grp=addLE))
addO <- names(which(aa[,2]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="O",grp=addO))

#alternative figure
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelInteractionBodySizeVelocity_allgps_full.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,
       pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),xlab="",
       ylab="Standardized coefficients",cex.lab=1.5,axes=F)
#plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15),pt.cex=2)
legend("bottomleft",pch=c(17,15),legend=c("Trailing","Leading"),pt.cex=2.2,bty="n",cex=1.4)

mtext(outer=T,side=3,"Velocity:BodySize",cex=1.5,line=-3,font=3,adj=0.1)
dev.off()

#----------------------------------
#with the lags
library(plotrix); library(randomcoloR)

dd <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps_lags.csv")
ci <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps_lags.csv")
R2 <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps_lags.csv")

n <- length(factor(c("AAA",sort(unique(as.character(dd$grp[dd$grpvar=="class"]))))))
palette <- data.frame(col=distinctColorPalette(n))
dd$kingdom <- trait$kingdom[match(dd$grp,trait$class)]
palette$grp <- factor(c("Global",sort(unique(as.character(dd$grp[dd$grpvar=="class" & dd$kingdom=="Animalia"]))),sort(unique(as.character(dd$grp[dd$grpvar=="class" & dd$kingdom=="Plantae"])))))

pgp <- table(dd$grp[dd$grpvar=="class" & dd$Type %in% c("LE","TE")],dd$Gradient[dd$grpvar=="class"& dd$Type %in% c("LE","TE")])
palette$ELE <- ifelse(palette$grp %in% rownames(pgp)[which(pgp[,1]>0)],1,0)
palette$LAT <- ifelse(palette$grp %in% rownames(pgp)[which(pgp[,2]>0)],1,0)
palette$ELE[1] <- 1
palette$LAT[1] <- 1
palette$grp[1] <- "Global"

#----------------------------------
#latitude
#----------------------------------
#body size
fixed <- ci[ci$Variable == 'BodyLengthMean' & ci$Gradient=="LAT",c(1:3,6)]
fixed$class <- "Global"
names(fixed) <- c("condval","lwr","upr","Type","grp")
random <- dd[dd$term == 'BodyLengthMean'  & dd$Gradient=="LAT",c(4,6:7,9,3)]

datplot <- rbind(fixed[fixed$Type=="O" ,],random[random$Type=="O",],
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",],
                 fixed[fixed$Type=="TE" ,],random[random$Type=="TE",]
)

#clean table
aa <- table(datplot$grp,datplot$Type) 
rem <- names(which(apply(aa,1,sum)==1 & aa[,2]==1))
datplot <- datplot[-which(datplot$grp %in% rem),]
aa <- table(datplot$grp,datplot$Type) 
addTE <- names(which(aa[,3]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="TE",grp=addTE))
addLE <- names(which(aa[,1]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="LE",grp=addLE))
addO <- names(which(aa[,2]==0))
datplot <- rbind(datplot,data.frame(condval=NA,lwr=NA,upr=NA,Type="O",grp=addO))

#alternative figure
datplot$grp <- factor(datplot$grp,levels=palette$grp)
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelBodySize_allgps_lags.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.5,axes=F)
#plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15),pt.cex=2)
legend("topright",pch=c(16,17),legend=c("Leading","Trailing"),pt.cex=2,bty="n",cex=1.2)

mtext(outer=T,side=3,"BodySize",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()


