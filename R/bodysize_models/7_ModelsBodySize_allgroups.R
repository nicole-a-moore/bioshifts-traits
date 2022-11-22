bod<-read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/CompilationBodySizeBioshiftsv1harmonized.csv")
trait <- read.csv('G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Bioshiftv1_checked.csv')
trait$spnameoriginal <- trait$new_name
trait$SpeciesChecked[is.na(trait$SpeciesChecked)] <- trait$new_name[is.na(trait$SpeciesChecked)] #a few ID problem

#----------------------
#merge different classes of body sizes
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

#----------------------
#transform body size
BodyLengthMean <- bod[bod$Code == "BodyLengthMean",]
BodyLengthMean$BodyLength <- as.numeric(as.character(BodyLengthMean$BodyLength))
BodyLengthMean <- BodyLengthMean[which(is.na(BodyLengthMean$BodyLength)==F),]
hist(BodyLengthMean$BodyLength)
BodyLengthMean$BodyLength <- log(BodyLengthMean$BodyLength)
hist(BodyLengthMean$BodyLength)

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


#prepare datasets
#----------------------------------
##select only classes with >5 species
dat.MeanBodySize.ele <- data.frame(ele[which(is.na(ele$BodyLengthMean)==F),])
dat.MeanBodySize.ele <- dat.MeanBodySize.ele[which(is.na(dat.MeanBodySize.ele$v.ele.mean)==F),]

dsp <- unique(dat.MeanBodySize.ele[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >10))
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
##select only classes with >5 species
dat.MeanBodySize.lat <- data.frame(lat[which(is.na(lat$BodyLengthMean)==F),])

dsp <- unique(dat.MeanBodySize.lat[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >10))
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

#-----------------------------------
lat_le <- lat_le[which(is.na(lat_le$BodyLengthMean)==F),]
dsp <- unique(lat_le[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >10))
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
nn <- names(which(cl >10))
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
nn <- names(which(cl >10))
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
nn <- names(which(cl >10))
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
nn <- names(which(cl >10))
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
nn <- names(which(cl >10))
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


#--------------------------
#models
#--------------------------
library(lme4);library(lmerTest); library(corrplot); library(MuMIn)

#simple models with lags
#article_id sucks all the variance
m <- lmer(lags ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=dat.MeanBodySize.ele)
summary(m)

m <- lmer(lags ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=ele_le)
summary(m)

m <- lmer(lags ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=ele_te)
summary(m)

m <- lmer(lags ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=dat.MeanBodySize.lat)
summary(m)

m <- lmer(lags ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=lat_le)
summary(m)

m <- lmer(lags ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=lat_te)
summary(m)


#----------------------------------------

#simple models with shifts
m <- lmer(shift ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=dat.MeanBodySize.ele)
summary(m)

m <- lmer(shift ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=ele_le)
summary(m)

m <- lmer(shift ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=ele_te)
summary(m)

m <- lmer(shift ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=dat.MeanBodySize.lat)
summary(m)

m <- lmer(shift ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=lat_le)
summary(m)

m <- lmer(shift ~ BodyLengthMean + (1|article_id)+(BodyLengthMean|class),data=lat_te)
summary(m)


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
table(dat.MeanBodySize.lat$article_id)#not enough repetitions

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
# mod.full.ele <- lmer(shift ~ v.ele.mean * BodyLengthMean + (1|article_id)+(v.ele.mean * BodyLengthMean|param/class),data=dat.MeanBodySize.ele,REML=TRUE,na.action="na.fail")
# summary(mod.full.ele)
# rr <- ranef(mod.full.ele, condVar=TRUE)
# dd <- as.data.frame(rr)
# ddmod.full.ele <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
# lattice::dotplot(rr)
# R2.mod.full.ele<-r.squaredGLMM(mod.full.ele)
# R2.mod.full.ele
# 
# mod.full.ele2 <- lmer(shift ~ v.ele.mean * BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter+uncertainty_distribution+(1|article_id)+(v.ele.mean * BodyLengthMean|param/class),data=dat.MeanBodySize.ele,REML=TRUE,na.action="na.fail")
# summary(mod.full.ele2)
# rr <- ranef(mod.full.ele2, condVar=TRUE)
# dd <- as.data.frame(rr)
# ddmod.full.ele2 <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
# lattice::dotplot(rr)
# R2.mod.full.ele2<-r.squaredGLMM(mod.full.ele2)
# R2.mod.full.ele2


mod.full.ele.TE <- lmer(shift ~ v.ele.mean * BodyLengthMean + (v.ele.mean * BodyLengthMean|class)+(1|article_id),data=ele_te,REML=TRUE,na.action="na.fail")
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

mod.full.ele.O <- lmer(shift ~ v.ele.mean * BodyLengthMean + (v.ele.mean * BodyLengthMean|class)+(1|article_id),data=ele_o,REML=TRUE,na.action="na.fail")
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


# #latitude
par(mfrow=c(1,3))
hist(dat.MeanBodySize.lat$shift)
hist(dat.MeanBodySize.lat$v.lat.mean)
hist(dat.MeanBodySize.lat$BodyLengthMean)
# 
# table(dat.MeanBodySize.lat$param,dat.MeanBodySize.lat$class)
# mod.full.lat <- lmer(shift ~ v.lat.mean * BodyLengthMean + (1|article_id)+(v.lat.mean * BodyLengthMean|param/class),data=dat.MeanBodySize.lat,REML=TRUE,na.action="na.fail")
# summary(mod.full.lat)
# rr <- ranef(mod.full.lat, condVar=TRUE)
# dd <- as.data.frame(rr)
# ddmod.full.lat <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
# lattice::dotplot(rr)
# cimod.full.lat <- confint(mod.full.lat)
# cimod.full.lat <- cbind(fixef(mod.full.lat),na.omit(cimod.full.lat))

# lattice::dotplot(cimod.full.lat)
# R2.mod.full.lat<-r.squaredGLMM(mod.full.lat)
# R2.mod.full.lat
# 
# mod.full.lat2 <- lmer(shift ~ v.lat.mean * BodyLengthMean + dur + id.area + data + sampling + uncertainty_parameter+uncertainty_distribution+(1|article_id)+(v.lat.mean * BodyLengthMean|param/class),data=dat.MeanBodySize.lat,REML=TRUE,na.action="na.fail")
# summary(mod.full.lat2)
# rr <- ranef(mod.full.lat2, condVar=TRUE)
# dd <- as.data.frame(rr)
# ddmod.full.lat2 <- transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
# lattice::dotplot(rr)
# cimod.full.lat2 <- confint(mod.full.lat2)
# cimod.full.lat2 <- cbind(fixef(mod.full.lat2),na.omit(cimod.full.lat2))

# lattice::dotplot(cimod.full.lat2)
# R2.mod.full.lat2<-r.squaredGLMM(mod.full.lat2)
# R2.mod.full.lat2


mod.full.lat.TE <- lmer(shift ~ v.lat.mean * BodyLengthMean + (v.lat.mean * BodyLengthMean|class)+(1|article_id),data=lat_te,REML=TRUE,na.action="na.fail")
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

mod.full.lat.LE <- lmer(shift ~ v.lat.mean * BodyLengthMean + (v.lat.mean * BodyLengthMean|class)+(1|article_id),data=lat_le,REML=TRUE,na.action="na.fail")
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

mod.full.lat.O <- lmer(shift ~ v.lat.mean * BodyLengthMean + (v.lat.mean * BodyLengthMean|class)+(1|article_id),data=lat_o,REML=TRUE,na.action="na.fail")
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

write.csv(ci,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps.csv",row.names=F)
write.csv(dd,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps.csv",row.names=F)
write.csv(R2,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps.csv",row.names=F)

#---------------------------------------------------------------------------------------------
library(plotrix); library(randomcoloR)

dd <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps.csv")
ci <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps.csv")
R2 <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps.csv")

n <- length(factor(c("AAA",sort(unique(as.character(dd$grp[dd$grpvar=="class"]))))))
palette <- data.frame(col=distinctColorPalette(n))
palette$grp <- factor(c("Global",sort(unique(as.character(dd$grp[dd$grpvar=="class"])))))

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
                 fixed[fixed$Type=="LE",],random[random$Type=="LE",],
                 fixed[fixed$Type=="TE",],random[random$Type=="TE",]
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


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelVelocity_allgps.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))
legend("topright",pch=c(16,15,17),legend=c("Center","Leading","Trailing"),pt.cex=1.5,bty="n")

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

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelBodySize_allgps.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))
legend("topright",pch=c(16,15,17),legend=c("Center","Leading","Trailing"),pt.cex=1.5,bty="n")

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

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelInteractionBodySizeVelocity_allgps.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02,las=2)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col[palette$LAT==1],legend=palette$grp[palette$LAT==1],bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))
legend("topright",pch=c(16,15,17),legend=c("Center","Leading","Trailing"),pt.cex=1.5,bty="n")

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


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelVelocity_allgps_ele.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

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

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelBodySize_allgps_ele.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

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

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelInteractionBodySizeVelocity_allgps_ele.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

box()
abline(h=0,lty=2)

par(xpd=NA)
legend("topright",pch=c(16,15,17),legend=c("Center","Leading","Trailing"),pt.cex=1.5,bty="n")

mtext(outer=T,side=3,"Velocity:BodySize",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()
