bod<-read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/CompilationBodySizeBioshiftsv1harmonized.csv")
trait <- read.csv('G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Bioshiftv1_checked.csv')
trait$spnameoriginal <- trait$new_name
trait$SpeciesChecked[is.na(trait$SpeciesChecked)] <- trait$new_name[is.na(trait$SpeciesChecked)] #a few ID problem

#----------------------
#correct a few body sizes
bod$Code[bod$class=="Phaeophyceae"] <- 'VegetativeHeight'

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
BodyLengthMax <- bod[bod$Code == "BodyLengthMax",]
VegetativeHeight <- bod[bod$Code == "VegetativeHeight",]
Max_Height <- bod[bod$Code == "Max_Height",]
WingspanMean <- bod[bod$Code == "WingspanMean",]
WingspanMax <- bod[bod$Code == "WingspanMax",]

BodyLengthMean$BodyLength <- as.numeric(as.character(BodyLengthMean$BodyLength))
BodyLengthMax$BodyLength <- as.numeric(as.character(BodyLengthMax$BodyLength))
VegetativeHeight$BodyLength <- as.numeric(as.character(VegetativeHeight$BodyLength))
Max_Height$BodyLength <- as.numeric(as.character(Max_Height$BodyLength))
WingspanMean$BodyLength <- as.numeric(as.character(WingspanMean$BodyLength))
WingspanMax$BodyLength <- as.numeric(as.character(WingspanMax$BodyLength))

BodyLengthMean <- BodyLengthMean[which(is.na(BodyLengthMean$BodyLength)==F),]
BodyLengthMax <- BodyLengthMax[which(is.na(BodyLengthMax$BodyLength)==F),]
VegetativeHeight <- VegetativeHeight[which(is.na(VegetativeHeight$BodyLength)==F),]
Max_Height <- Max_Height[which(is.na(Max_Height$BodyLength)==F),]
WingspanMean <- WingspanMean[which(is.na(WingspanMean$BodyLength)==F),]
WingspanMax <- WingspanMax[which(is.na(WingspanMax$BodyLength)==F),]

#par(mfrow=c(3,2))
hist(BodyLengthMean$BodyLength)
hist(BodyLengthMax$BodyLength)
hist(VegetativeHeight$BodyLength)
hist(Max_Height$BodyLength)
hist(WingspanMean$BodyLength)
hist(WingspanMax$BodyLength)

BodyLengthMean$BodyLength <- log(BodyLengthMean$BodyLength)
BodyLengthMax$BodyLength <- log(BodyLengthMax$BodyLength)
VegetativeHeight$BodyLength <- log(VegetativeHeight$BodyLength)
Max_Height$BodyLength <- log(Max_Height$BodyLength)
WingspanMean$BodyLength <- log(WingspanMean$BodyLength)
WingspanMax$BodyLength <- log(WingspanMax$BodyLength)

#par(mfrow=c(3,2))
hist(BodyLengthMean$BodyLength)
hist(BodyLengthMax$BodyLength)
hist(VegetativeHeight$BodyLength)
hist(Max_Height$BodyLength)
hist(WingspanMean$BodyLength)
hist(WingspanMax$BodyLength)

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

# ele$v.ele.mean <- scale(ele$v.ele.mean)
# ele$shift <- scale(ele$shift)
# ele$lags <- scale(ele$lags)
# 
# lat$v.lat.mean <- scale(lat$v.lat.mean)
# lat$shift <- scale(lat$shift)
# lat$lags <- scale(lat$lags)


#add body size
ele$BodyLengthMean <- BodyLengthMean$BodyLength[match(ele$SpeciesChecked,BodyLengthMean$SpeciesChecked)]
ele$BodyLengthMax <- BodyLengthMax$BodyLength[match(ele$SpeciesChecked,BodyLengthMax$SpeciesChecked)]
ele$VegetativeHeight <- VegetativeHeight$BodyLength[match(ele$SpeciesChecked,VegetativeHeight$SpeciesChecked)]
ele$Max_Height <- Max_Height$BodyLength[match(ele$SpeciesChecked,Max_Height$SpeciesChecked)]
ele$WingspanMean <- WingspanMean$BodyLength[match(ele$SpeciesChecked,WingspanMean$SpeciesChecked)]
ele$WingspanMax <- WingspanMax$BodyLength[match(ele$SpeciesChecked,WingspanMax$SpeciesChecked)]

lat$BodyLengthMean <- BodyLengthMean$BodyLength[match(lat$SpeciesChecked,BodyLengthMean$SpeciesChecked)]
lat$BodyLengthMax <- BodyLengthMax$BodyLength[match(lat$SpeciesChecked,BodyLengthMax$SpeciesChecked)]
lat$VegetativeHeight <- VegetativeHeight$BodyLength[match(lat$SpeciesChecked,VegetativeHeight$SpeciesChecked)]
lat$Max_Height <- Max_Height$BodyLength[match(lat$SpeciesChecked,Max_Height$SpeciesChecked)]
lat$WingspanMean <- WingspanMean$BodyLength[match(lat$SpeciesChecked,WingspanMean$SpeciesChecked)]
lat$WingspanMax <- WingspanMax$BodyLength[match(lat$SpeciesChecked,WingspanMax$SpeciesChecked)]

# #scale variables
# ele$BodyLengthMean <- scale(ele$BodyLengthMean)
# ele$BodyLengthMax <- scale(ele$BodyLengthMax)
# ele$VegetativeHeight <- scale(ele$VegetativeHeight)
# ele$Max_Height <- scale(ele$Max_Height)
# ele$WingspanMean <- scale(ele$WingspanMean)
# ele$WingspanMax <- scale(ele$WingspanMax)
# ele$n <- scale(ele$n)
# ele$id.area <- scale(ele$id.area)
# ele$start <- scale(ele$start)
# ele$dur <- scale(ele$dur)
# ele$shift <- scale(ele$shift)
# 
# 
# lat$BodyLengthMean <- scale(lat$BodyLengthMean)
# lat$BodyLengthMax <- scale(lat$BodyLengthMax)
# lat$VegetativeHeight <- scale(lat$VegetativeHeight)
# lat$Max_Height <- scale(lat$Max_Height)
# lat$WingspanMean <- scale(lat$WingspanMean)
# lat$WingspanMax <- scale(lat$WingspanMax)
# lat$n <- scale(lat$n)
# lat$id.area <- scale(lat$id.area)
# lat$start <- scale(lat$start)
# lat$dur <- scale(lat$dur)
# 
# lat$shift <- scale(lat$shift)
# 

#----------------------
#model with body size
library(lme4); library(lmerTest); library(corrplot); library(MuMIn)

ele_te <- ele[ele$param=="TE",]
ele_le <- ele[ele$param=="LE",]
ele_o <- ele[ele$param=="O",]
lat_te <- lat[lat$param=="TE",]
lat_le <- lat[lat$param=="LE",]
lat_o <- lat[lat$param=="O",]

#animals
par(mfrow=c(2,2))
plot(ele_le$lags ~ ele_le$BodyLengthMean,col=factor(ele_le$class))
plot(ele_te$lags ~ ele_te$BodyLengthMean,col=factor(ele_te$class))

plot(lat_le$lags ~ lat_le$BodyLengthMean,col=factor(lat_le$class))
plot(lat_te$lags ~ lat_te$BodyLengthMean,col=factor(lat_te$class))


par(mfrow=c(2,2))
plot(ele_le$shift ~ ele_le$BodyLengthMean,col=factor(ele_le$class))
plot(ele_te$shift ~ ele_te$BodyLengthMean,col=factor(ele_te$class))

plot(lat_le$shift ~ lat_le$BodyLengthMean,col=factor(lat_le$class))
plot(lat_te$shift ~ lat_te$BodyLengthMean,col=factor(lat_te$class))

#plants

par(mfrow=c(2,2))
plot(ele_le$lags ~ ele_le$VegetativeHeight,col=factor(ele_le$class))
plot(ele_te$lags ~ ele_te$VegetativeHeight,col=factor(ele_te$class))

plot(lat_le$lags ~ lat_le$VegetativeHeight,col=factor(lat_le$class))
plot(lat_te$lags ~ lat_te$VegetativeHeight,col=factor(lat_te$class))


par(mfrow=c(2,2))
plot(ele_le$shift ~ ele_le$VegetativeHeight,col=factor(ele_le$class))
plot(ele_te$shift ~ ele_te$VegetativeHeight,col=factor(ele_te$class))

plot(lat_le$shift ~ lat_le$VegetativeHeight,col=factor(lat_le$class))
plot(lat_te$shift ~ lat_te$VegetativeHeight,col=factor(lat_te$class))

#prepare datasets
#----------------------------------
##select only classes with >5 species
dat.MeanBodySize.ele <- data.frame(ele[which(is.na(ele$BodyLengthMean)==F),])
dat.MeanBodySize.ele <- dat.MeanBodySize.ele[which(is.na(dat.MeanBodySize.ele$v.ele.mean)==F),]

dsp <- unique(dat.MeanBodySize.ele[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
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
nn <- names(which(cl >5))
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

# #----------------------------------
# ##select only classes with >5 species
# dat.MaxBodySize.ele <- data.frame(ele[which(is.na(ele$BodyLengthMax)==F),])
# dat.MaxBodySize.ele <- dat.MaxBodySize.ele[which(is.na(dat.MaxBodySize.ele$v.ele.mean)==F),]
# 
# dsp <- unique(dat.MaxBodySize.ele[,c(57,66)])
# cl <- table(dsp$class)
# nn <- names(which(cl >5))
# dat.MaxBodySize.ele <- dat.MaxBodySize.ele[which(dat.MaxBodySize.ele$class %in% nn),]
# dat.MaxBodySize.ele <-droplevels(dat.MaxBodySize.ele)
# 
# #scale variables
# dat.MaxBodySize.ele$BodyLengthMax <- scale(dat.MaxBodySize.ele$BodyLengthMax)
# dat.MaxBodySize.ele$n <- scale(dat.MaxBodySize.ele$n)
# dat.MaxBodySize.ele$id.area <- scale(dat.MaxBodySize.ele$id.area)
# dat.MaxBodySize.ele$start <- scale(dat.MaxBodySize.ele$start)
# dat.MaxBodySize.ele$dur <- scale(dat.MaxBodySize.ele$dur)
# dat.MaxBodySize.ele$shift <- scale(dat.MaxBodySize.ele$shift)
# dat.MaxBodySize.ele$v.ele.mean <- scale(dat.MaxBodySize.ele$v.ele.mean)

# #----------------------------------
# ##select only classes with >5 species
# dat.MaxBodySize.lat <- data.frame(lat[which(is.na(lat$BodyLengthMax)==F),])
# dsp <- unique(dat.MaxBodySize.lat[,c(57,66)])
# cl <- table(dsp$class)
# nn <- names(which(cl >5))
# dat.MaxBodySize.lat <- dat.MaxBodySize.lat[which(dat.MaxBodySize.lat$class %in% nn),]
# dat.MaxBodySize.lat <-droplevels(dat.MaxBodySize.lat)
# 
# #scale variables
# dat.MaxBodySize.lat$BodyLengthMax <- scale(dat.MaxBodySize.lat$BodyLengthMax)
# dat.MaxBodySize.lat$n <- scale(dat.MaxBodySize.lat$n)
# dat.MaxBodySize.lat$id.area <- scale(dat.MaxBodySize.lat$id.area)
# dat.MaxBodySize.lat$start <- scale(dat.MaxBodySize.lat$start)
# dat.MaxBodySize.lat$dur <- scale(dat.MaxBodySize.lat$dur)
# dat.MaxBodySize.lat$shift <- scale(dat.MaxBodySize.lat$shift)
# dat.MaxBodySize.lat$v.lat.mean <- scale(dat.MaxBodySize.lat$v.lat.mean)

#----------------------------------
##select only classes with >5 species
dat.VegetativeHeight.ele <- data.frame(ele[which(is.na(ele$VegetativeHeight)==F),])
dat.VegetativeHeight.ele <- dat.VegetativeHeight.ele[which(is.na(dat.VegetativeHeight.ele$v.ele.mean)==F),]

dsp <- unique(dat.VegetativeHeight.ele[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
dat.VegetativeHeight.ele <- dat.VegetativeHeight.ele[which(dat.VegetativeHeight.ele$class %in% nn),]
dat.VegetativeHeight.ele <-droplevels(dat.VegetativeHeight.ele)

#scale variables
dat.VegetativeHeight.ele$VegetativeHeight <- scale(dat.VegetativeHeight.ele$VegetativeHeight)
dat.VegetativeHeight.ele$n <- scale(dat.VegetativeHeight.ele$n)
dat.VegetativeHeight.ele$id.area <- scale(dat.VegetativeHeight.ele$id.area)
dat.VegetativeHeight.ele$start <- scale(dat.VegetativeHeight.ele$start)
dat.VegetativeHeight.ele$dur <- scale(dat.VegetativeHeight.ele$dur)
dat.VegetativeHeight.ele$shift <- scale(dat.VegetativeHeight.ele$shift)
dat.VegetativeHeight.ele$v.ele.mean <- scale(dat.VegetativeHeight.ele$v.ele.mean)

#----------------------------------
##select only classes with >5 species
dat.VegetativeHeight.lat <- data.frame(lat[which(is.na(lat$VegetativeHeight)==F),])
dsp <- unique(dat.VegetativeHeight.lat[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
dat.VegetativeHeight.lat <- dat.VegetativeHeight.lat[which(dat.VegetativeHeight.lat$class %in% nn),]
dat.VegetativeHeight.lat <-droplevels(dat.VegetativeHeight.lat)

#scale variables
dat.VegetativeHeight.lat$VegetativeHeight <- scale(dat.VegetativeHeight.lat$VegetativeHeight)
dat.VegetativeHeight.lat$n <- scale(dat.VegetativeHeight.lat$n)
dat.VegetativeHeight.lat$id.area <- scale(dat.VegetativeHeight.lat$id.area)
dat.VegetativeHeight.lat$start <- scale(dat.VegetativeHeight.lat$start)
dat.VegetativeHeight.lat$dur <- scale(dat.VegetativeHeight.lat$dur)
dat.VegetativeHeight.lat$shift <- scale(dat.VegetativeHeight.lat$shift)
dat.VegetativeHeight.lat$v.lat.mean <- scale(dat.VegetativeHeight.lat$v.lat.mean)

# #----------------------------------
# ##select only classes with >5 species
# dat.Max_Height.ele <- data.frame(ele[which(is.na(ele$Max_Height)==F),])
# dat.Max_Height.ele <- dat.Max_Height.ele[which(is.na(dat.Max_Height.ele$v.ele.mean)==F),]
# dsp <- unique(dat.Max_Height.ele[,c(57,66)])
# cl <- table(dsp$class)
# nn <- names(which(cl >5))
# dat.Max_Height.ele <- dat.Max_Height.ele[which(dat.Max_Height.ele$class %in% nn),]
# dat.Max_Height.ele <-droplevels(dat.Max_Height.ele)
# 
# #scale variables
# dat.Max_Height.ele$Max_Height <- scale(dat.Max_Height.ele$Max_Height)
# dat.Max_Height.ele$n <- scale(dat.Max_Height.ele$n)
# dat.Max_Height.ele$id.area <- scale(dat.Max_Height.ele$id.area)
# dat.Max_Height.ele$start <- scale(dat.Max_Height.ele$start)
# dat.Max_Height.ele$dur <- scale(dat.Max_Height.ele$dur)
# dat.Max_Height.ele$shift <- scale(dat.Max_Height.ele$shift)
# dat.Max_Height.ele$v.ele.mean <- scale(dat.Max_Height.ele$v.ele.mean)

# #----------------------------------
# ##select only classes with >5 species
# dat.Max_Height.lat <- data.frame(lat[which(is.na(lat$Max_Height)==F),])
# dsp <- unique(dat.Max_Height.lat[,c(57,66)])
# cl <- table(dsp$class)
# nn <- names(which(cl >5))
# dat.Max_Height.lat <- dat.Max_Height.lat[which(dat.Max_Height.lat$class %in% nn),]
# dat.Max_Height.lat <-droplevels(dat.Max_Height.lat)
# 
# #scale variables
# dat.Max_Height.lat$Max_Height <- scale(dat.Max_Height.lat$Max_Height)
# dat.Max_Height.lat$n <- scale(dat.Max_Height.lat$n)
# dat.Max_Height.lat$id.area <- scale(dat.Max_Height.lat$id.area)
# dat.Max_Height.lat$start <- scale(dat.Max_Height.lat$start)
# dat.Max_Height.lat$dur <- scale(dat.Max_Height.lat$dur)
# dat.Max_Height.lat$shift <- scale(dat.Max_Height.lat$shift)
# dat.Max_Height.lat$v.lat.mean <- scale(dat.Max_Height.lat$v.lat.mean)

# #----------------------------------
# ##select only classes with >5 species
# dat.WingspanMean.ele <- data.frame(ele[which(is.na(ele$WingspanMean)==F),])
# dat.WingspanMean.ele <- dat.WingspanMean.ele[which(is.na(dat.WingspanMean.ele$v.ele.mean)==F),]
# dsp <- unique(dat.WingspanMean.ele[,c(57,66)])
# cl <- table(dsp$class)
# nn <- names(which(cl >5))
# dat.WingspanMean.ele <- dat.WingspanMean.ele[which(dat.WingspanMean.ele$class %in% nn),]
# dat.WingspanMean.ele <-droplevels(dat.WingspanMean.ele)
# 
# 
# #scale variables
# dat.WingspanMean.ele$WingspanMean <- scale(dat.WingspanMean.ele$WingspanMean)
# dat.WingspanMean.ele$n <- scale(dat.WingspanMean.ele$n)
# dat.WingspanMean.ele$id.area <- scale(dat.WingspanMean.ele$id.area)
# dat.WingspanMean.ele$start <- scale(dat.WingspanMean.ele$start)
# dat.WingspanMean.ele$dur <- scale(dat.WingspanMean.ele$dur)
# dat.WingspanMean.ele$shift <- scale(dat.WingspanMean.ele$shift)
# dat.WingspanMean.ele$v.ele.mean <- scale(dat.WingspanMean.ele$v.ele.mean)

# #----------------------------------
# ##select only classes with >5 species
# dat.WingspanMean.lat <- data.frame(lat[which(is.na(lat$WingspanMean)==F),])
# dsp <- unique(dat.WingspanMean.lat[,c(57,68)])
# cl <- table(dsp$class)
# nn <- names(which(cl >5))
# dat.WingspanMean.lat <- dat.WingspanMean.lat[which(dat.WingspanMean.lat$class %in% nn),]
# dat.WingspanMean.lat <-droplevels(dat.WingspanMean.lat)
# 
# #scale variables
# dat.WingspanMean.lat$WingspanMean <- scale(dat.WingspanMean.lat$WingspanMean)
# dat.WingspanMean.lat$n <- scale(dat.WingspanMean.lat$n)
# dat.WingspanMean.lat$id.area <- scale(dat.WingspanMean.lat$id.area)
# dat.WingspanMean.lat$start <- scale(dat.WingspanMean.lat$start)
# dat.WingspanMean.lat$dur <- scale(dat.WingspanMean.lat$dur)
# dat.WingspanMean.lat$shift <- scale(dat.WingspanMean.lat$shift)
# dat.WingspanMean.lat$v.lat.mean <- scale(dat.WingspanMean.lat$v.lat.mean)

#----------------------------------
#per parameter
#-----------------------------------
#----------------------------------
lat_leVH <- lat_le[which(is.na(lat_le$VegetativeHeight)==F),]
dsp <- unique(lat_leVH[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
lat_leVH <- lat_leVH[which(lat_leVH$class %in% nn),]
lat_leVH <-droplevels(lat_leVH)

#scale variables
lat_leVH$VegetativeHeight <- scale(lat_leVH$VegetativeHeight)
lat_leVH$n <- scale(lat_leVH$n)
lat_leVH$id.area <- scale(lat_leVH$id.area)
lat_leVH$start <- scale(lat_leVH$start)
lat_leVH$dur <- scale(lat_leVH$dur)
lat_leVH$shift <- scale(lat_leVH$shift)
lat_leVH$v.lat.mean <- scale(lat_leVH$v.lat.mean)

#----------------------------------
lat_teVH <- lat_te[which(is.na(lat_te$VegetativeHeight)==F),]
dsp <- unique(lat_teVH[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
lat_teVH <- lat_teVH[which(lat_teVH$class %in% nn),]
lat_teVH <-droplevels(lat_teVH)

#scale variables
lat_teVH$VegetativeHeight <- scale(lat_teVH$VegetativeHeight)
lat_teVH$n <- scale(lat_teVH$n)
lat_teVH$id.area <- scale(lat_teVH$id.area)
lat_teVH$start <- scale(lat_teVH$start)
lat_teVH$dur <- scale(lat_teVH$dur)
lat_teVH$shift <- scale(lat_teVH$shift)
lat_teVH$v.lat.mean <- scale(lat_teVH$v.lat.mean)

#----------------------------------
lat_oVH <- lat_o[which(is.na(lat_o$VegetativeHeight)==F),]
dsp <- unique(lat_oVH[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
lat_oVH <- lat_oVH[which(lat_oVH$class %in% nn),]
lat_oVH <-droplevels(lat_oVH)

#scale variables
lat_oVH$VegetativeHeight <- scale(lat_oVH$VegetativeHeight)
lat_oVH$n <- scale(lat_oVH$n)
lat_oVH$id.area <- scale(lat_oVH$id.area)
lat_oVH$start <- scale(lat_oVH$start)
lat_oVH$dur <- scale(lat_oVH$dur)
lat_oVH$shift <- scale(lat_oVH$shift)
lat_oVH$v.lat.mean <- scale(lat_oVH$v.lat.mean)

#----------------------------------
ele_teVH <- ele_te[which(is.na(ele_te$VegetativeHeight)==F),]
dsp <- unique(ele_teVH[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
ele_teVH <- ele_teVH[which(ele_teVH$class %in% nn),]
ele_teVH <-droplevels(ele_teVH)

#scale variables
ele_teVH$VegetativeHeight <- scale(ele_teVH$VegetativeHeight)
ele_teVH$n <- scale(ele_teVH$n)
ele_teVH$id.area <- scale(ele_teVH$id.area)
ele_teVH$start <- scale(ele_teVH$start)
ele_teVH$dur <- scale(ele_teVH$dur)
ele_teVH$shift <- scale(ele_teVH$shift)
ele_teVH$v.ele.mean <- scale(ele_teVH$v.ele.mean)

#----------------------------------
ele_leVH <- ele_le[which(is.na(ele_le$VegetativeHeight)==F),]
dsp <- unique(ele_leVH[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
ele_leVH <- ele_leVH[which(ele_leVH$class %in% nn),]
ele_leVH <-droplevels(ele_leVH)

#scale variables
ele_leVH$VegetativeHeight <- scale(ele_leVH$VegetativeHeight)
ele_leVH$n <- scale(ele_leVH$n)
ele_leVH$id.area <- scale(ele_leVH$id.area)
ele_leVH$start <- scale(ele_leVH$start)
ele_leVH$dur <- scale(ele_leVH$dur)
ele_leVH$shift <- scale(ele_leVH$shift)
ele_leVH$v.ele.mean <- scale(ele_leVH$v.ele.mean)

#----------------------------------
ele_oVH <- ele_o[which(is.na(ele_o$VegetativeHeight)==F),]
dsp <- unique(ele_oVH[,c(57,66)])
cl <- table(dsp$class)
nn <- names(which(cl >5))
ele_oVH <- ele_oVH[which(ele_oVH$class %in% nn),]
ele_oVH <-droplevels(ele_oVH)

#scale variables
ele_oVH$VegetativeHeight <- scale(ele_oVH$VegetativeHeight)
ele_oVH$n <- scale(ele_oVH$n)
ele_oVH$id.area <- scale(ele_oVH$id.area)
ele_oVH$start <- scale(ele_oVH$start)
ele_oVH$dur <- scale(ele_oVH$dur)
ele_oVH$shift <- scale(ele_oVH$shift)
ele_oVH$v.ele.mean <- scale(ele_oVH$v.ele.mean)

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


#--------------------------
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

m <- lmer(lags ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=dat.VegetativeHeight.ele)
summary(m)

m <- lmer(lags ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=ele_leVH)
summary(m)

m <- lmer(lags ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=ele_teVH)
summary(m)

m <- lmer(lags ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=dat.VegetativeHeight.lat)
summary(m)

m <- lmer(lags ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=lat_leVH)
summary(m)

m <- lmer(lags ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=lat_teVH)
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

m <- lmer(shift ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=dat.VegetativeHeight.ele)
summary(m)

m <- lmer(shift ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=ele_leVH)
summary(m)

m <- lmer(shift ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=ele_teVH)#***
summary(m)

m <- lmer(shift ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=dat.VegetativeHeight.lat)
summary(m)

m <- lmer(shift ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=lat_leVH)
summary(m)

m <- lmer(shift ~ VegetativeHeight + (1|article_id)+(VegetativeHeight|class),data=lat_teVH)
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

write.csv(ci,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed.csv",row.names=F)
write.csv(dd,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random.csv",row.names=F)
write.csv(R2,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2.csv",row.names=F)

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
library(plotrix); library(randomcoloR)

dd <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random.csv")
ci <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed.csv")
R2 <- read.csv("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2.csv")

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


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelVelocity.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

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

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelBodySize.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

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

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelInteractionBodySizeVelocity.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

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


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelVelocity_ele.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

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

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelBodySize_ele.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

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

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelInteractionBodySizeVelocity_ele.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

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
