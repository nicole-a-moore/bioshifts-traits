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

#remove freshwater fishes
trait <- trait[-which(trait$class == "Actinopterygii" & trait$eco=="T"),]

#----------------------
#keep only marine
trait <- trait[which(trait$eco=="M"),]

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

#latitude
lat <- trait[trait$type == "LAT",]
lat$lags <- lat$shift - lat$v.lat.mean
hist(lat$lags)

#add body size
lat$BodyLengthMean <- BodyLengthMean$BodyLength[match(lat$SpeciesChecked,BodyLengthMean$SpeciesChecked)]


#prepare datasets
#----------------------------------
##select only classes with >5 species

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

#-------------------------------------------
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

#cannot run as only fish in the model
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

#---------------------------------------------
#concatenate datasets
#--------------------------------------------
R2 <- rbind(R2.mod.full.lat.LE,R2.mod.full.lat.TE,R2.mod.full.lat.O)
R2 <- data.frame(R2,Gradient = c(rep("LAT",3)),Type = rep(c("LE","TE","O"),1))


dd <- rbind(ddmod.full.lat.LE,ddmod.full.lat.TE,ddmod.full.lat.O)
dd <- data.frame(dd,Gradient = c(rep("LAT",nrow(ddmod.full.lat.LE)+nrow(ddmod.full.lat.TE)+nrow(ddmod.full.lat.O))),
                 Type = c(rep("LE",nrow(ddmod.full.lat.LE)),rep("TE",nrow(ddmod.full.lat.TE)),rep("O",nrow(ddmod.full.lat.O)))
                 
)
ci <- rbind(cimod.full.lat.LE,cimod.full.lat.TE,cimod.full.lat.O)

ci <- data.frame(ci,Variable = rownames(ci),Gradient = c(rep("LAT",nrow(cimod.full.lat.LE)+nrow(cimod.full.lat.TE)+nrow(cimod.full.lat.O))),
                 Type = c(rep("LE",nrow(cimod.full.lat.LE)),rep("TE",nrow(cimod.full.lat.TE)),rep("O",nrow(cimod.full.lat.O)))
)

ci <- na.omit(ci)

write.csv(ci,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_fixed_allgps_full_M.csv",row.names=F)
write.csv(dd,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/Coefficents_ModelsBodySize_random_allgps_full_M.csv",row.names=F)
write.csv(R2,"G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/R2_allgps_full_M.csv",row.names=F)

#---------------------------------------------------------------------------------------------
library(plotrix); library(randomcoloR)

n <- length(factor(c("Aal",sort(unique(as.character(dd$grp[dd$grpvar=="class"]))))))
palette <- data.frame(col=distinctColorPalette(n))
palette$grp <- factor(c("Aal",sort(unique(as.character(dd$grp[dd$grpvar=="class"])))))

#----------------------------------
#latitude
#----------------------------------
#velocity
fixed <- ci[ci$Variable == 'v.lat.mean',c(1:3,6)]
fixed$class <- "Aal"

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
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]


jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelVelocity_allgps_full_M.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col,legend=palette$grp,bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))

mtext(outer=T,side=3,"Velocity",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()



#body size
fixed <- ci[ci$Variable == 'BodyLengthMean' & ci$Gradient=="LAT",c(1:3,6)]
fixed$class <- "Aal"
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
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelBodySize_allgps_full_M.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col,legend=palette$grp,bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))
mtext(outer=T,side=3,"BodySize",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()


#interaction body size:velocity
fixed <- ci[ci$Variable == 'v.lat.mean:BodyLengthMean' & ci$Gradient=="LAT",c(1:3,6)]
fixed$class <- "Aal"
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
datplot <- datplot[order(datplot$grp),]
colo <- palette[match(datplot$grp,palette$grp),1]

jpeg("G:/My Drive/BIOSHIFTS Working Group/Preliminary_analyses/traits_coverage/ModelInteractionBodySizeVelocity_allgps_full_M.jpeg",width = 16*900/cm(1),height = 10*900/cm(1),res=900,pointsize=6)

plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=16,cex=2,col=ifelse(datplot$Type=="O",colo,"transparent"),xlab="",ylab="Standardized coefficients",cex.lab=1.3,axes=F)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=15,cex=2,col=ifelse(datplot$Type=="LE",colo,"transparent"),add=T)
plotCI(x = 1:nrow(datplot),y = datplot$condval,li=datplot$lwr,ui=datplot$upr,pch=17,cex=2,col=ifelse(datplot$Type=="TE",colo,"transparent"),add=T)

axis(2,tck=-0.02)
box()
abline(h=0,lty=2)

par(xpd=NA)
legend("bottomleft",col=palette$col,legend=palette$grp,bty="n",pch=16,cex=1.1,ncol=8,inset=c(0,-0.15))

mtext(outer=T,side=3,"Velocity:BodySize",cex=1.3,line=-3,font=3,adj=0.1)
dev.off()

