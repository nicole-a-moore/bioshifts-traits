
################################################################
# RETRIEVE ACCEPTED NAME AND CLASSIFICATION FOR A LIST OF TAXA #
################################################################


dat = read.table("...")
#--> dat= import table of names to verify with 2 columns named Genus and Species

#keep only species level (remove varieties or subspecies)
species = sapply(strsplit(as.character(dat$Species),"_"),'[',1)

#identify genus level (if no species name)
species[grep("spp.",species,fixed=T)] = NA
species[grep("sp.",species,fixed=T)] = NA
species[grep("cf",species,fixed=T)] = NA
species[grep("sect.",species,fixed=T)] = NA

#create vector name
mynames <-na.omit(unique(ifelse(is.na(species),as.character(dat$Genus),paste(dat$Genus,species))))      

#check number of names
length(mynames)
mynames_species <-unique(paste(dat$Genus,species)[-which(is.na(species))])      
length(mynames_species)

#---- FUNCTIONS ####
library(rgbif)
library(taxize)

GB_classif = function(x){
  gb = name_backbone(x)
  phylum = ifelse(length(gb$phylum) > 0, gb$phylum, NA)
  class = ifelse(length(gb$class) > 0, gb$class, NA)
  order = ifelse(length(gb$order) > 0, gb$order, NA)
  family = ifelse(length(gb$family) > 0, gb$family, NA)
  genus = ifelse(length(gb$genus) > 0, gb$genus, NA)
  species = ifelse(length(gb$species) > 0, gb$species, NA)
  classif = c(phylum,class,order,family,genus,species,x)
  return(classif)
}

NCBI_classif = function(x,db){
  gb = try(classification(x,db = db, rows = 1),silent=T)
  if(class(gb) == "try-error"){
    classif = c(NA,NA,NA,NA,NA,NA,as.character(x))
    return(classif)
  } else{    
    phylum = ifelse(length(gb[[1]]) ==1, NA,as.character(gb[[1]][,1][gb[[1]][,2] %in% c('phylum','Phylum')]))
    class = ifelse(length(gb[[1]]) ==1, NA,as.character(gb[[1]][,1][gb[[1]][,2] %in% c('class','Class')]))
    order = ifelse(length(gb[[1]]) ==1, NA,as.character(gb[[1]][,1][gb[[1]][,2] %in% c('order','Order')]))
    family = ifelse(length(gb[[1]]) ==1, NA,as.character(gb[[1]][,1][gb[[1]][,2] %in% c('family','Family')]))
    genus = ifelse(length(gb[[1]]) ==1, NA,as.character(gb[[1]][,1][gb[[1]][,2] %in% c('genus','Genus')]))
    species = ifelse(length(gb[[1]]) ==1, NA,as.character(gb[[1]][,1][gb[[1]][,2] %in% c('species','Species')]))
    classif = c(phylum,class,order,family,genus,species,as.character(x))
    return(classif)
  }
}

#---- RETRIEVE NAMES AND CLASSIFICATION ####

#I - retrive names using NCBI ####
#*****************************

gbiftestout_df <- llply(mynames, NCBI_classif,db='ncbi',.progress = "text")
dfnames.ncbi  = sapply(gbiftestout_df,as.data.frame)
dfnames.ncbi  = sapply(dfnames.ncbi,unlist)
dfnames.ncbi = t(dfnames.ncbi)
dfnames.ncbi  = data.frame(dfnames.ncbi)
colnames(dfnames.ncbi) <- c("phylum", "class", "order", "family","genus","species","old_name")
dfnames.ncbi  = apply(dfnames.ncbi,2,unlist)
dfnames.ncbi  = apply(dfnames.ncbi,2,as.character)
dfnames.ncbi  = data.frame(dfnames.ncbi)

## 1 -Check for changes in species names
na.omit(cbind(as.character(dfnames.ncbi$old_name[which(!dfnames.ncbi$old_name %in% dfnames.ncbi$species)]),as.character(dfnames.ncbi$species[which(!dfnames.ncbi$old_name %in% dfnames.ncbi$species)])))

## 2 -Select unidentified species & remove from the table
unid.ncbi = dfnames.ncbi$old_name[which(sapply(strsplit(as.character(dfnames.ncbi$old_name)," "),length)== 2 & is.na(dfnames.ncbi$species))]

#3 - save file
if (length(unid.ncbi)>0){
dfnames.ncbi = dfnames.ncbi[-which(sapply(strsplit(as.character(dfnames.ncbi$old_name)," "),length)== 2 & is.na(dfnames.ncbi$species)),] 
}

dfnames.ncbi = apply(dfnames.ncbi,2,as.character)
write.table(dfnames.ncbi,"dfnames.ncbi.txt")

#4 - compile
dfnames.ncbi = read.table("dfnames.ncbi.txt",h=T)
source.ncbi = rep('ncbi',nrow(dfnames.ncbi))
taxo_level = ifelse(sapply(strsplit(as.character(dfnames.ncbi$old_name)," "),length) == 1,"Genus","Species")
finaltab = data.frame(dfnames.ncbi,source.ncbi,taxo_level)
colnames(finaltab)=c("phylum","class","order","family","genus","species","old_name","source","taxo_level")


#/!\ STOP HERE IF length(unid.ncbi)=0
#/!\ all species have been retrieved no need to go further
length(unid.ncbi)


#II - retrive names using ITIS ####
#***************************** 

gbiftestout_df.itis <- llply(unid.ncbi, NCBI_classif,db='itis',.progress = "text")

dfnames.itis  = sapply(gbiftestout_df.itis,as.data.frame)
dfnames.itis  = sapply(dfnames.itis,unlist)
dfnames.itis = t(dfnames.itis)
dfnames.itis  = data.frame(dfnames.itis)
colnames(dfnames.itis) <- c("phylum", "class", "order", "family","genus","species","old_name")
dfnames.itis  = apply(dfnames.itis,2,unlist)
dfnames.itis  = data.frame(dfnames.itis)

## 1 -Check for changes in species names
na.omit(cbind(as.character(dfnames.itis$old_name[which(!dfnames.itis$old_name %in% dfnames.itis$species)]),as.character(dfnames.itis$species[which(!dfnames.itis$old_name %in% dfnames.itis$species)])))
#
## 2 - Select unidentified species
unid.itis = dfnames.itis$old_name[which(sapply(strsplit(as.character(dfnames.itis$old_name)," "),length)== 2 & is.na(dfnames.itis$species))]
#3 - save file
if (length(unid.itis)>0){
dfnames.itis = dfnames.itis[-which(sapply(strsplit(as.character(dfnames.itis$old_name)," "),length)== 2 & is.na(dfnames.itis$species)),] 
}

write.table(dfnames.itis,"dfnames.itis.txt")

#4 - compile
dfnames.itis = read.table("dfnames.itis.txt",h=T)
source.itis = rep('itis',nrow(dfnames.itis))
taxo_level = ifelse(sapply(strsplit(as.character(dfnames.itis$old_name)," "),length) == 1,"Genus","Species")
dfnames.itis=data.frame(dfnames.itis,source.itis,taxo_level)
colnames(dfnames.itis)=c("phylum","class","order","family","genus","species","old_name","source","taxo_level")
finaltab = rbind(finaltab,dfnames.itis)

#/!\ STOP HERE IF length(unid.ncbi)= 0
#/!\ all species have been retrieved no need to go further

length(unid.itis)
 
#III - retrive names using GBIF ####
#but not the classif (names are easier to retrieve using GBIF but the clasif has many problems)
#*****************************

gbiftestout_df.gb <- llply(as.character(unid.itis), GB_classif,.progress = "text")

dfnames.gb  = sapply(gbiftestout_df.gb,as.data.frame)
dfnames.gb  = sapply(dfnames.gb,unlist)
dfnames.gb = t(dfnames.gb)
dfnames.gb  = data.frame(dfnames.gb)
colnames(dfnames.gb) <- c("phylum", "class", "order", "family","genus","species","old_name")
dfnames.gb  = apply(dfnames.gb,2,unlist)
dfnames.gb  = apply(dfnames.gb,2,as.character)
dfnames.gb  = data.frame(dfnames.gb)


## 1 - Check for changes in species names
na.omit(cbind(as.character(dfnames.gb$old_name[which(!dfnames.gb$old_name %in% dfnames.gb$species)]),as.character(dfnames.gb$species[which(!dfnames.gb$old_name %in% dfnames.gb$species)])))

## 2 - Save file
write.table(dfnames.gb,"dfnames.gb.txt")

dfnames.gb = read.table("dfnames.gb.txt",h=T)
source.gb = rep('gb',nrow(dfnames.gb))

## 3 -check the info against ncbi (with new names from GBIF)
check.ncbi <- llply(dfnames.gb$species, NCBI_classif,db='ncbi',.progress = "text")
dfnames.check.ncbi  = sapply(check.ncbi,as.data.frame)
dfnames.check.ncbi  = sapply(dfnames.check.ncbi,unlist)
dfnames.check.ncbi = t(dfnames.check.ncbi)
dfnames.check.ncbi  = data.frame(dfnames.check.ncbi)
colnames(dfnames.check.ncbi) <- c("phylum", "class", "order", "family","genus","species","old_name")
dfnames.check.ncbi  = apply(dfnames.check.ncbi,2,unlist)
dfnames.check.ncbi  = apply(dfnames.check.ncbi,2,as.character)
dfnames.check.ncbi  = data.frame(dfnames.check.ncbi)

## 4 - Check for changes in species names
na.omit(cbind(as.character(dfnames.check.ncbi$old_name[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species)]),as.character(dfnames.check.ncbi$species[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species)])))

##5 - change the names in DB
dfnames.gb = apply(dfnames.gb,2,as.character)
dfnames.gb[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),1] = as.character(dfnames.check.ncbi[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),1])
dfnames.gb[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),2] =  as.character(dfnames.check.ncbi[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),2])
dfnames.gb[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),3] =  as.character(dfnames.check.ncbi[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),3])
dfnames.gb[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),4] =  as.character(dfnames.check.ncbi[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),4])
dfnames.gb[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),5] =  as.character(dfnames.check.ncbi[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),5])
dfnames.gb[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),6] =  as.character(dfnames.check.ncbi[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F),5])

source.gb[which(!dfnames.check.ncbi$old_name %in% dfnames.check.ncbi$species & is.na(dfnames.check.ncbi$species) == F)]  = 'ncbi'

##6 -  check unidentified species  : mainly mistakes in names
dfnames.gb=as.data.frame(dfnames.gb)
dfnames.gb$old_name[which(is.na(dfnames.gb$species) & is.na(dfnames.gb$genus))]

#7 - compile
taxo_level = ifelse(sapply(strsplit(as.character(dfnames.gb$old_name)," "),length) == 1,"Genus","Species")
dfnames.gb=data.frame(dfnames.gb,source.gb,taxo_level)
colnames(dfnames.gb)=c("phylum","class","order","family","genus","species","old_name","source","taxo_level")
finaltab = rbind(finaltab,dfnames.gb)

#finaltab gives for each name entered in dat (finaltab$old_name) the full classification from phylum to species
#species=NA have not been recognized as valid species and need to be search elsewhere or removed
#classification sometimes has to be completed

