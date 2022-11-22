## seeing whether empirical dispersal scale explains variation in range shifts
library(tidyverse)
source("R/taxonomic-harmonization/clean_taxa_functions.R")

#----------------------
## read in dispersal scale data 
dscale = read.csv("data-processed/dispersal-distance-collated_temp.csv") %>%
  filter(!is.na(class)) ## get rid of species without class

unique(dscale$Unit)

val = dscale$DispersalDistance
new_val = as.numeric(as.character(dscale$DispersalDistance))
val[which(is.na(new_val))]
## only one that isn't a number is a range 
## get rid of it for now 

## convert all to Km
dscale <- dscale %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(DispersalDistanceKm = ifelse(Unit == "m",
                                     DispersalDistance/1000, 
                                     DispersalDistance)) %>%
  filter(!is.na(DispersalDistanceKm)) 

## plot distribution
dscale %>%
  ggplot(aes(x = log(DispersalDistanceKm), fill = class)) + geom_histogram() 

#----------------------
## read in list of all bioshifts species 
sp <- read_csv("data-raw/splist.csv")

#----------------------
#### get shift data for species with dispersal scale ####
## read in bioshifts v1
v1 = read.table("data-raw/bioshiftsv1/Shifts2018_checkedtaxo.txt",
                header = T,
                encoding="latin1") 

## clean names to make them match reported names in species list
v1$reported_name = Clean_Names(gsub("_", " ", v1$Publi_name), return_gen_sps = F)

spv1 <- filter(sp, v1 == 1) %>%
  select(reported_name, scientificName) %>%
  unique()

## yay! all are there
which(!v1$reported_name %in% spv1$reported_name)

## join to add scientific name column
v1 = left_join(v1, spv1)

## subset to species with dispersal scale 
v1 <- filter(v1, scientificName %in% dscale$scientificName)
length(unique(v1$scientificName)) #582 species 


## read in bioshifts v2
v2 = read.csv("data-raw/bioshiftsv2/Bioshifts.v2.final.csv")
## clean names to make them match reported names in species list
v2$reported_name = Clean_Names(gsub("_", " ", v2$Scientific.Name), return_gen_sps = F)

spv2 <- filter(sp, v2 == 1) %>%
  select(reported_name, scientificName)

## yay! all are there
which(!v2$reported_name %in% spv2$reported_name)

## join to add scientific name column
v2 = left_join(v2, spv2)

## subset to species with dispersal scale 
v2 <- filter(v2, scientificName %in% dscale$scientificName)
length(unique(v2$scientificName)) #541 species 

## make sure all the species are there:
length(unique(c(v2$scientificName, v1$scientificName))) # 604 unique species 
length(unique(dscale$scientificName)) # 604 unique species 

#----------------------
### relate velocity of shift to dispersal scale ###
colnames(v1)
colnames(v2)

#----------------------
## clean bioshifts v1 

#--------------------------
#use insects orders to avoid mixing butterflies with others
dscale$class[dscale$class=="Insecta"] <- dscale$order[dscale$class=="Insecta"]
v1$Class[v1$Class=="Insecta"] <- v1$Order[v1$Class=="Insecta"]

#----------------------
#correct and reorganize methodological variables
v1$Data[v1$Data=="occurence-based"] <- "occurrence-based"
v1$Sampling <- ifelse(v1$Sampling == "TWO","TWO","MULTIPLE")
v1$Grain_size <- ifelse(v1$Grain_size %in% c("large","very_large"),"large",v1$Grain_size)
v1$Uncertainty_Distribution <- ifelse(v1$Uncertainty_Distribution %in% c("RESAMPLING","RESAMPLING(same)"),"RESAMPLING",
                                         ifelse(v1$Uncertainty_Distribution %in% c("MODEL","MODEL+RESAMPLING(same)","RESAMPLING+MODEL"),"MODEL",
                                                ifelse(v1$Uncertainty_Distribution %in% c("DETECTABILITY","RESAMPLING(same)+DETECTABILITY"),"DETECTABILITY",
                                                       v1$Uncertainty_Distribution)))
v1$Uncertainty_Distribution <- ifelse(v1$Uncertainty_Distribution == "RAW","OPPORTUNISTIC","PROCESSED")  

#transform study area
v1$ID.area <- log(v1$ID.area)


#----------------------
#remove freshwater fishes
#v1 <- v1[-which(v1$Class == "Actinopterygii" & v1$ECO=="T"),]

#remove marine birds
#v1 <- v1[-which(v1$Class == "Aves" & v1$ECO=="M"),]

#----------------------
## add dispersal scale 
## get rid of old taxonomy columns from v1 (they aren't right)
v1 <- select(v1, -c("Kingdom", "Phylum", "Class", "Order", "Family"))
v1_saved = v1

## after everything, there should be:
length(which(v1$Type == "ELE")) ## 1337 elevation shifts
length(which(v1$Type == "LAT")) ## 1902 latitude shifts

## get rid of columns that will cause duplication in dispersal scale database 
dscale <- select(dscale, -c("reported_name", "reported_name_fixed", "db", "db_code")) %>%
  unique()

v1 = left_join(v1, dscale) 

## check on merge
length(which(is.na(v1$DispersalDistanceKm))) #0 missing dispersal scale
length(unique(v1$scientificName)) #still 582 species


#----------------------
## calculate one dispersal distance value per species per metric
## ex. one mean, one max
# dscale %>%
#   ggplot(aes(x = Code)) + geom_bar() + coord_flip()

v1 = v1 %>% 
  group_by(scientificName, Code) %>%
  mutate(Code = ifelse(Code %in% c("99thPercentileDispersalDistance",
                                    "90thPercentileDispersalDistance"),
                        "MaxDispersalDistance", 
                        Code)) %>% ## recode 99th and 90th percentile to be "max"
  mutate(DispersalDistanceKm_unique = ifelse(Code == "MaxDispersalDistance",
                                             max(DispersalDistanceKm),
                                             ifelse(Code == "MedianDispersalDistance",
                                                    median(DispersalDistanceKm), 
                                                    ifelse(Code %in% c("MeanDispersalDistance", "DispersalDistance"),
                                                           mean(DispersalDistanceKm),
                                                           NA)))) %>%
  ungroup() %>%
  group_by(scientificName) %>%
  mutate(DispersalDistanceKm_max = max(DispersalDistanceKm)) %>%
  select(Code, scientificName, DispersalDistanceKm, DispersalDistanceKm_unique, DispersalDistanceKm_max, everything()) %>%
  ungroup()

## look at the metrics that were larger than max 
max = v1 %>%
  filter(Code == "MaxDispersalDistance") %>%
  mutate(matches = ifelse(DispersalDistanceKm_unique == DispersalDistanceKm_max, "Y", "N")) %>%
  select(matches, scientificName) %>%
  unique() %>%
  filter(matches == "N")

length(unique(max$scientificName))
## 8 species have other dispersal distance metrics that are larger than their so-called maximum 

v1 %>%
  filter(scientificName %in% max$scientificName) %>%
  filter(DispersalDistanceKm == DispersalDistanceKm_max) %>% 
  select(scientificName, Code) %>% 
  unique() %>% View
## most are mean values, some median
## mostly fish, 1 plant, 1 tree, some small mammals

## correct max dispersal distance for these species
v1 = v1 %>%
  group_by(scientificName, Code) %>% 
  mutate(DispersalDistanceKm_unique = ifelse(DispersalDistanceKm_unique != DispersalDistanceKm_max &
                                               Code == "MaxDispersalDistance", 
                                             DispersalDistanceKm_max,
                                             DispersalDistanceKm_unique)) %>%
  ungroup() %>%
  select(-DispersalDistanceKm_max)

v1 %>%
  select(Code, scientificName) %>%
  unique() %>%
  group_by(Code) %>%
  tally()

## make sure we don't drop any data
maxsp = unique(v1$scientificName[which(v1$Code == "MaxDispersalDistance")])
nomaxsp <- unique(v1$scientificName[which(!v1$scientificName %in% maxsp)])

## 561 species have max dispersal distance
## start with max!
v1 <- v1 %>%
  filter(Code == "MaxDispersalDistance") %>%
  select(-DispersalDistanceKm, -DispersalDistance, -Unit, -Field, -Source, -Database, -ObservationType, -Sex) %>%
  unique()

which(nomaxsp %in% v1$scientificName)
which(maxsp %in% v1$scientificName)

check = v1_saved %>%
  filter(scientificName %in% maxsp)
nrow(check) == nrow(v1)
## yay = all observations are there without any duplication

v1 <- rename(v1, "MaxDispersalDistanceKm" = DispersalDistanceKm_unique)

## save dataset 
#write.csv(v1, "data-processed/bioshiftsv1_max-dispersal-distance.csv", row.names = FALSE)
v1 = read.csv("data-processed/bioshiftsv1_max-dispersal-distance.csv")

#----------------------
#take log of dispersal scale
v1$MaxDispersalDistanceKm = log(v1$MaxDispersalDistanceKm)
hist(v1$MaxDispersalDistanceKm)

#----------------------
#prepare a dataset for each gradient

#elevation
ele <- v1[v1$Type == "ELE",]
ele <- ele[which(is.na(ele$v.ele.mean)==F),]
ele$lags <- ele$SHIFT - ele$v.ele.mean
hist(ele$lags)

#latitude
lat <- v1[v1$Type == "LAT",]
lat$lags <- lat$SHIFT - lat$v.lat.mean
hist(lat$lags)

## how many shifts for species with maximum dispersal distance in v1? 
nrow(ele) + nrow(lat) #3151

## how many leading edge shifts for species with maximum dispersal distance in v1? 
length(which(ele$Param == "LE")) + length(which(lat$Param == "LE")) #973

#----------------------
##select only classes with >5 species
## elevation
dat.MaxDispDist.ele <- data.frame(ele[which(is.na(ele$MaxDispersalDistanceKm)==F),])

dsp <- unique(dat.MaxDispDist.ele[,c(which(colnames(dat.MaxDispDist.ele) == "class"),
                                     which(colnames(dat.MaxDispDist.ele) == "scientificName"))])
cl <- table(dsp$class)
nn <- names(which(cl >10))
dat.MaxDispDist.ele <- dat.MaxDispDist.ele[which(dat.MaxDispDist.ele$class %in% nn),]
dat.MaxDispDist.ele <-droplevels(dat.MaxDispDist.ele)
length(unique(dat.MaxDispDist.ele$class))

#scale variables
dat.MaxDispDist.ele$max_disp_dist <- scale(dat.MaxDispDist.ele$MaxDispersalDistanceKm)
dat.MaxDispDist.ele$n <- scale(dat.MaxDispDist.ele$N)
dat.MaxDispDist.ele$id.area <- scale(dat.MaxDispDist.ele$ID.area)
dat.MaxDispDist.ele$start <- scale(dat.MaxDispDist.ele$START)
dat.MaxDispDist.ele$dur <- scale(dat.MaxDispDist.ele$DUR)
dat.MaxDispDist.ele$shift <- scale(dat.MaxDispDist.ele$SHIFT)
dat.MaxDispDist.ele$v.ele.mean_scaled <- scale(dat.MaxDispDist.ele$v.ele.mean)

## latitude 
dat.MaxDispDist.lat <- data.frame(lat[which(is.na(lat$MaxDispersalDistanceKm)==F),])

dsp <- unique(dat.MaxDispDist.lat[,c(which(colnames(dat.MaxDispDist.lat) == "class"),
                                     which(colnames(dat.MaxDispDist.lat) == "scientificName"))])
cl <- table(dsp$class)
nn <- names(which(cl >10))
dat.MaxDispDist.lat <- dat.MaxDispDist.lat[which(dat.MaxDispDist.lat$class %in% nn),]
dat.MaxDispDist.lat <-droplevels(dat.MaxDispDist.lat)
length(unique(dat.MaxDispDist.lat$class))

#scale variables
dat.MaxDispDist.lat$max_disp_dist <- scale(dat.MaxDispDist.lat$MaxDispersalDistanceKm)
dat.MaxDispDist.lat$n <- scale(dat.MaxDispDist.lat$N)
dat.MaxDispDist.lat$id.area <- scale(dat.MaxDispDist.lat$ID.area)
dat.MaxDispDist.lat$start <- scale(dat.MaxDispDist.lat$START)
dat.MaxDispDist.lat$dur <- scale(dat.MaxDispDist.lat$DUR)
dat.MaxDispDist.lat$shift <- scale(dat.MaxDispDist.lat$SHIFT)
dat.MaxDispDist.lat$v.lat.mean_scaled <- scale(dat.MaxDispDist.lat$v.lat.mean)


#----------------------
#prepare a dataset for each range shift parameter
ele_te <- dat.MaxDispDist.ele[dat.MaxDispDist.ele$Param=="TE",]
ele_le <- dat.MaxDispDist.ele[dat.MaxDispDist.ele$Param=="LE",]
ele_o <- dat.MaxDispDist.ele[dat.MaxDispDist.ele$Param=="O",]
lat_te <- dat.MaxDispDist.lat[dat.MaxDispDist.lat$Param=="TE",]
lat_le <- dat.MaxDispDist.lat[dat.MaxDispDist.lat$Param=="LE",]
lat_o <- dat.MaxDispDist.lat[dat.MaxDispDist.lat$Param=="O",]

## simple simple model to warm up:
lm(shift ~ max_disp_dist, data = lat_le)
lm(shift ~ max_disp_dist, data = ele_le)

lat_le %>%
  ggplot(aes(x = MaxDispersalDistanceKm, y = SHIFT, colour = class)) + geom_point()
## interesting - so birds have larger dispersal distances AND more variable shifts

## what does relationship between shift and velocity look like?
lat_le %>%
  ggplot(aes(x = v.lat.mean, y = SHIFT, colour = class)) + geom_point()
## ah - same study must each have same mean velocity 

lat_le %>%
  ggplot(aes(x = v.lat.mean, y = SHIFT, colour = Article_ID)) + geom_point()
## yep

## look within a study
lat_le %>%
  filter(Article_ID == "A32") %>%
  ggplot(aes(x = MaxDispersalDistanceKm, y = SHIFT, colour = class)) + geom_point()

lat_le %>%
  filter(Article_ID == "A138") %>%
  ggplot(aes(x = MaxDispersalDistanceKm, y = SHIFT, colour = class)) + geom_point()

lat_le %>%
  filter(Article_ID == "A9") %>%
  ggplot(aes(x = MaxDispersalDistanceKm, y = SHIFT, colour = class)) + geom_point()
lat_le %>%
  filter(Article_ID == "A9") %>%
  ggplot(aes(x = MaxDispersalDistanceKm, y = lags, colour = class)) + geom_point()
## ah - good, learning a lot about the data 
## because there is one mean velocity per study, using lags doesn't change the structure of variation within studies 


m_lat = lmer(shift ~ v.lat.mean_scaled*max_disp_dist + (1|Article_ID) + (max_disp_dist|class), 
             data = lat_le)
summary(m_lat)

m_lat = lmer(shift ~ v.lat.mean_scaled*max_disp_dist + dur + id.area + n + (1|Article_ID) + (max_disp_dist|class), 
             data = lat_le)
summary(m_lat)





## things to ask: 
## is dispersal scale more important at leading edge or trailing edge?
## how much variation does dispersal scale explain at leading edge?



### thoughts:
## reproductive/dispersal frequency also matters!!!! need to know this too

## we have multiple observations of the same species, and the traits for each species are the same - include as a random effect?
## (Lise did not)
## count species within same article 
lat_le %>% 
  group_by(Article_ID, scientificName) %>%
  arrange(-n) %>%
  tally() %>% View




#----------------------
## model building notes from Lise:
# shifts ~ velocity x bodySize + StudyArea + StudyDuration + dataType + 
#   Nsamplings + dataProcessing + ( velocity x bodySize | Class) + (1|ArticleID)
#  To avoid singular fit (and shrinkage), I only selected taxonomic groups with > 10 or 15 species


## StudyArea = ID.area 
## StudyDuration = DUR
## dataType = Type
## Nsamplings = N
## dataProcessing = Uncertainty_Distribution
## Type = LE, TE, O
