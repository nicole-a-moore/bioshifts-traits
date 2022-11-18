## seeing whether empirical dispersal scale explains variation in range shifts
library(tidyverse)

## read in dispersal scale data 
dscale = read.csv("data-processed/dispersal-distance-collated_temp.csv")

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


#### get shift data for species with dispersal scale ####
## read in bioshifts v1
v1 = read.table("data-raw/bioshiftsv1/Shifts2018_checkedtaxo.txt",
                header = T,
                encoding="latin1") 
## clean names to make them match reported names in species list
v1$reported_name = Clean_Names(gsub("_", " ", v1$Publi_name), return_gen_sps = F)

spv1 <- filter(sp, v1 == 1) %>%
  select(reported_name, scientificName)

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


### relate velocity of shift to dispersal scale 
colnames(v1)
colnames(v2)



