## make database with 1 clipped study polygon per species x study area
## include study polygons clipped with species ranges from sources in the following order of priority:
## IUCN, BOTW, GARD, BIEN, Fishmap, Aquamaps, Butterfly Maps, GIFT, GBIF occurrence
library(tidyverse)
library(ggplot2)
library(sf)
library(terra)

###################################################
##      READ IN ALL CLIPPED STUDY POLYGONS       ##
###################################################
## read in clipped species-specific polygons and combine
files <- list.files("data-processed/spp-specific-study-polygons", pattern = ".shp", full.names = TRUE)
sa_vec <- list.files("data-processed/spp-specific-study-polygons", pattern = ".shp")

sas <- c()
for(i in 1:length(files)) {
  temp <- read_sf(files[i]) 
  temp$ID = str_split_fixed(sa_vec[i], "\\.shp", 2)[,1]
  sas <- rbind(sas, temp)
  print(paste0("file num: ", i))
}

st_write(sas, "data-processed/spp-specific-study-polygons/clipped-study-polygons_ALL.shp", append = FALSE)

##########################################################
##      FILTER TO ONE POLYGON PER SPECIES x STUDY       ##
##########################################################
## read in all clipped study area polygons 
sas <- st_read("data-processed/spp-specific-study-polygons/clipped-study-polygons_ALL.shp")

## make column names match Bioshifts column names 
## select only necessary columns
colnames(sas)

sas <- sas %>%
  rename("sp_name_checked" = binomil, "RangeSource" = rng_src, "species_studyid" = spcs_st) %>%
  select(Name, Article, DOI, ID, sp_name_checked, RangeSource, species_studyid) 

View(sas[1:10, ])

## select one range source per species_studyid 
## list in order of highest priority to lowest:
range_sources <- c("IUCN", "BOTW", "GARD", "BIEN", "Fishmap", "Aquamaps", "Butterfly Maps", "GIFT", "GBIF occurrence")

sas <- sas %>%
  mutate(RangeSource = factor(.$RangeSource, ordered = TRUE, 
                              levels = range_sources))  %>%
  group_by(species_studyid) %>% ## group by species_studyid
  arrange(RangeSource, .by_group = TRUE) %>% ## sort from high to low priority within group
  filter(row_number() == 1) ## keep first row of each group (highest priority)

## check that there is one row per species x study area combination
length(unique(sas$species_studyid)) == nrow(sas)

## write geopackage 
st_write(sas, "data-processed/spp-specific-study-polygons/FINAL_species-specific-study-polygons.gpkg",
         driver = "gpkg", drop = FALSE)

## try reading it in
sas = st_read("data-processed/spp-specific-study-polygons/FINAL_species-specific-study-polygons.gpkg")
colnames(sas)

##########################################
##      ADD DESCRIPTIVE VARIABLES       ##
##########################################
## add same variables as are present in full study area polygons
## "LatCentDeg" "LonCentDeg" "Areakm2"    "LatExtentk" "EleExtentm" "Eleminm"  "Elemaxm"  "Elemeanm"


############################################
##      VISUALIZE AND SUMMARIZE DATA      ##
############################################




