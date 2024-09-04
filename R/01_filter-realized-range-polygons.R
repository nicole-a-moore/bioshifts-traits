## calculate ratio of overlap between species' primary range maps and study polygons  
## calculate a species-level climate velocity 
library(tidyverse)
library(sf)
library(stringr)

## get ranges of bioshifts species
## read in bioshifts 
v1 = read.table("data-raw/bioshiftsv1/Shifts2018_checkedtaxo.txt",
                header = T,
                encoding="latin1") 
v1$species_name = str_replace_all(v1$New_name, "\\_", " ")
v1$species_name

##########################################
##          IUCN RANGE MAPS             ##
##########################################
## read in IUCN maps 
folders = list.dirs("data-raw/large-data/IUCN")[-1]
## for each folder
f = 1
while(f <= length(folders)) {
  fol = folders[f]
  
  ## get list of shapefiles in folder 
  files = list.files(path = fol, pattern="shp$", recursive=TRUE, full.names=TRUE)
  
  i = 1
  while(i <= length(files)) {
    file = files[i]
    
    ## read in shapefiles 
    sf <- st_read(file)
    
    ## filter to only species in bioshifts
    if("sci_name" %in% colnames(sf)) {
      ## rename 
      sf = rename(sf, "binomial" = sci_name)
    }
    sf <- filter(sf, sf$binomial %in% v1$species_name)
  
    ## write
    new_filename = str_replace(str_split_fixed(file, ".shp", 2)[,1], "data-raw", "data-processed")
    st_write(sf, paste0(new_filename, "_filtered.shp"), append = FALSE)
             
    i = i + 1
  }
  
  f = f + 1
}

## make one huge range map shp file for all IUCN filtered ranges
iucn <- c()
## read in IUCN maps 
folders = list.dirs("data-processed/large-data/IUCN")[-1]
## for each folder
f = 1
while(f <= length(folders)) {
  fol = folders[f]
  
  ## get list of shapefiles in folder 
  files = list.files(path = fol, pattern="shp$", recursive=TRUE, full.names=TRUE)
  
  i = 1
  while(i <= length(files)) {
    file = files[i]
    
    ## read in shapefiles 
    sf <- st_read(file)
    
    ## bind
    if(is.null(iucn)) {
      iucn <- sf
    }
    else {
      iucn <- rbind(iucn, sf)
    }
    
    i = i + 1
  }
  
  f = f + 1
}

## filter to only extant ranges for now
iucn = iucn[which(iucn$legend == "Extant (resident)"),] 

## get rid of duplicated shapes 
iucn = iucn %>%
  group_by(id_no, SHAPE_Area, SHAPE_Leng) %>%
  mutate(count = 1:n()) %>%
  filter(count == 1) %>%
  ungroup() %>%
  select(-count)
  
length(unique(iucn$binomial)) ## 1067 - some species have multiple still
## combine them later?

## write
st_write(iucn, "data-processed/large-data/IUCN/IUCN_bioshifts.shp", append = FALSE)

##########################################
##         BIRDS OF THE WORLD           ##
##########################################

## read BOTW
botw = st_layers("data-raw/large-data/BirdsOfTheWorld/BOTW/BOTW.gdb")
botw = st_read("data-raw/large-data/BirdsOfTheWorld/BOTW/BOTW.gdb", "All_Species")

## filter to bioshifts species
botw_bs <- filter(botw, sci_name %in% v1$species_name)

## write
st_write(botw_bs, "data-processed/large-data/BirdsOfTheWorld/BOTW_filtered.shp")


############################
##         GARD           ##
############################

gard <- st_read("data-raw/large-data/GARD/GARDdoi_10.5061_dryad.83s7k__v1/GARD1.1_dissolved_ranges/modeled_reptiles.shp")

## filter to bioshifts species
gard_bs <- filter(gard, Binomial %in% v1$species_name)

st_write(gard_bs, "data-processed/large-data/GARD/modeled_reptiles_filtered.shp")


View(iucn[which(duplicated(iucn$binomial)),])

iucn %>% 
  filter(id_no == 76317591) %>%  View()
  ggplot(aes(fill = source)) + 
  geom_sf()

