## calculate ratio of overlap between species' primary range maps and study polygons  
## calculate a species-level climate velocity 
library(tidyverse)
library(sf)
library(stringr)

## get ranges of bioshifts v3 species
## read in bioshifts 
v3 = read.csv("data-raw/bioshiftsv3/BIOSHIFTS_v3.csv")
v3$species_name = str_replace_all(v3$sp_name_checked, "\\_", " ")
v3$species_name

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
    sf <- filter(sf, sf$binomial %in% v3$species_name)
  
    ## write
    new_filename = str_replace(str_split_fixed(file, ".shp", 2)[,1], "data-raw", "data-processed")
    st_write(sf, paste0(new_filename, "_filtered_v3.shp"), append = FALSE)
             
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
  
  ## get list of shapefiles in folder that are from v3
  files = list.files(path = fol, pattern="shp$", recursive=TRUE, full.names=TRUE)
  files = files[which(str_detect(files, "v3"))]
  
  i = 1
  while(i <= length(files)) {
    file = files[i]
    
    ## read in shapefiles 
    sf <- st_read(file)
    
    ## bind
    if(nrow(sf) > 0) {
      if(is.null(iucn)) {
        iucn <- sf
      }
      else {
        iucn <- rbind(iucn, sf)
      }
    }
    
    i = i + 1
  }
  
  f = f + 1
}

iucn$legend

## get rid of duplicated shapes 
iucn = iucn %>%
  group_by(id_no, SHAPE_Area, SHAPE_Leng) %>%
  mutate(count = 1:n()) %>%
  filter(count == 1) %>%
  ungroup() %>%
  select(-count)

## save object that has multiple ranges for each species (extant, extinct, breeding, non-breeding, etc.)
st_write(iucn, "data-processed/large-data/IUCN/IUCN_bioshifts_v3_notunion.shp", append = FALSE)

## keep all ranges (extant, extinct, breeding, non-breeding, etc.)
## combine all ranges for a single species into one
sf_use_s2(FALSE)

iucn <- iucn %>%
  group_by(binomial) %>% 
  summarize(geometry = st_union(st_make_valid(geometry))) %>%
  ungroup()

length(unique(iucn$binomial)) ## 1266 - now each species has only 1

## write
st_write(iucn, "data-processed/large-data/IUCN/IUCN_bioshifts_v3.shp", append = FALSE)

##########################################
##         BIRDS OF THE WORLD           ##
##########################################
## read BOTW
botw = st_layers("data-raw/large-data/BirdsOfTheWorld/BOTW/BOTW.gdb")
botw = st_read("data-raw/large-data/BirdsOfTheWorld/BOTW/BOTW.gdb", "All_Species")

## filter to bioshifts species
botw_bs <- filter(botw, sci_name %in% v3$species_name)

## save object that has multiple seasonal ranges for each species
st_write(iucn, "data-processed/large-data/BirdsOfTheWorld/BOTW_filtered_v3_notunion.shp", append = FALSE)

## combine all ranges for a single species into one
botw_bs <- botw_bs %>%
  group_by(sci_nam) %>% 
  summarize(geometry = st_union(st_make_valid(geometry))) %>%
  ungroup()

length(unique(botw_bs$sci_nam)) ## 1332 - now each species has only 1

## write
st_write(botw_bs, "data-processed/large-data/BirdsOfTheWorld/BOTW_filtered_v3.shp", append = FALSE)

############################
##         GARD           ##
############################
gard <- st_read("data-raw/large-data/GARD/GARDdoi_10.5061_dryad.83s7k__v1/GARD1.1_dissolved_ranges/modeled_reptiles.shp")

## filter to bioshifts species
gard_bs <- filter(gard, Binomial %in% v3$species_name)

st_write(gard_bs, "data-processed/large-data/GARD/modeled_reptiles_filtered_v3.shp")


############################
##         BIEN           ##
############################
library(BIEN)

## filter bioshifts to only plant species to search BIEN
plants <- filter(v3, kingdom == "Plantae")
bien <- NULL
p = 1
while(p <= length(unique(plants$species_name))) {
  
  pl_range <- BIEN_ranges_load_species(unique(plants$species_name)[p])
  
  if(nrow(pl_range) != 0) {
    if(is.null(bien)) {
      bien <- pl_range
    } 
    else {
      bien <- rbind(bien, pl_range)
    }
  }
  
  print(paste0("On plant no. ", p))
  p = p + 1
}

## save 
st_write(bien, "data-processed/large-data/BIEN/BIEN_v3.shp")


###############################
##         FISHMAP           ##
###############################
fshmap <- st_read("data-raw/large-data/Fishmap/CAAB_FISHMAP/CAAB_FISHMAPPolygon.shp")

## filter to bioshifts species
fshmap_bs <- filter(fshmap, SCIENTIFIC %in% v3$species_name)

st_write(fshmap_bs, "data-processed/large-data/Fishmap/Fishmap_v3.shp")


######################################
##          BUTTERFLIES             ##
######################################
## read all together
files = list.files(path = "data-raw/large-data/ButterflyMaps/butterfly_maps_bioshifts", 
          recursive=TRUE, full.names=TRUE)

bfs = NULL
for(i in 1:length(files)) {
  cur = st_read(files[i])
 
  if(is.null(bfs)) {
    bfs <- cur
  }
  else {
    bfs <- rbind(bfs, cur)
  }
}

## filter to bioshifts species
bfs <- filter(bfs, binomial %in% v3$species_name)

st_write(bfs, "data-processed/large-data/ButterflyMaps/ButterflyMaps_v3.shp", append = FALSE)



