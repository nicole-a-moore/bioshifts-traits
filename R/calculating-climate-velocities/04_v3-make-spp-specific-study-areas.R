## calculate species-specific climate velocities for species with primary range maps 
library(tidyverse) 
library(sf)
library(stringr)
library(terra)
library(tidyterra)


##########################
##      LOAD DATA       ##
##########################
## read in bioshifts study area polygons
polys <- st_read("data-processed/bioshiftsv3_studypolygons.shp")

## read in bioshifts v3
v3 = read.csv("data-raw/bioshiftsv3/BIOSHIFTS_v3.csv")
v3$species_name = str_replace_all(v3$sp_name_checked, "\\_", " ")
v3$species_name

## make sure all study areas have polygons
length(unique(v3$ID)) # 319
length(unique(polys$Name)) # 344
length(which(!unique(v3$ID) %in% unique(polys$Name))) # 0

## create unique ID for each species x study polygon combination
id <- v3 %>%
  select(species_name, ID) %>%
  distinct() %>%
  mutate(id_cv = paste(species_name, ID, sep = "_"))

length(id$id_cv) # 24005

## read in collated range maps
ranges <- readRDS("data-processed/large-data/collated-ranges.rds")
nrow(ranges)

##############################################
##       CROP STUDY AREA BY RANGE MAP       ##
##############################################
id <- id %>%
  filter(species_name %in% ranges$binomial) 

polys <- filter(polys, Name %in% id$ID)

length(unique(id$id_cv)) # 15926 unique species-study polygon combos

## for each study area 
sf_use_s2(FALSE)
cropped_polys = NULL
sa = 1
while(sa <= length(unique(id$ID))) {
  curr_sa = unique(id$ID)[sa]
  poly <- filter(polys, Name == curr_sa)
  
  ## get list of species in study
  species_list <- filter(id, ID == curr_sa)
  
  sp = 1
  while(sp <= length(unique(species_list$species_name))) {
    species = unique(species_list$species_name)[sp]
    
    ## crop study polygon by species range map
    range <- filter(ranges, binomial == species)
    
    r = 1
    while(r <= nrow(range)) {
      
      # ggplot(range[r,]) +
      #   geom_sf(fill = "orange")
      
      crop <- st_intersection(st_make_valid(poly), st_make_valid(range[r,])) 
      
      if(nrow(crop) == 0) {
        
      }
      else if(st_geometry_type(crop) == "GEOMETRYCOLLECTION") {
        
        
        # note: shorebirds become weird - e.g., Thalasseus sandvicensis in study area 50 (A141_P1)
        # fix later
        # ggplot(poly) +
        #   geom_sf(fill = "orange") +
        #   geom_sf(data = range, fill = "blue") +
        #   geom_sf(data = crop, fill = "red")
      }
      else {
        if(is.null(cropped_polys)) {
          cropped_polys = crop %>%
            mutate(species_studyid = paste(species, curr_sa, sep = "_"), 
                   range_source = range$range_source[r])
        } else {
          cropped_polys = crop %>%
            mutate(species_studyid = paste(species, curr_sa, sep = "_"),
                   range_source = range$range_source[r]) %>%
            rbind(cropped_polys, .)
        }
      }
      
      r = r + 1
    }
    
    print(paste0("On species no. ", sp, " of study no. ", sa))
    sp = sp + 1
  }
  
  ## save the cropped polygons for each study area in a separate file 
  if(!is.null(cropped_polys)) {
    ## make sure no duplicates
    cropped_polys = unique(cropped_polys)
    
    st_write(cropped_polys, paste0("data-processed/spp-specific-study-polygons/", curr_sa, ".shp"), 
             append = FALSE)
  }
  
  ## empty object
  cropped_polys <- NULL
  
  sa = sa + 1
}
