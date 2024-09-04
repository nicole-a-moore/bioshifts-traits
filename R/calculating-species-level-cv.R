## calculate ratio of overlap between species' primary range maps and study polygons  
## calculate a species-level climate velocity 
library(tidyverse)
library(sf)
library(stringr)

###############################################
##       BIOSHIFTS STUDY AREA POLYGONS       ##
###############################################
## read in bioshifts 
v1 = read.table("data-raw/bioshiftsv1/Shifts2018_checkedtaxo.txt",
                header = T,
                encoding="latin1") 
v1$species_name = str_replace_all(v1$New_name, "\\_", " ")
v1$species_name

## filter to only species with latitudinal range shifts  
v1 <- filter(v1, Type == "LAT")

## now, get bioshifts study polygons and match to species
layers <- st_layers("data-raw/bioshifts-download/Bioshifts/Study_Areas.gdb")
sf_use_s2(FALSE)

polys <- c()
names = c()
i=1
while(i <= length(layers$name)) {
  
  curr = st_read("data-raw/bioshifts-download/Bioshifts/Study_Areas.gdb", layer = layers$name[i]) %>%
    select(Shape) %>%
    st_make_valid(.) 
  
  names = append(names, rep(layers$name[i], nrow(curr)))
  
  polys <- rbind(polys, curr)
  
  i=i+1
}
polys$ID = names

##########################################
##          IUCN RANGE MAPS             ##
##########################################
iucn <- c()

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
    
    ## join to larger object
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

## transform
iucn <- st_transform(iucn, crs = st_crs(polys))

## join study area polys to important parts of v1 
key = select(v1, species_name, ID, Type, Param) %>% distinct()

polys_join <-left_join(key, polys) %>%
  filter(species_name %in% iucn$binomial) %>%
  st_as_sf(.)

## now all together 
all <- left_join(polys_join, iucn, by = c("species_name" = "binomial"))


#saveRDS(all, "data-processed/iucn_ranges_v1_mammals.rds")
all <- readRDS("data-processed/iucn_ranges_v1_mammals.rds")


###############################################
##      PLOT RANGE MAPS AND STUDY AREAS      ##
###############################################
## plot a study area
library(rnaturalearth)
worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')

## write plotting function
plot_sp <- function(data, studyid_species) {
  test = data %>%
    group_by(studyid_species) %>% 
    mutate(geometry = st_combine(geometry)) %>%
    select(studyid_species, ID, Shape, geometry) %>%
    distinct()
  
  bbox = st_bbox(test$geometry)
  cropped_map <- st_crop(worldmap, c(bbox$xmin - 10, bbox$ymin - 10, 
                                     bbox$xmax + 10, bbox$ymax + 10))
  ## plot 
  ## range map = blue, study polygon = red
  
  sp_range <- select(test, studyid_species, ID, geometry) %>% distinct()
  study_poly <- select(test, studyid_species, ID, Shape) %>% distinct()
  
  plot = sp_range %>%
    ggplot(data = ., aes(geometry = geometry)) + 
    geom_sf(data = cropped_map, aes(geometry = geometry), inherit.aes = FALSE) + 
    geom_sf(fill = "blue", alpha = 0.1) + # range map 
    geom_sf(data = study_poly, 
            aes(geometry = Shape), fill = "red", alpha = 0.25,
            inherit.aes = FALSE) + # study area 
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    labs(title = test$studyid_species)
  
  ggsave(plot, path = "figures/range-polygon-study-area-comparison",
         filename = paste0(studyid_species, ".png"),
         width = 5, height = 4)
}

## save plot for each species x study area 
## calculate ratio of overlapping area:non overlapping area
## (crop study poly by species range, crop species range by study poly + calculate areas)


## pick out one species 
curr = all %>%
  filter(species_name == unique(iucn$binomial)[1]) 

## pick out one study
curr = curr %>%
  filter(ID == unique(curr$ID)[1]) %>%
  distinct()

plot_sp(curr, curr$species_name)

## combine the multipolygons for each species together 


## loop through species-study pairs:
all$studyid_species <- paste(all$ID, all$species_name, sep = "_")

for(i in 1:length(unique(all$studyid_species))) {
  curr = all %>%
    filter(studyid_species == unique(all$studyid_species)[i]) 
  
  plot_sp(curr, unique(curr$studyid_species))
}


curr = all %>%
  filter(species_name == "Triturus cristatus") 
data=curr
studyid_species = unique(curr$studyid_species)



#################################
##      CALCULATE OVERLAP      ##
#################################

polys_join
iucn$species_name = iucn$binomial

iucn  <- iucn %>%
  filter(legend == "Extant (resident)") 

## for each unique species study pair:
polys_join$studyid_species <- paste(polys_join$ID, polys_join$species_name, sep = "_")

i=1
while(i <= length(unique(polys_join$studyid_species))) {
  studyid_species = polys_join$studyid_species[i]
  species = polys_join$species_name[i]
  
  curr = filter(polys_join, studyid_species == unique(polys_join$studyid_species)[i]) %>%
    st_union(.)
  
  ## get the species IUCN range
  iucn_range <- iucn[which(iucn$species_name == species),] %>%
    st_union(.)
  
  ## calculate overlap
  intersection <- st_intersection(curr, iucn_range)

  a_intersection <- as.numeric(st_area(intersection))
  
  if(length(nchar(a_intersection)) == 0) {
    a_intersection = 0
  }
  a_range <- as.numeric(st_area(iucn_range))
  a_study <- as.numeric(st_area(curr))
  
  # study_area_nonoverlap <- a_study - a_intersection 
  # range_nonoverlap <- a_range - a_intersection 
  
  if(i == 1) {
    df <- data.frame("studyid_species" = studyid_species,
                               "range_area_m2" = a_range,
                               "study_area_m2" = a_study,
                               "overlap_area_m2" = as.numeric(a_intersection))
  } 
  else {
    df <- rbind(df, data.frame("studyid_species" = studyid_species,
                               "range_area_m2" = a_range,
                               "study_area_m2" = a_study,
                               "overlap_area_m2" = a_intersection))
  }
  
  i = i + 1
}

df$study_area_nonoverlap <- df$study_area_m2 - df$overlap_area_m2
df$range_nonoverlap <- df$range_area_m2 - df$overlap_area_m2

## calculate proportion of sp range that is not covered
df$prop_range_studied <- df$overlap_area_m2/df$range_area_m2
  
hist(df$prop_range_studied)
  

## for the precision of climate velocity issue, what matters is the size of the species range inside the study area relative to the study area size 

df$prop_inside_study_area <- df$overlap_area_m2 / df$study_area_m2

hist(df$prop_inside_study_area)
## lots are close to 0 meaning that often species occupy a small part of the study area 


