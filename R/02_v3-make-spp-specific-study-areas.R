## calculate species-specific climate velocities for species with primary range maps 
library(tidyverse) 
library(sf)
library(stringr)


###############################################
##       BIOSHIFTS STUDY AREA POLYGONS       ##
###############################################
## get bioshifts study polygons and match to species in the study
files = list.files(path = "data-raw/bioshiftsv3/ShapefilesBioShiftsv3", 
                   pattern="shp$", recursive=TRUE, full.names=TRUE)

polys <- c()
i=1
while(i <= length(files)) {
  
  curr = st_read(files[i]) %>%
    st_make_valid(.) 
  
  polys <- rbind(polys, curr)
  
  i=i+1
}

## save 
st_write(polys, "data-processed/bioshiftsv3_studypolygons.shp", append = FALSE)
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


############################
##       RANGE MAPS       ##
############################
## make one object that has range maps from all sources 

## IUCN
iucn <- st_read("data-processed/large-data/IUCN/IUCN_bioshifts_v3.shp")

iucn <- iucn[, which(colnames(iucn) %in% c("binomial", "geometry"))]
iucn$range_source = "IUCN"

## BOTW (birds of the world)
botw <- st_read("data-processed/large-data/BirdsOfTheWorld/BOTW_filtered_v3.shp")

botw <- botw[, which(colnames(botw) %in% c("sci_nam", "geometry"))]
botw <- rename(botw, "binomial" = sci_nam)
botw$range_source = "BOTW"

## GARD
gard <- st_read("data-processed/large-data/GARD/modeled_reptiles_filtered_v3.shp")

gard <- gard[, which(colnames(gard) %in% c("Binomial", "geometry"))]
gard <- rename(gard, "binomial" = Binomial)
gard$range_source = "GARD"

## BIEN
bien <- st_read("data-processed/large-data/BIEN/BIEN_v3.shp")

bien <- bien[, which(colnames(bien) %in% c("species", "geometry"))]
bien <- rename(bien, "binomial" = species)
bien$range_source = "BIEN"

## Fishmap
fshmap <- st_read("data-processed/large-data/Fishmap/Fishmap_v3.shp")

fshmap <- fshmap[, which(colnames(fshmap) %in% c("SCIENTIFIC", "geometry"))]
fshmap <- rename(fshmap, "binomial" = SCIENTIFIC)
fshmap$range_source = "Fishmap"

## GBIF 
## read all together and save 
files = list.files(path = "data-processed/large-data/GBIF_polygons", 
                   pattern="shp$", recursive=TRUE, full.names=TRUE)
i = 1
gbif = NULL
while(i <= length(files)) {
  file = files[i]
  
  ## read in shapefiles 
  sf <- st_read(file)
  
  ## bind
  if(nrow(sf) > 0) {
    if(is.null(gbif)) {
      gbif <- sf
    }
    else {
      gbif <- rbind(gbif, sf)
    }
  }
  
  i = i + 1
}
st_write(gbif, "data-processed/large-data/GBIF/GBIF_v3.shp", append = FALSE)

gbif <- gbif[, which(colnames(gbif) %in% c("spcs_nm", "geometry"))]
gbif <- rename(gbif, "binomial" = spcs_nm)
gbif$range_source = "GBIF occurrence"

## combine
ranges <- rbind(iucn, gard, botw, bien, fshmap, gbif)
length(unique(ranges$binomial)) # 6583 spp

## for now:
# files = list.files(path = "/Volumes/NIKKI/Bioshfits_GBIF/GBIF_data")
# gbif_list = str_split_fixed(files, ".qs", 2)[,1]
# gbif_list = str_replace_all(gbif_list, "\\_", " ")
# 
# ranges_for <- unique(c(gbif_list, ranges$binomial))
# length(unique(ranges_for))

## what kinds of species are they?
v3 %>%
  filter(species_name %in% ranges) %>%
  select(species_name, class) %>%
  distinct() %>% 
  select(-species_name) %>%
  table() %>% View

v3 %>%
  filter(species_name %in% ranges) %>%
  select(species_name, class) %>%
  distinct() %>% 
  ggplot(aes(x = class)) +
  geom_bar() +
  coord_flip()

v3 %>%
  filter(!species_name %in% ranges) %>%
  select(species_name, class) %>%
  distinct() %>% 
  select(-species_name) %>%
  table() %>% View
## magnoliopsida, Liliopsida, bryopsida, 
## Insecta, arachnidia

#dups <- filter(ranges, binomial %in% ranges$binomial[which(duplicated(ranges$binomial))])
#View(st_drop_geometry(dups))
## note: some species more than 1 range polygon - but none have multiple from the same source  


## try to get GIFT ranges for missing species 
############################
##         GIFT           ##
############################
## GIFT does not have species range polyongs per se, but maps species to polygons they are present in using checklists for different regions
## keep regions where species is native, non-native or unknown 
library(GIFT)

all_sp <- GIFT_species(api = "https://gift.uni-goettingen.de/api/extended/",
                       GIFT_version = "latest")

length(which(unique(v3$species_name) %in% all_sp$work_species))## 6008 spp in bioshifts 
sp <- unique(v3$species_name)[unique(v3$species_name) %in% all_sp$work_species]

sp <- sp[which(!sp %in% ranges_for)]
length(unique(sp)) #1416 species that we don't already have ranges for 

## get their range maps
get_gift_range <- function(x) {
  
  filename = paste0("data-processed/large-data/GIFT/", str_replace_all(sp[x], " ", "\\_"), ".shp")
  
  if(!file.exists(filename)) {
    
    dist = GIFT_species_distribution(genus = str_split_fixed(sp, "\\ ", 2)[x,1],
                                     epithet = str_split_fixed(sp, "\\ ", 2)[x,2],
                                     aggregation = TRUE)  %>%
      filter(!is.na(native) |
               is.na(native) & is.na(naturalized)) # filter to native/naturalized 
    ## if native = NA and naturalized = NA, status is unknown in that polygon (KEEP)
    ## non-native and non-naturalized, denotes unstable cases where the species is in the process of becoming naturalized
    
    if(nrow(dist) > 0) {
      
      shapes <- GIFT_shapes(entity_ID = dist$entity_ID)
      
      map <- left_join(shapes, dist, by = "entity_ID")
      
      # map %>%
      #   ggplot() +
      #   geom_sf(aes(fill = work_species)) +
      #   scale_fill_brewer("Status", palette = "Set2") +
      #   theme_bw()
      
      ## combine into a multipoly
      map <- st_union(map) %>%
        st_as_sf() %>%
        mutate(species = unique(dist$work_species),
               range_source = "GIFT") %>%
        rename("geometry" = x)
      
      st_geometry(map) <- "geometry"
      
      st_write(map, filename, append = FALSE)
    }
  }
  print(paste0("On species number ", x))
}

lapply(1:length(sp), FUN = get_gift_range)

dead_files <- sp[which(sp %in% ranges_for)]
dead_files = paste(paste0("data-processed/large-data/GIFT/", str_replace_all(dead_files, " ", "\\_"), ".shx"))

i= 1
while(i <= length(dead_files)) {
  fn = dead_files[i]
  if (file.exists(fn)) {
    #Delete file if it exists
    file.remove(fn)
  }
  
  i = i + 1
}


## read all together and save 
files = list.files(path = "data-processed/large-data/GIFT", pattern="shp$", recursive=TRUE, full.names=TRUE)
i = 1
gift = NULL
while(i <= length(files)) {
  file = files[i]
  
  ## read in shapefiles 
  sf <- st_read(file)
  
  ## bind
  if(nrow(sf) > 0) {
    if(is.null(gift)) {
      gift <- sf
    }
    else {
      gift <- rbind(gift, sf)
    }
  }
  
  i = i + 1
}
st_write(gift, "data-processed/large-data/GIFT/GIFT_v3.shp", append = FALSE)

gift <- gift[, which(colnames(gift) %in% c("species", "geometry"))]
gift <- rename(gift, "binomial" = species)
gift$range_source = "GIFT"

## combine
ranges <- rbind(ranges, gift)


##############################################
##       CROP STUDY AREA BY RANGE MAP       ##
##############################################
id <- id %>%
  filter(species_name %in% ranges$binomial) 

polys <- filter(polys, Name %in% id$ID)

length(unique(id$id_cv)) # 6549 unique species-study polygon combos

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
      
      ggplot(range[r,]) +
        geom_sf(fill = "orange")
      
      crop <- st_intersection(st_make_valid(poly), st_make_valid(range[r,])) 
      
      if(nrow(crop) == 0) {
        
      }
      else if(st_geometry_type(crop) == "GEOMETRYCOLLECTION") {
        
        
        ## note: shorebirds become weird - e.g., Thalasseus sandvicensis in study area 50 (A141_P1)
        ## fix later
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
    st_write(cropped_polys, paste0("data-processed/spp-specific-study-polygons/", curr_sa, ".shp"), append = TRUE)
  }
  
  ## empty object
  cropped_polys <- NULL
  
  sa = sa + 1
}


## crop study polygon by range map of each species in study
## and save

ggplot(cropped_polys) +
  geom_sf() +
  geom_sf(data = crop, fill = "red") +
  geom_sf(range)