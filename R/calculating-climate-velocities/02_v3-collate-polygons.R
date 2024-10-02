## collate study area polygons + species range map polygons from different sources
library(tidyverse) 
library(sf)
library(stringr)
library(terra)
library(tidyterra)


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
bien$binomial <- str_replace_all(bien$binomial, "\\_", " ")
bien$range_source = "BIEN"

## Fishmap
fshmap <- st_read("data-processed/large-data/Fishmap/Fishmap_v3.shp")

fshmap <- fshmap[, which(colnames(fshmap) %in% c("SCIENTIFIC", "geometry"))]
fshmap <- rename(fshmap, "binomial" = SCIENTIFIC)
fshmap$range_source = "Fishmap"

## Butterfly Maps 
bfs <- st_read("data-processed/large-data/ButterflyMaps/ButterflyMaps_v3.shp")

bfs <- bfs[, which(colnames(bfs) %in% c("binomial", "geometry"))]
bfs$range_source = "Butterfly Maps"


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
  nrow <- append(nrow, nrow(sf))
  
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
saveRDS(gbif, "data-processed/large-data/GBIF_polygons/GBIF_v3.rds")

gbif <- readRDS("data-processed/large-data/GBIF_polygons/GBIF_v3.rds")
gbif <- gbif[, which(colnames(gbif) %in% c("spcs_nm", "geometry"))]
gbif <- rename(gbif, "binomial" = spcs_nm)
gbif$range_source = "GBIF occurrence"

## combine
ranges <- rbind(iucn, gard, botw, bien, fshmap, bfs, gbif)
length(unique(ranges$binomial)) # 9363 spp

## for now:
# files = list.files(path = "/Volumes/NIKKI/Bioshfits_GBIF/GBIF_data")
# gbif_list = str_split_fixed(files, ".qs", 2)[,1]
# gbif_list = str_replace_all(gbif_list, "\\_", " ")
# 
# ranges_for <- unique(c(gbif_list, ranges$binomial))
# ranges_for = str_replace_all(ranges_for, "\\_", " ")
# length(unique(ranges_for))

## what kinds of species are they?
# v3 %>%
#   filter(species_name %in% ranges) %>%
#   select(species_name, class) %>%
#   distinct() %>% 
#   select(-species_name) %>%
#   table() %>% View
# 
# v3 %>%
#   filter(species_name %in% ranges) %>%
#   select(species_name, class) %>%
#   distinct() %>% 
#   ggplot(aes(x = class)) +
#   geom_bar() +
#   coord_flip()
# 
# v3 %>%
#   filter(!species_name %in% ranges) %>%
#   select(species_name, class) %>%
#   distinct() %>% 
#   select(-species_name) %>%
#   table() %>% View
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
gift <- unique(gift)
st_write(gift, "data-processed/large-data/GIFT/GIFT_v3.shp", append = FALSE)

gift <- st_read("data-processed/large-data/GIFT/GIFT_v3.shp")
gift <- gift[, which(colnames(gift) %in% c("species", "geometry"))]
gift <- rename(gift, "binomial" = species)
gift$range_source = "GIFT"
gift = unique(gift)

## combine
ranges <- rbind(ranges, gift)
length(unique(ranges$binomial)) # 10710 spp

## try to get Aquamaps ranges for missing species 
################################
##         Aquamaps           ##
################################
am_key <- read.csv("/Volumes/NIKKI/Aquamaps/Speciesoccursum_v10_2019-2.csv")

am_key$Binomial = paste(am_key$Genus_valid, am_key$Species_valid)

## filter to bioshifts species 
bs_am <- filter(am_key, Binomial %in% v3$species_name)

## filter to missing bioshfits species 
bs_am <- bs_am[which(!bs_am$Binomial %in% ranges$binomial),]
length(unique(bs_am$Binomial)) ## 127

## read in the aquamaps 
aqua <- read.csv("/Volumes/NIKKI/Aquamaps/GTE10_HSPEC_NATIVE/hcaf_species_native_gte10.csv")

## filter to aquamaps for the missing species 
aqua <- filter(aqua, SpeciesID %in% bs_am$SpeciesID)

length(unique(bs_am$SpeciesID))
length(unique(aqua$SpeciesID)) ## note: 7 missing 

# filter to values above threshold ----------------------------------------
# usually, >0.5 occurrence probability works for the range
# inner join species to points in canadian EEZ ----------------------------
threshold <- .05
aqua_thresh <- aqua %>%
  # first, cut ranges to specified threshold
  # since we only consider ranges present above this
  filter(Probability >= threshold)

# split by species --------------------------------------------------------
aqua_split <- aqua_thresh %>%
  split(.$SpeciesID)

# convert to raster -------------------------------------------------------
ranges_split_rast <- purrr::map(
  .x = aqua_split,
  .f = ~.x %>%
    select(-SpeciesID, -CsquareCode) %>%
    relocate(CenterLong, CenterLat) %>%
    rast(crs = crs("EPSG: 4326")),
  .progress = T
)

ggplot() +
  geom_spatraster(data = ranges_split_rast[[4]] ) +
  theme_bw() +
  scale_fill_continuous(na.value = "transparent")

# convert to polygons -----------------------------------------------------
ranges_split_polygons <- purrr::map(
  .x = ranges_split_rast,
  .f = ~.x %>%
    as.polygons()
)

# ggplot() +
#   geom_spatvector(data = ranges_split_polygons[[2]], fill = "red") +
#   theme_bw() 

## convert to sf
sf_use_s2(FALSE)
ranges_split_sf <- purrr::map(
  .x = ranges_split_polygons,
  .f = ~.x %>%
    st_as_sf() %>%
    st_make_valid() %>%
    st_union()
)

ggplot(data = ranges_split_sf[[2]]) +
  geom_sf() 

## make into one sfc 
aquamaps <- st_sf(data.frame(geometry = do.call(rbind, ranges_split_sf)))
aquamaps$SpeciesID = names(ranges_split_sf)

## add species name column 
aquamaps <- left_join(aquamaps, bs_am) %>%
  select(Binomial, geometry) %>%
  rename("binomial" = Binomial) %>%
  mutate(range_source = "Aquamaps")

## change crs
st_crs(aquamaps) <- st_crs(ranges)

## save 
st_write(aquamaps, "data-processed/large-data/Aquamaps/Aquamaps_v3.shp", append = FALSE)

## add to ranges collection
aquamaps <- st_read("data-processed/large-data/Aquamaps/Aquamaps_v3.shp")
aquamaps <- aquamaps[, which(colnames(aquamaps) %in% c("binomil", "geometry"))]
aquamaps <- rename(aquamaps, "binomial" = binomil)
aquamaps$range_source = "Aquamaps"

ranges = rbind(ranges, aquamaps)
length(unique(ranges$binomial)) #10833

## save collated ranges as rds object 
saveRDS(ranges, "data-processed/large-data/collated-ranges.rds")




## make a breakdown of there they are from 
ranges_sp = data.frame(species_name = c(iucn$binomial, gard$binomial, botw$binomial, bien$binomial, fshmap$binomial,
                                        bfs$binomial,aquamaps$binomial, gift$binomial, gbif_list), 
                  range_source = c(iucn$range_source, gard$range_source, botw$range_source, bien$range_source, 
                                   fshmap$range_source, bfs$range_source,
                                   aquamaps$range_source, gift$range_source, rep("GBIF convex hull", length(gbif_list))))
length(unique(ranges_sp$species_name)) #10838
key = select(v3, species_name, class) %>%
  distinct()
ranges_sp <- left_join(ranges_sp, key)

ranges_sp %>%
  ggplot(aes(x = range_source)) +
  geom_bar() +
  coord_flip() +
  labs(y = "Number of species", x = "Range polygon source")

ranges_sp %>%
  ggplot(aes(x = class)) +
  geom_bar() +
  coord_flip() +
  labs(y = "Number of species", x = "Class")

missing_sp = filter(key, !species_name %in% ranges_sp$species_name)

missing_sp %>%
  ggplot(aes(x = class)) +
  geom_bar() +
  coord_flip() +
  labs(y = "Number of species missing", x = "Class")

## save list of missing species to search for in GBIF
write.csv(missing_sp, "data-processed/v3_spp-missing-range-polygons.csv", row.names = FALSE)


