## comparing IUCN range maps to study areas 
library(tidyverse)
library(sf)
library(stringr)
## start with plants

## start with first half
plants1 <- st_read("data-raw/large-data/PLANTS/PLANTS_PART1.shp")
plants1$sci_name

## filter to ones overlapping with bioshifts 
## (do taxonomic harmonization later)
v1 = read.table("data-raw/bioshiftsv1/Shifts2018_checkedtaxo.txt",
                header = T,
                encoding="latin1") 

v1$species_name = str_replace_all(v1$New_name, "\\_", " ")
v1$species_name

## filter to only latitudinal ranges 
v1 <- filter(v1, Type == "LAT")

length(which(unique(v1$species_name) %in% plants1$sci_name)) ## 67 unique lat bioshifts spp in range maps

## filter 
plants1 <- plants1[which(plants1$sci_name %in% v1$species_name),]

## now second half 
plants2 <- st_read("data-raw/large-data/PLANTS/PLANTS_PART2.shp")
plants2$sci_name

length(which(unique(v1$species_name) %in% plants2$sci_name)) ## another 51 unique lat bioshifts spp 

## filter 
plants2 <- plants2[which(plants2$sci_name %in% v1$species_name),]

## combine 
plants = rbind(plants1, plants2)
remove(plants1)
remove(plants2)

length(unique(plants$sci_name)) ## 118 

plants = st_transform(plants, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


## now, get study polygons and match to species
layers <- st_layers("data-raw/bioshifts-download/Bioshifts/Study_Areas.gdb")

sf_use_s2(FALSE)

polys <- c()
names = c()
i=1
while(i <= length(layers$name)) {
  
  curr = st_read("data-raw/bioshifts-download/Bioshifts/Study_Areas.gdb", layer = layers$name[i]) %>%
    select(Shape)
  
  names = append(names, rep(layers$name[i], nrow(curr)))
  
  polys <- rbind(polys, curr)
  
  i=i+1
}

polys$ID = names


## left join study areas with plant polygons based on species in each study 
key = select(v1, species_name, ID, Type, Param) %>% distinct()

polys_join <-left_join(key, polys) %>%
  filter(species_name %in% plants$sci_name)



## plot a study area
library(rnaturalearth)
worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')

# 
# ## for each:
# ## crop map to species range
# ## plot study area on top
# curr = polys_join %>%
#   filter(ID == unique(polys_join$ID)[8]) 
# 
# bbox = st_bbox(curr$Shape)
# 
# cropped_map <- st_crop(worldmap, c(bbox$xmin - 10, bbox$ymin - 10, 
#                                    bbox$xmax + 10, bbox$ymax + 10))
# 
# curr %>%
#   ggplot(data = ., aes(geometry = Shape)) +
#   geom_sf(data = cropped_map, aes(geometry = geometry), inherit.aes = FALSE) + 
#   geom_sf(fill = "red") +
#   theme_bw() + 
#   theme(panel.grid = element_blank())
# 
# 
# 
# ## iucn
# curr = plants %>%
#   filter(sci_name == unique(plants$sci_name)[8]) 
# 
# bbox = st_bbox(curr$geometry)
# 
# cropped_map <- st_crop(worldmap, c(bbox$xmin - 10, bbox$ymin - 10, 
#                                    bbox$xmax + 10, bbox$ymax + 10))
# 
# curr %>%
#   ggplot(data = ., aes(geometry = geometry)) +
#   geom_sf(data = cropped_map, aes(geometry = geometry), inherit.aes = FALSE) + 
#   geom_sf(fill = "red") +
#   theme_bw() + 
#   theme(panel.grid = element_blank())







## now together 
all <- left_join(polys_join, plants, by = c("species_name" = "sci_name"))

#saveRDS(all, "data-processed/temp_plantpolys.rds")
all <- readRDS("data-processed/temp_plantpolys.rds")

## pick out one species 
curr = all %>%
  filter(species_name == unique(plants$sci_name)[115]) 

## pick out one study
curr = curr %>%
  filter(ID == unique(curr$ID)[1]) %>%
  distinct()

plot_sp(curr)
unique(curr$Param)

## combine the multipolygons for each species together 

plot_sp <- function(data) {
  test = data %>%
    group_by(species_name) %>% 
    mutate(geometry = st_combine(geometry)) %>%
    select(species_name, ID, Shape, geometry) %>%
    distinct()
  
  bbox = st_bbox(test$geometry)
  cropped_map <- st_crop(worldmap, c(bbox$xmin - 10, bbox$ymin - 10, 
                                     bbox$xmax + 10, bbox$ymax + 10))
  
  test %>%
    ggplot(data = ., aes(geometry = geometry)) + 
    geom_sf(data = cropped_map, aes(geometry = geometry), inherit.aes = FALSE) + 
    geom_sf(fill = "blue", alpha = 0.5) + # range map 
    geom_sf(aes(geometry = Shape), fill = "red", alpha = 0.25) + # study area 
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    labs(title = test$species_name)
}



# Ulex europaeus studied in its invasive range in America? sp 6
## article A138_P1 seems to include many invasives from europe?


test$species_name



## plot one study + all the species in the study's maps

## pick out one study
curr = all %>%
  filter(ID == "A234_P1") %>%
  distinct()

unique(curr$Param)
## all leading edge 

## combine the multipolygons for each species together 
test = curr %>%
  group_by(species_name) %>% 
  mutate(geometry = st_combine(geometry)) %>%
  select(species_name, ID, Shape, geometry) %>%
  distinct()

bbox = st_bbox(test$geometry)
cropped_map <- st_crop(worldmap, c(bbox$xmin - 10, bbox$ymin - 10, 
                                   bbox$xmax + 10, bbox$ymax + 10))

test %>%
  ggplot(data = ., aes(geometry = geometry)) + 
  geom_sf(data = cropped_map, aes(geometry = geometry), inherit.aes = FALSE) + 
  geom_sf(fill = "blue", alpha = 0.5) + # range map 
  geom_sf(aes(geometry = Shape), fill = "red", alpha = 0.25) + # study area 
  theme_bw() + 
  theme(panel.grid = element_blank())

unique(test$species_name)






## pick out another study
curr = all %>%
  filter(ID == "A32_P1") %>%
  distinct()

## pick out one species 
curr = curr %>%
  filter(species_name == unique(curr$species_name)[1]) 

unique(curr$Param)
## all leading edge 

## combine the multipolygons for each species together 
test = curr %>%
  group_by(species_name) %>% 
  mutate(geometry = st_combine(geometry)) %>%
  select(species_name, ID, Shape, geometry) %>%
  distinct()

bbox = st_bbox(test$geometry)
cropped_map <- st_crop(worldmap, c(bbox$xmin - 10, bbox$ymin - 10, 
                                   bbox$xmax + 10, bbox$ymax + 10))

test %>%
  ggplot(data = ., aes(geometry = geometry)) + 
  geom_sf(data = cropped_map, aes(geometry = geometry), inherit.aes = FALSE) + 
  geom_sf(fill = "blue", alpha = 0.5) + # range map 
  geom_sf(aes(geometry = Shape), fill = "red", alpha = 0.25) + # study area 
  theme_bw() + 
  theme(panel.grid = element_blank()) 


