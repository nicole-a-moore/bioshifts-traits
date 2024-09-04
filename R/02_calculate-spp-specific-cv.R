## calculate species-specific climate velocities for species with primary range maps 
library(tidyverse)
library(sf)


###############################################
##       BIOSHIFTS STUDY AREA POLYGONS       ##
###############################################
## get bioshifts study polygons and match to species in the stud
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

## read in bioshifts 
v1 = read.table("data-raw/bioshiftsv1/Shifts2018_checkedtaxo.txt",
                header = T,
                encoding="latin1") 
v1$species_name = str_replace_all(v1$New_name, "\\_", " ")
v1$species_name

## make sure all study areas have polygons
length(unique(v1$ID)) # 325
length(unique(polys$ID)) # 337
length(which(!unique(v1$ID) %in% unique(polys$ID))) # 0

## create unique ID for each species x study polygon combination
id <- v1 %>%
  select(species_name, ID) %>%
  distinct() %>%
  mutate(id_cv = paste(species_name, ID, sep = "_"))

length(id$id_cv) # 23385


############################
##       RANGE MAPS       ##
############################
## make one object that has range maps from all sources 

## IUCN
iucn <- st_read("data-processed/large-data/IUCN/IUCN_bioshifts.shp")
iucn <- iucn[, which(colnames(iucn) %in% c("binomial", "geometry"))]
iucn$range_source = "IUCN"

## birds of the world 
botw <- st_read("data-processed/large-data/BirdsOfTheWorld/BOTW_filtered.shp")

botw <- botw[, which(colnames(botw) %in% c("sci_nam", "geometry"))]
botw <- rename(botw, "binomial" = sci_nam)
botw$range_source = "BOTW"

## GARD
gard <- st_read("data-processed/large-data/GARD/modeled_reptiles_filtered.shp")

gard <- gard[, which(colnames(gard) %in% c("Binomial", "geometry"))]
gard <- rename(gard, "binomial" = Binomial)
gard$range_source = "GARD"

## combine
ranges <- rbind(iucn, gard, botw)

length(unique(ranges$binomial)) # 2043 spp

## what kinds of species are they?
v1 %>%
  select(species_name, Class) %>%
  distinct() %>%
  filter(species_name %in% ranges$binomial) %>%
  ggplot(aes(x = Class)) +
  geom_bar() +
  coord_flip()

## TO DO
## make plot of bioshifts v1 species with & without range maps 
## - which groups are we missing maps for?





