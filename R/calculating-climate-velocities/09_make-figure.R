## plot the species-specific study area polygons on top of each other
library(tidyverse)
library(ggplot2)
library(sf)
library(terra)


sas <- st_read("data-processed/spp-specific-study-polygons/clipped-study-polygons_ALL.shp")

#####################################################
##      RASTERIZE SP SPECIFIC STUDY POLYGONS       ##
#####################################################
## for each study, rasterize each species' clipped range into a binary presence / absence grid

## loop through studies
study = 1 
while(study <= length(unique(sas$ID))) {
  ## filter to one study
  sa <- filter(sas, ID == unique(sas$ID)[study])
  
  ## for each species within the study, rasterize clipped range and add to raster representing sum
  i=1
  while(i <= nrow(sa)) {
    ## rasterize the species-specific poly
    v <- vect(sa[i,])
    r <- rast(res=1)
    ## return cover so that even cells with low overlap are counted 
    rast = rasterize(v, r, background=0, cover = TRUE)
    ## make binary presence absence grid 
    rast[rast > 0] = 1
    
    ## add oresence absence raster to other rasters in study 
    if(i==1) {
      sum_rast <- rast
    }
    else {
      sum_rast <- mosaic(sum_rast, rast, fun = "sum")
    }
    i=i+1
  }
  names(sum_rast) <- unique(sa$ID)
  
  ## save raster
  writeRaster(sum_rast, paste0("data-processed/spp-specific-study-rasters/", unique(sa$ID), ".tif"), overwrite = TRUE)
  
  print(paste0("study num: ", study))
  
  study = study + 1
}

## read in all study-level rasters and sum all
files <- list.files("data-processed/spp-specific-study-rasters", full.names = TRUE)
rast_files <- list.files("data-processed/spp-specific-study-rasters")

i=1
for(i in 1:length(files)) {
  ## read in the study-level richness raster 
  temp <- rast(files[i]) 
  
  ## sum to rasters from other studies 
  if(i==1) {
    sum_rast = temp 
  }
  else {
    sum_rast <- sum(sum_rast, temp)
  }
  print(paste0("study num: ", i))
}
sum_rast[sum_rast == 0] <- NA
plot(sum_rast)

sum_rast_copy = sum_rast
sum_rast_copy[sum_rast <= 5] <- NA
plot(sum_rast_copy)

#####################################################
##      RASTERIZE NON-SPECIFIC STUDY POLYGONS      ##
#####################################################
## IMPORTANT NOTE: need to add 1 to each cell for every species in each study without a species specific shape file




write
sas %>%
  filter(ID == "A10_P1") %>%
  ggplot(., aes()) + 
  theme_minimal() + 
  geom_sf(aes(), size = 0.5, colour = "transparent", fill = "skyblue", alpha = 0.01) +
  theme(legend.position = "none")

sa %>%
  ggplot(countries, aes(x = long, y = lat, group = group)) + 
  theme_minimal() + 
  geom_polygon(fill = "grey") +
  geom_polygon(data = ., aes())





