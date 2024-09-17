## make range map polygons by drawing convex hulls around filtered GBIF occurences
library(gbif.range)
library(ggplot2)
library(terra)
library(stringr)
library(tidyverse)
library(sf)

## read in bioshifts 
v3 = read.csv("data-raw/bioshiftsv3/BIOSHIFTS_v3.csv")
v3$species_name = str_replace_all(v3$sp_name_checked, "\\_", " ")
v3$species_name

## read GBIF occurrence filenames
files = list.files("/Volumes/NIKKI/Bioshfits_GBIF/GBIF_data")

## make sure all are in v3
which(!str_split_fixed(files, "\\.qs", 2)[,1] %in% v3$sp_name_checked)
## 530 and 2602 not there 
## DEAL WITH LATER

## get realm from v3
key = v3 %>%
  filter(sp_name_checked %in% str_split_fixed(files, "\\.qs", 2)[,1]) %>%
  select(sp_name_checked, Eco, species_name) %>%
  distinct()

## get the ecoregions 
eco_terra = read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)
eco_fresh = read_bioreg(bioreg_name = "eco_fresh", save_dir = NULL)
eco_marine = read_bioreg(bioreg_name = "eco_marine", save_dir = NULL)

## for each species with occurrence data 
i = 1
ranges = NULL
while(i <= length(files)) {
  ## read in filtered occurrence data
  occ = qs::qread(paste0("/Volumes/NIKKI/Bioshfits_GBIF/GBIF_data/", files[i]))
  sp_name = str_split_fixed(files, "\\.qs", 2)[i,1]
  
  realm = key$Eco[which(key$sp_name_checked == sp_name)]
  
  ## select correct ecoregion
  if("Ter" %in% realm) {
    ecoregion = eco_terra
    name = "ECO_NAME"
    r = "Ter"
  }
  else if("Mar" %in% realm) {
    ecoregion = eco_marine
    name = "REALM"
    r = "Mar"
  }
  else(
    ## for now, if other realm, don't make a file
    ecoregion = NA
  )
  
  if(!is.na(ecoregion)[1]) {
    ## draw the range 
    range = get_range(occ_coord = occ,
                      bioreg = ecoregion,
                      bioreg_name = name)
    
    ## plot the range + the occurrence records
    plot <- ggplot() +
      tidyterra::geom_spatraster(data = range, aes(fill = layer)) +
      scale_fill_continuous(na.value = "transparent") +
      geom_point(data = occ, aes(x = decimalLongitude, y = decimalLatitude), size = 0.1,
                 inherit.aes = FALSE) +
      labs(x = "Longitude", y = "Latitude", title = sp_name) +
      theme_bw() +
      theme(legend.position = "none")
    
    ## save plot
    ggsave(plot, path = "figures/gbif_range_polygons/", filename = paste0(sp_name, ".png"),
           width = 6, height = 4)
    
    range = st_as_sf(as.polygons(range))
    
    range <- range %>%
          mutate(species_name = unique(key$species_name[which(key$sp_name_checked == sp_name)]),
                 range_source = "GBIF occurrence",
                 realm = r) %>%
          select(-layer)
    
    st_write(range, paste0("data-processed/large-data/GBIF_polygons/", sp_name, ".shp"), append = FALSE)
    
    ## save the range 
    # if(is.null(ranges)) {
    #   ranges = range %>%
    #     mutate(species_name = unique(key$species_name[which(key$sp_name_checked == sp_name)]),
    #            range_source = "GBIF occurrence",
    #            realm = r) %>%
    #     select(-layer)
    # }
    # else {
    #   ranges <- range %>%
    #     mutate(species_name = unique(key$species_name[which(key$sp_name_checked == sp_name)]),
    #            range_source = "GBIF occurrence", 
    #            realm = r) %>%
    #     select(-layer) %>%
    #     rbind(ranges, .)
    # }
  }
  
  print(paste0("On species no. ", i, " out of ", length(unique(files))))
  i = i + 1
}

#388, 773, 793, 2000, 2184, 3048, 3828, 3836, 3845, 3988, 4141

##4474


## save ranges shapefile
#st_write(ranges, "data-processed/large-data/GBIF_polygons_v3.shp", append = FALSE)






