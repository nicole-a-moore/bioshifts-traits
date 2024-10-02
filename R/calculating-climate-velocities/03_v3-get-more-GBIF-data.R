## script to get GBIF data for species without range polygon from another source 
## BO script adapted by NM
library(tidyverse)
library(rgbif)
library(parallel)
library(data.table)

# set the settings
computer = "matrics"
work_dir <- getwd()
source("R/calculating-climate-velocities/settings.R")
source("R/calculating-climate-velocities/decompress_file.R")
range_env_data <- range(c(temporal_range_env_data("Ter"),temporal_range_env_data("Mar")))


####################################
#    GET GBIF CODES FOR SPECIES    #
####################################
## get list of species missing range polygons 
missing <- read.csv("data-processed/v3_spp-missing-range-polygons.csv")

## get their gbif codes
codes <- read.csv("data-raw/bioshiftsv3/splist_v3.csv") %>%
  select(species, kingdom, phylum, class, order, family, db, db_code) %>%
  unique()

missing <- left_join(missing, codes, by = c("species_name"= "species", "class" = "class")) %>%
  filter(!is.na(db_code))

missing$db_code <- gsub("GBIF:","", missing$db_code)

missing <- filter(missing, !str_detect(missing$db_code, "\\:"))


##############################
#    SUBMIT GBIF REQUESTS    #
##############################
## request occurrence downloads for each species 
## not: 200,1,3,4,206 throw error 404
requests <- data.frame()
inforequest <- list()
for(i in 1264:nrow(missing)){ 
    
    print(paste0("Requesting species ", i))
    
    inforequest[[i]] <- occ_download(
          #pred("taxonKey", 8417931),
          pred_in("speciesKey", missing$db_code[i]),
          pred_in("basisOfRecord",
                  basisOfRecord),
          pred_gte("year", range_env_data[1]),
          pred_lte("year", range_env_data[2]),
          pred("hasCoordinate", TRUE),
          pred("hasGeospatialIssue", FALSE),
          pred("occurrenceStatus", "PRESENT"),
          format = "SIMPLE_CSV",
          user = Sys.getenv('GBIF_USER'), 
          pwd = Sys.getenv('GBIF_PWD'),
          email = "nicole.moore@mail.mcgill.ca")
        
        tmp <- inforequest[[i]]
        
        # save information about request
        tmp.save <- attributes(tmp)
        tmp.save <- data.frame(download_key = tmp[1],
                               created = tmp.save$created,
                               download_link = tmp.save$downloadLink,
                               doi = tmp.save$doi,
                               citation = tmp.save$citation,
                               format = tmp.save$format,
                               user = tmp.save$user,
                               email = tmp.save$email)
        requests <- rbind(requests,tmp.save)
    
        occ_download_wait(inforequest[[i]])
}

## check  
library(purr)
df <-
  inforequest %>%
  map_df(
    ~ data.frame(
      download_key = .[1],
      created = attr(., "created"),
      citation = attr(., "citation"),
      format = attr(., "format"),
      user = attr(., "user"),
      email = attr(., "email"),
      doi    =  attr(., "doi"),
      download_link    =  attr(., "downloadLink")
    )
  )

requests <- requests[-which(!requests$doi %in% df$doi),]

requests <- rbind(requests, df[which(!df$doi %in% requests$doi),])

write.csv(requests, "data-processed/gbif_requests.csv", row.names = F)

############################################
#    DOWNLOAD OCCURRENCES FROM REQUESTS    #
############################################
# Download requests at the HPC ----
requests <- read.csv("data-processed/gbif_requests.csv")
head(requests)

GBIF_zip_dir <- "/Volumes/NIKKI/Bioshfits_GBIF/GBIF_data_NM"
if(!dir.exists(GBIF_zip_dir)){
    dir.create(GBIF_zip_dir,recursive = T)
}
ncores = parallelly::availableCores()

check <- "Error"
attempt <- 0
errors <- TRUE

keystogo <- requests$download_key
while( any(errors) ) {
    attempt <- attempt + 1
    
    cat("Attempt", attempt, "\n")
    
    check <- mclapply(1:length(keystogo), function(i){
        
        keystogo_i = as.character(keystogo[i])
        
        test = occ_download_meta(keystogo_i)
        test = test$status
        
        while(test == "RUNNING"){
            
            test = occ_download_meta(keystogo_i)
            test = test$status
            
            Sys.sleep(60)
            
        }
        
        check[[i]] <- try (
            {
                rgbif::occ_download_get(key = keystogo_i,
                                        path = GBIF_zip_dir,
                                        overwrite = TRUE)
                return(paste("OK", i))
            },
            silent = TRUE
        )
        
    }
    , mc.cores = ifelse(ncores>length(keystogo),length(keystogo),ncores))
    
    errors <- sapply(check, function(x) class(x)=="try-error")
    
    if(any(errors)){
        errors <- requests$download_key[which(errors)]
        keystogo <- errors
    }
    cat("Error on keys:", errors)
} 


#################################################
#    TEST IF EVERYTHING DOWNLOADED CORRECTLY    #
#################################################
# Test if everything downloaded correctly
# Compare file sizes of remote and local files
download_size <- function(url) as.numeric(httr::HEAD(url)$headers$`content-length`)

for(i in 1:length(keystogo)){ cat("\rChecking file", i, "from", length(keystogo))
    # local file size
    size_local <- file.size(here::here(GBIF_zip_dir, paste0(keystogo[i],".zip")))
    # remote file size
    size_remote <- download_size(requests$download_link[i])
    # download again if files dont have the same size
    if(!size_local == size_remote){
        cat("\rFiles don't have the same size\nDownloading file", i)
        rgbif::occ_download_get(key = keystogo[i],
                                path = GBIF_zip_dir,
                                overwrite = TRUE)
    }
}

#####################################
#    FILTER AND SAVE OCCURRENCES    #
#####################################
## add ecosystem 
v3 = read.csv("data-raw/bioshiftsv3/BIOSHIFTS_v3.csv")
v3$species_name = str_replace_all(v3$sp_name_checked, "\\_", " ")
v3$species_name

v3 <- v3 %>%
  mutate(Eco = ifelse(Eco %in% c("Mar"), "Mar", "Ter")) %>%
  select(species_name, Eco) %>%
  distinct()

missing <- left_join(missing, v3)

GBIF_zip_dir <- "/Volumes/NIKKI/Bioshfits_GBIF/GBIF_data_NM"
if(!dir.exists(GBIF_zip_dir)){
    dir.create(GBIF_zip_dir,recursive = T)
}

myselection <- c("basisOfRecord", "speciesKey", "decimalLongitude", "decimalLatitude", "year", "month","species")
basisOfRecord_selected <- basisOfRecord

terrestrials <- missing$species_name[which(missing$Eco == "Ter")]
marines <- missing$species_name[which(missing$Eco == "Mar")]

# my terrestrial raster
ter.ras <- terra::rast("data-raw/model_raster_ter_1km.tif")
# my marine raster
mar.ras <- terra::rast("data-raw/model_raster_mar.tif")

# test
any(terrestrials %in% marines)

# zipped occ files
occs <- list.files(here::here(GBIF_zip_dir), pattern = ".zip")

# create dir to save sps occurrences
occ_dir = "/Volumes/NIKKI/Bioshfits_GBIF/GBIF_occurrences_NM"
if(!dir.exists(occ_dir)){
    dir.create(occ_dir, recursive = T)
}

# create temp dir to decompress zip file
tmp.dir <- here::here(occ_dir,"tmp")
if(!dir.exists(tmp.dir)){
    dir.create(tmp.dir)
}


ncores = parallelly::availableCores()

test_if_work <- mclapply(1:length(occs), function(i){
    
    
    # decompress zip file
    zippedfile <- here::here(GBIF_zip_dir, occs[i])
    decompress_file(directory = tmp.dir, file = zippedfile)
    unzippedfile <- gsub(".zip",".csv",occs[i])
    unzippedfile <- gsub(".zip",".csv",here::here(tmp.dir,unzippedfile))
    
    # read in
    tmp <- fread(unzippedfile, select = myselection, nThread = 10)
    tmp <- na.omit(tmp) 
    
    # filter basisOfRecord
    tmp <- tmp %>% dplyr::filter(basisOfRecord %in% basisOfRecord_selected)
    tmp_species <- unique(tmp$species)
    
    # remove species I already have data
    # saved species
    saved_sps <- list.files(occ_dir,pattern = ".qs")
    saved_sps <- gsub(".qs","",saved_sps)
    saved_sps <- gsub("_"," ",saved_sps)
    
    if(any(tmp_species %in% saved_sps)){
        tmp <- tmp %>% dplyr::filter(!species %in% saved_sps)
        tmp_species <- unique(tmp$species)
    }
    
    if(length(tmp_species)>0){
        
        tmp <- tmp %>% dplyr::filter(species %in% tmp_species)
        
        # subset terrestrials
        ter. <- tmp[which(tmp$species %in% terrestrials),]
        # get cells
        cells. <-  terra::extract(ter.ras,
                                  ter.[,c("decimalLongitude", "decimalLatitude")],
                                  cells=TRUE)
        ter.$cell <- cells.$cell
        ter.$layer <- cells.[,2]
        # remove cells falling in the ocean
        ter. <- ter.[which(ter.$layer==1),] # keep land / remove ocean
        ter. = ter.[,-"layer"]
        # remove occurrences outside the temporal range of the env data
        ter. <- ter. %>% filter(year >= temporal_range_env_data("Ter")[1]+1 + n_yr_bioclimatic # because we cannot calculate bioclimatics for the year 1
                                & year <= temporal_range_env_data("Ter")[2])
        
        # subset marines
        mar. <- tmp[which(tmp$species %in% marines),] 
        # get cells
        cells. <-  terra::extract(mar.ras, 
                                  mar.[,c("decimalLongitude", "decimalLatitude")],
                                  cells=TRUE)
        mar.$cell <- cells.$cell
        mar.$layer <- cells.[,2]
        # remove cells falling in land
        mar. <- mar.[which(mar.$layer==1),] # remove land / keep ocean
        mar. = mar.[,-"layer"]
        # remove occurrences outside the temporal range of the env data
        mar. <- mar. %>% filter(year > temporal_range_env_data("Mar")[1] + n_yr_bioclimatic # because we cannot calculate bioclimatics for the year 1
                                & year <= temporal_range_env_data("Mar")[2])
        
        # Group all
        sps. <- rbind(ter.,mar.)
        # remove temporal ter. and mar. data to save memory space
        rm(ter., mar.);gc()
        
        # Remove duplicates: same species in the same location and date
        rm <- duplicated(sps.[,c("species", "cell", "year", "month")])
        if(any(rm)){
            sps. <- sps.[-which(rm),]
        }
        
        # Remove other potential issues with the package CoordinateCleaner
        sps. = CoordinateCleaner::clean_coordinates(
            x = sps.,
            lon = "decimalLongitude",
            lat = "decimalLatitude",
            tests = c("capitals", "centroids", "equal",
                      "gbif", "institutions", "zeros"),
            value = "clean")
        
        # keep species with more than 30 obs
        my_sps_i <- table(sps.$species)
        my_sps_i <- names(my_sps_i[which(my_sps_i >= 30)])
        sps. <- sps.[which(sps.$species %in% my_sps_i),]
        
        # save data per species
        if(length(my_sps_i)==1){
            for(j in 1:length(my_sps_i)){ 
                cat("\rsaving sps", j, "from", length(my_sps_i))
                tmp_sps <- subset(sps., species == my_sps_i[j])
                
                qs::qsave(tmp_sps, 
                          here::here(occ_dir, paste0(gsub(" ","_",my_sps_i[j]),".qs")))
            }
        }
        
        # delete unzipped file
        unlink(unzippedfile)
        
    }
    
    
}
, mc.cores = ncores)

# delete tmp dir used to decompress zipfiles
unlink(tmp.dir, recursive = TRUE)
unlink(GBIF_zip_dir, recursive = TRUE)


################################
#    MAKE MORE CONVEX HULLS    #
################################
library(gbif.range)
library(sf)
library(terra)
occ_dir = "/Volumes/NIKKI/Bioshfits_GBIF/GBIF_occurrences_NM"

## read in bioshifts 
v3 = read.csv("data-raw/bioshiftsv3/BIOSHIFTS_v3.csv")
v3$species_name = str_replace_all(v3$sp_name_checked, "\\_", " ")
v3$species_name

## read GBIF occurrence filenames
files = list.files(occ_dir)

## get realm from v3
key = v3 %>%
  filter(sp_name_checked %in% str_split_fixed(files, "\\.qs", 2)[,1]) %>%
  mutate(Eco = ifelse(Eco %in% c("Mar"), "Mar", "Ter")) %>%
  select(sp_name_checked, Eco, species_name) %>%
  distinct()

## get the ecoregions 
eco_terra = read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)
eco_fresh = read_bioreg(bioreg_name = "eco_fresh", save_dir = NULL)
eco_marine = read_bioreg(bioreg_name = "eco_marine", save_dir = NULL)

## for each species with occurrence data 
# 721
i = 1
while(i <= length(files)) {
  ## read in filtered occurrence data
  occ = qs::qread(paste0(occ_dir, "/", files[i]))
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
    
    ## save the range
    st_write(range, paste0("data-processed/large-data/GBIF_polygons_new/", sp_name, ".shp"), append = FALSE)
  }
  
  print(paste0("On species no. ", i, " out of ", length(unique(files))))
  i = i + 1
}


## read all together and save 
files = list.files(path = "data-processed/large-data/GBIF_polygons_new", 
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
st_write(gbif, "data-processed/large-data/GBIF_polygons_new/GBIF_v3_2.shp", append = FALSE)

gbif2 <- st_read("data-processed/large-data/GBIF_polygons_new/GBIF_v3_2.shp")
gbif2 <- gbif2[, which(colnames(gbif2) %in% c("spcs_nm", "geometry"))]
gbif2 <- rename(gbif2, "binomial" = spcs_nm)
gbif2$range_source = "GBIF occurrence"

## add to ranges RDS
ranges <- readRDS("data-processed/large-data/collated-ranges.rds")
nrow(ranges)

ranges <- rbind(ranges, gbif2)
length(unique(ranges$binomial)) ## 11666

saveRDS(ranges, "data-processed/large-data/collated-ranges_ALL.rds")
