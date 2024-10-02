## recalculate climate velocity per species per study area
library(tidyverse)
library(terra)
library(sf)

## read in each study area climate velocity tif
files = list.files(path = "data-raw/bioshiftsv3/Velocity_SA", pattern="tif$", recursive=TRUE, full.names=FALSE)
files

study_ids <- unique(paste(str_split_fixed(files, "\\_", 3)[,1], str_split_fixed(files, "\\_", 3)[,2], sep = "_"))

## make sure we have velocities for all study areas
v3 = read.csv("data-raw/bioshiftsv3/BIOSHIFTS_v3.csv")
unique(v3$ID)[which(!unique(v3$ID) %in% study_ids)]

cv_data = NULL

i = 1
while(i <= length(study_ids)) {
  study = study_ids[i]
  
  ## if not already computed
  if(!file.exists(paste0("data-processed/new-climate-velocities/cv_", study, ".csv")))  {
    ## read in cropped study polygons for this study
    if(file.exists(paste0("data-processed/spp-specific-study-polygons/", study, ".shp"))) {
      study_polys <- st_read(paste0("data-processed/spp-specific-study-polygons/", study, ".shp"))
      
      ## get all cv layers for this study area
      layers <- paste0("data-raw/bioshiftsv3/Velocity_SA/", files[which(str_detect(files, study))])
      
      if(length(layers)!=0) {
        layers <- layers[which(!str_detect(layers, "gVel.tif"))] # get rid of non-northward layers
        
        if(length(layers)!=0) {
          if(str_detect(layers[1], "mat")) {
            layers <- layers[which(str_detect(layers, "gVelLat") | str_detect(layers, "gVelEle"))] # and get rid of velocity layer that isn't northward
          }
          
          cv_type = str_split_fixed((str_split_fixed(layers, paste0(study, "_"), 2)[,2]), ".tif", 2)[,1]
          
          ## calculate cv metrics across entire study area
          l = 1
          while(l <= length(layers)) {
            r = rast(layers[l]) %>%
              project(., crs(study_polys))
            
            mean_cv_studylevel = mean(values(r), na.rm = TRUE)
            sd_cv_studylevel = sd(values(r), na.rm = TRUE)
            
            ## calculate species-specific cv metrics 
            
            ## crop cv raster by the sp specific study polgyons
            cropped = lapply(1:nrow(study_polys), function(x) {(crop(r, study_polys[x,], mask = TRUE))})
            
            ## calculate mean per sp specific polygon 
            cropped_means = unlist(lapply(1:length(cropped), function(x) {mean(values(cropped[[x]]), na.rm = TRUE)}))
            cropped_sds = unlist(lapply(1:length(cropped), function(x) {sd(values(cropped[[x]]), na.rm = TRUE)}))
            
            ## save
            if(is.null(cv_data)) {
              cv_data <- data.frame(study_id = rep(study, length(cropped_means)),
                                    cv_type = rep(cv_type[l], length(cropped_means)),
                                    mean_cv_studylevel = rep(mean_cv_studylevel, length(cropped_means)),
                                    sd_cv_studylevel = rep(sd_cv_studylevel, length(cropped_means)),
                                    species_studyid = study_polys$spcs_st,
                                    range_source = study_polys$rng_src,
                                    mean_cv_sppspecific = cropped_means,
                                    sd_cv_sppspecific = cropped_sds)
            } else {
              cv_data <- rbind(cv_data, 
                               data.frame(study_id = rep(study, length(cropped_means)),
                                          cv_type = rep(cv_type[l], length(cropped_means)),
                                          mean_cv_studylevel = rep(mean_cv_studylevel, length(cropped_means)),
                                          sd_cv_studylevel = rep(sd_cv_studylevel, length(cropped_means)),
                                          species_studyid = study_polys$spcs_st,
                                          range_source = study_polys$rng_src,
                                          mean_cv_sppspecific = cropped_means,
                                          sd_cv_sppspecific = cropped_sds))
            }
            
            l = l + 1
          }
        }
        
      }
    }
    
    ## save for study 
    ## if not empty
    if(!is.null(cv_data)) {
      write.csv(cv_data, paste0("data-processed/new-climate-velocities/cv_", study, ".csv"), row.names = FALSE)
    }
    
    ## reset dataframe 
    cv_data = NULL
  }
  
  i = i + 1

  print(paste0("On study area ", i, " of ", length(study_ids)))
}

cv_data %>%
  filter(cv_type == "mat_gVelLat") %>%
  ggplot(aes(x = mean_cv_sppspecific, fill = study_id)) +
  geom_histogram() +
  facet_wrap(cv_type~study_id) +
  geom_vline(aes(xintercept = mean_cv_studylevel)) +
  labs(x = "Mean climate velocity", y = "No. species")

cv_data %>%
  filter(cv_type == "mat_gVelEle") %>%
  ggplot(aes(x = mean_cv_sppspecific, fill = study_id)) +
  geom_histogram() +
  facet_wrap(cv_type~study_id) +
  geom_vline(aes(xintercept = mean_cv_studylevel)) +
  labs(x = "Mean climate velocity", y = "No. species")


## try plotting against range shift 
v3_sub <- cv_data %>%
  mutate(species_name = str_split_fixed(.$species_studyid, "\\_", 2)[,1],
         ID =  str_split_fixed(.$species_studyid, "\\_", 2)[,2]) %>%
  left_join(., v3)

v3_sub %>%
  filter(cv_type == "mat_gVelLat") %>%
  ggplot(aes(x = mean_cv_studylevel, y = Rate)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_x_continuous(limits = c(0,2)) +
  labs(x = "Climate velocity", y = 'Range shift rate')

v3_sub %>%
  filter(cv_type == "mat_gVelLat") %>%
  ggplot(aes(x = mean_cv_sppspecific, y = Rate)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_x_continuous(limits = c(0,2))  +
  labs(x = "Climate velocity", y = 'Range shift rate')


v3_sub %>%
  filter(cv_type == "mat_gVelEle") %>%
  ggplot(aes(x = mean_cv_studylevel, y = Rate)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()  +
  labs(x = "Climate velocity", y = 'Range shift rate')


v3_sub %>%
  filter(cv_type == "mat_gVelEle") %>%
  ggplot(aes(x = mean_cv_sppspecific, y = Rate)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  labs(x = "Climate velocity", y = 'Range shift rate')


## PLOTS TO MAKE FOR MEETING
## - realized range source histogram
## - realized range bar plot by Class
## - example of study area A10_P1 (North America) + cropped study areas 
## - histrogram of spp-specific cv, colour = study, with vertical line showing study mean 
## - histrogram of sd of spp-specific climate velocity with vertical line showing study sd (expect it to be high relative to spp-specific)

## if time: range shift versus cv (study mean versus species-specific)


test <- read.csv("data-raw/bioshiftsv3/Velocity_SA/A1_P1.csv") 

## read one
test <- rast("data-raw/bioshiftsv3/Velocity_SA/A10_P1_mat_gVel.tif") 
testlat <- rast("data-raw/bioshiftsv3/Velocity_SA/A10_P1_mat_gVelLat.tif") ## this one is velocity north 
test

plot(test)
plot(testlat)


## get all climate velocity layers fo

## calculate climate velocity metrics for whole study area
## check against Brunno's calculation

## crop cv tif by each species x study area combo