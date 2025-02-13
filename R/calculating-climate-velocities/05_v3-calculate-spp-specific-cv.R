## recalculate climate velocity measurements per species per study area
library(tidyverse)
library(terra)
library(sf)

## read in each study area climate velocity tif
files = list.files(path = "data-raw/bioshiftsv3/Velocity_SA", pattern="tif$", recursive=TRUE, full.names=FALSE)
files

## get rid of layers that are not the latitudinal component
files <- files[-which(str_detect(files, "_gVel.tif"))]
files <- files[-which(str_detect(files, "spatgrad"))]
files <- files[-which(str_detect(files, "avg_climate_layers"))]
files <- files[-which(str_detect(files, "trend"))]

study_ids <- unique(paste(str_split_fixed(files, "\\_", 3)[,1], str_split_fixed(files, "\\_", 3)[,2], sep = "_"))

## make sure we have velocities for all study areas
v3 = read.csv("data-raw/bioshiftsv3/BIOSHIFTS_v3.csv")
unique(v3$ID)[which(!unique(v3$ID) %in% study_ids)]
## some are missing
## A100, A171, A191 too small to calculate climate velocity

sf_use_s2(FALSE)

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

          ## get cv type for each layer
          cv_type = str_split_fixed((str_split_fixed(layers, paste0(study, "_"), 2)[,2]), ".tif", 2)[,1]
          
          ## for each cv layer, calculate cv metrics across entire study area
          l = 1
          while(l <= length(layers)) {
            r = rast(layers[l]) %>%
              project(., crs(study_polys))
            
            mean_cv_studylevel = mean(values(r), na.rm = TRUE)
            sd_cv_studylevel = sd(values(r), na.rm = TRUE)
            
            ## calculate species-specific cv metrics 
            ## make sure all extents overlap
            ol <- is.related(vect(study_polys), r, "intersects")

            if(any(ol == FALSE)) {
              ## remove the non-overlapping poly
              index = c(1:nrow(study_polys))[which(ol == TRUE)]
            }
            else{
              index = c(1:nrow(study_polys))
            }
            
            ## crop cv raster by the sp specific study polgyons
            cropped = lapply(index, function(x) {(crop(r, study_polys[x,], mask = TRUE))})
            
            ## calculate mean per sp specific polygon 
            cropped_means = unlist(lapply(1:length(cropped), function(x) {mean(values(cropped[[x]]), na.rm = TRUE)}))
            cropped_sds = unlist(lapply(1:length(cropped), function(x) {sd(values(cropped[[x]]), na.rm = TRUE)}))
            cropped_p_10 = unlist(lapply(1:length(cropped), function(x) {quantile(values(cropped[[x]]), c(0.10), na.rm = TRUE)}))
            cropped_p_90 = unlist(lapply(1:length(cropped), function(x) {quantile(values(cropped[[x]]), c(0.90), na.rm = TRUE)}))

            ## save
            if(is.null(cv_data)) {
              cv_data <- data.frame(study_id = rep(study, length(cropped_means)),
                                    cv_type = rep(cv_type[l], length(cropped_means)),
                                    mean_cv_studylevel = rep(mean_cv_studylevel, length(cropped_means)),
                                    sd_cv_studylevel = rep(sd_cv_studylevel, length(cropped_means)),
                                    species_studyid = study_polys$spcs_st[index],
                                    range_source = study_polys$rng_src[index],
                                    mean_cv_sppspecific = cropped_means,
                                    sd_cv_sppspecific = cropped_sds,
                                    p10_cv_sppspecific = cropped_p_10,
                                    p90_cv_sppspecific = cropped_p_90)
            } else {
              cv_data <- rbind(cv_data, 
                               data.frame(study_id = rep(study, length(cropped_means)),
                                          cv_type = rep(cv_type[l], length(cropped_means)),
                                          mean_cv_studylevel = rep(mean_cv_studylevel, length(cropped_means)),
                                          sd_cv_studylevel = rep(sd_cv_studylevel, length(cropped_means)),
                                          species_studyid = study_polys$spcs_st[index],
                                          range_source = study_polys$rng_src[index],
                                          mean_cv_sppspecific = cropped_means,
                                          sd_cv_sppspecific = cropped_sds,
                                          p10_cv_sppspecific = cropped_p_10,
                                          p90_cv_sppspecific = cropped_p_90))
            }
            
            l = l + 1
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

## check which ones have no overlap between species range x study area
unique(v3$ID)[which(!paste0("cv_", unique(v3$ID), ".csv") %in% list.files("data-processed/new-climate-velocities/"))]
