## inspect species-specific climate velocities
library(tidyverse)
library(sf)
library(ggpattern)
library(stringr)

## read in uncut study polygons 
polys <- st_read("data-processed/bioshiftsv3_studypolygons.shp")

## read in ranges 
ranges <- readRDS("data-processed/large-data/collated-ranges.rds")
ranges$id = paste(ranges$binomial, ranges$range_source, sep = "_")

## for each study area:
i = 1
while(i <= nrow(polys)) {
  
  study_id <- polys$Name[i]
  full_poly <- polys[which(polys$Name == study_id),]
  
  ## create folder to store figures for this study
  ifelse(!dir.exists(paste0("figures/new-climate-velocities-shapes/", study_id)), 
         dir.create(paste0("figures/new-climate-velocities-shapes/", study_id)), FALSE)
  
  ## if there is a cv file
  if(file.exists(paste0("data-processed/new-climate-velocities/cv_", study_id, ".csv"))) {
    ## read in cv file
    cvs <- read.csv(paste0("data-processed/new-climate-velocities/cv_", study_id, ".csv"))
    cvs$species_name <- str_split_fixed(cvs$species_studyid, "\\_", 2)[,1]
    
    ## read in cut study polygons 
    cut_polys <- st_read(paste0("data-processed/spp-specific-study-polygons/", study_id, ".shp"))
    cut_polys <- unique(cut_polys)
    
    ## make unique id
    cvs$cv_id = paste(cvs$species_studyid, cvs$range_source, sep = "_")
    cut_polys$cv_id = paste(cut_polys$spcs_st, cut_polys$rng_src, sep = "_")
    
    ## for each cut poly in study
    p = 1
    while(p <= length(unique(cut_polys$cv_id))) {
      ## match to full study poly, cut study polys, and range and plot
      cut_poly <- cut_polys[which(cut_polys$cv_id == cvs$cv_id[p]),]
      range <- ranges[which(ranges$id == paste(cut_poly$binomil, cut_poly$rng_src, sep = "_")),]
      
      ## get cv
      cv <- cvs[which(cvs$cv_id == cut_poly$cv_id),] %>%
        select(cv_type, mean_cv_studylevel, mean_cv_sppspecific) %>%
        unique() 
      
      label = paste(paste(cv$cv_type, rep(" ", nrow(cv)), 
                          rep("CV whole study area: ", nrow(cv)), cv$mean_cv_studylevel,
                          collapse = "\n"),
                    paste(cv$cv_type, rep(" ", nrow(cv)), 
                          rep("CV species-specific: ", nrow(cv)), cv$mean_cv_sppspecific,
                          collapse = "\n"),
                    sep = "\n")
      
      ## plot with cv values on plot
      plot <- ggplot(cut_poly) +
        geom_sf(fill = "red", alpha = 0.5, aes(colour = "Cropped study area")) +
        geom_sf(data = range, fill = "blue", alpha = 0.1, aes(colour = "Range polygon")) +
        geom_sf_pattern(data = full_poly, fill = "transparent", 
                        pattern_density = 0.01, 
                        aes(colour = "Original study area")) +
        labs(title = cut_poly$cv_id, subtitle = label,
             colour = NULL) +
        theme_bw() +
        theme(panel.grid = element_blank(), legend.position = "bottom") +
        scale_colour_manual(values = c("transparent", "black", "transparent"))
      
      ggsave(plot, path = paste0("figures/new-climate-velocities-shapes/", study_id), 
             filename = paste0(cut_poly$cv_id, ".png"), width = 8, height = 8)
      
      print(paste0("On species ", p, " of study no. ", i))
      
      p = p + 1
    }
  
  }
  
  i = i + 1
}
