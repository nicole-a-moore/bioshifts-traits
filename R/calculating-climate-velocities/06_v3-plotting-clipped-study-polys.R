## inspect species-specific climate velocities
library(tidyverse)
library(sf)
library(ggpattern)
library(stringr)
select = dplyr::select

## read in uncut study polygons 
polys <- st_read("data-processed/bioshiftsv3_studypolygons.shp")

## read in ranges 
ranges <- readRDS("data-processed/large-data/collated-ranges_ALL.rds")
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
      
      # label = paste(paste(cv$cv_type, rep(" ", nrow(cv)), 
      #                     rep("CV whole study area: ", nrow(cv)), cv$mean_cv_studylevel,
      #                     collapse = "\n"),
      #               paste(cv$cv_type, rep(" ", nrow(cv)), 
      #                     rep("CV species-specific: ", nrow(cv)), cv$mean_cv_sppspecific,
      #                     collapse = "\n"),
      #               sep = "\n")
       
      ## plot with cv values on plot
      plot <- ggplot(cut_poly) +
        geom_sf(fill = "red", alpha = 0.5, aes(colour = "Cropped study area")) +
        geom_sf(data = range, fill = "blue", alpha = 0.1, aes(colour = "Range polygon")) +
        geom_sf_pattern(data = full_poly, fill = "transparent", 
                        pattern_density = 0.01, 
                        aes(colour = "Original study area")) +
        labs(title = cut_poly$cv_id, #subtitle = label,
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












## garbage

## read in cut study polygons 
files <-list.files("data-processed/spp-specific-study-polygons/",  full.names=TRUE, pattern="shp$")
cut_polys <- c()
for(i in files) {
  cut_polys <- rbind(cut_polys, st_read(i))
}

length(unique(cut_polys$Name))


# ## calculate metrics about the shapefiles
# ## area, max and min latitude and longitude 
# sf_use_s2(FALSE)
# 
# area_study_pols <- polys %>%
#   select(Name) %>%
#   mutate(area_m2 = st_area(.)) 
# 
# area_cut_pols <- cut_polys %>%
#   select(Name, spcs_st, rng_src) %>%
#   mutate(area_m2 = st_area(.))
# 
# area_ranges <- ranges %>%
#   select(binomial, range_source) %>%
#   mutate(area_m2 = st_area(.))
# 
# 
# ## make sf dataframe into list
# area_study_pols_split <- area_study_pols %>%
#   split(f = .$Name)
# 
# # find bboxes
# bbox_list <- map(
#   .x = area_study_pols_split,
#   .f = ~st_bbox(.x) %>% unlist() %>% c()
# )
# 
# # unsplit
# df <- bbox_list %>% bind_rows(.id = "Name")
# 
# ## get max and min latitude of species range within min and max longitude of study area 
# ## make boxes
# # split by id
# df_spit <- df %>%
#   split(f = .$Name)
# 
# # create bbox polygons from xmin to xmax, -90 to 0 and 0 to 90
# boxes_S <- purrr::map(
#   .x = df_spit,
#   .f = ~st_bbox(c(xmin = .x$xmin,
#                   xmax = .x$xmax,
#                   ymin = -90,
#                   ymax = 0)) %>% 
#     st_as_sfc(crs = st_crs(area_study_pols)) %>%
#     st_set_crs(st_crs(area_study_pols)) %>%
#     st_as_sf()
# ) %>%
#   bind_rows(.id = "Name")
# 
# boxes_N <- purrr::map(
#   .x = df_spit,
#   .f = ~st_bbox(c(xmin = .x$xmin,
#                   xmax = .x$xmax,
#                   ymin = 0,
#                   ymax = 90)) %>% 
#     st_as_sfc(crs = st_crs(area_study_pols)) %>%
#     st_set_crs(st_crs(area_study_pols)) %>%
#     st_as_sf()
# ) %>%
#   bind_rows(.id = "Name")
# 
# # split ranges by study id
# # to make stacks of all the species ranges within one study
# ## add columns 
# test <- data.frame(st_drop_geometry(area_cut_pols)) %>%
#   mutate(binomial = str_split_fixed(spcs_st, "_", 2)[,1]) %>%
#   select(Name, binomial, rng_src) %>%
#   rename("range_source" = rng_src)
# 
# area_ranges_new <- left_join(area_ranges, test)
# 
# ranges_split <- area_ranges_new %>% split(f = .$Name)
# 
# # split study areas bounding box polygons by study id - should be just one each
# boxes_N <- boxes_N %>% split(f = .$Name)
# boxes_S <- boxes_S %>% split(f = .$Name)
# 
# # make vector of all study area ids
# # because this is how we'll call elements of each list
# study_area_names_N <- names(boxes_N)
# study_area_names_S <- names(boxes_S)
# 
# # crop full stack of ranges in each study id cluster 
# # by the bounding box for that study
# boxes_N <- boxes_N[which(names(boxes_N) %in% names(ranges_split))]
# boxes_S <- boxes_S[which(names(boxes_S) %in% names(ranges_split))]
# 
# study_area_names_N <- study_area_names_N[which(study_area_names_N %in% names(ranges_split))]
# study_area_names_S <- study_area_names_S[which(study_area_names_S %in% names(ranges_split))]
# 
# cropped_ranges_N <- purrr::map(
#   # map over ids of all study areas
#   .x = study_area_names_N,
#   .f = ~st_crop(ranges_split[[.x]],
#                 boxes_N[[.x]]),
#   .progress = TRUE) %>%
#   set_names(study_area_names_N)
# 
# cropped_ranges_S <- purrr::map(
#   # map over ids of all study areas
#   .x = study_area_names_S,
#   .f = ~st_crop(ranges_split[[.x]],
#                 boxes_S[[.x]]),
#   .progress = TRUE) %>%
#   set_names(study_area_names_S)
# 
# cropped_ranges_bound_N <- cropped_ranges_N %>%
#   bind_rows(.id = "Name")
# cropped_ranges_bound_S <- cropped_ranges_S %>%
#   bind_rows(.id = "Name")
# 
# cropped_ranges_bound_N$hemisphere = "N"
# cropped_ranges_bound_S$hemisphere = "S"
# 
# cropped_ranges_bound <- rbind(cropped_ranges_bound_N, cropped_ranges_bound_S)
# saveRDS(cropped_ranges_bound, "data-processed/ranges_cropped_by_study_bbox.rds")
# 
# ## write out study polys 
# saveRDS(area_study_pols, "data-processed/study_polygons.rds")
# 
# which(cropped_ranges_bound)
# 
# 
# ## check Sylvia communis_A133_P2_GBIF occurrence
# A133_P2_N <- filter(cropped_ranges_bound_N, Name == "A133_P2")
# A133_P2_S <- filter(cropped_ranges_bound_S, Name == "A133_P2")
# 
# A133_P2_N$species_studyid <- paste0(A133_P2_N$binomial, "_A133_P2_", A133_P2_N$range_source)
# A133_P2_S$species_studyid <- paste0(A133_P2_S$binomial, "_A133_P2_", A133_P2_S$range_source)
#   
# 
# A133_P2_N %>%
#   filter(species_studyid == "Sylvia communis_A133_P2_GBIF occurrence") %>%
#   ggplot(aes()) +
#   geom_sf(fill = "red") +
#   geom_sf(data = filter(A133_P2_S, species_studyid == "Sylvia communis_A133_P2_GBIF occurrence"), 
#           fill = "blue") +
#   geom_sf(data = filter(area_study_pols, Name == "A133_P2"), 
#           fill = "green") 
#   
# 
# ## now crop the IUCN and birds of the world seasonal ranges by study area bbox 
# 
# 
# v3 %>%
#   filter(ID == "A133_P2") %>%
#   select(Type, Param, everything())
# 
# 
# 
# 
# area_study_pols = left_join(area_study_pols, df)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## make sf dataframe into list
# area_cut_pols_split <- area_cut_pols %>%
#   split(f = .$spcs_st)
# 
# # find bboxes
# bbox_list <- map(
#   .x = area_cut_pols_split,
#   .f = ~st_bbox(.x) %>% unlist() %>% c()
# )
# 
# # unsplit
# df <- bbox_list %>% bind_rows(.id = "spcs_st")
# 
# area_cut_pols = left_join(area_cut_pols, df)
# 
# 
# ## make sf dataframe into list
# area_ranges$binomial_rangesource = paste0(area_ranges$binomial, "_", area_ranges$range_source)
# area_ranges_split <- area_ranges %>%
#   split(f = .$binomial_rangesource)
# 
# # find bboxes
# bbox_list <- map(
#   .x = area_ranges_split,
#   .f = ~st_bbox(.x) %>% unlist() %>% c()
# )
# 
# # unsplit
# df <- bbox_list %>% bind_rows(.id = "binomial_rangesource")
# 
# area_ranges = left_join(area_ranges, df)
# 
# 
# colnames(area_study_pols)
# colnames(area_cut_pols)
# colnames(area_ranges)
# 
# ## drop geo
# area_study_pols <- data.frame(st_drop_geometry(area_study_pols))
# area_cut_pols <- data.frame(st_drop_geometry(area_cut_pols))
# area_ranges <- data.frame(st_drop_geometry(area_ranges))
# 
# 
# ## combine into list 
# list <- list(area_study_pols, area_cut_pols, area_ranges)
# names(list) <- c("study_polys", "cut_study_polys", "range_polygons")
# 
# list[[1]]
# 
# saveRDS(list, "data-processed/polygon_metrics.RDS")
# 
# 
# ## make data frame the has one column for area of range, one for area of range within study poly, one for area of study poly
# head(area_ranges)
# 
# 
# sp <- area_study_pols %>%
#   rename("ID" = Name, "Area_m2_ID" = area_m2) %>%
#   select(ID, Area_m2_ID)
# 
# r <- area_ranges %>%
#   rename("Area_m2_range" = area_m2) %>%
#   select(binomial, range_source, Area_m2_range)
# 
# cp <- area_cut_pols %>%
#   rename("ID" = Name, "species_studyid" = spcs_st, "range_source" = rng_src, 
#          "Area_m2_overlap" = area_m2) %>%
#   select(ID, Area_m2_overlap, species_studyid, range_source) %>%
#   mutate(binomial = str_split_fixed(species_studyid, "\\_", 2)[,1])
# 
# 
# ## join
# df <- left_join(cp, r) %>%
#   left_join(., sp) %>%
#   select(ID, binomial, species_studyid, range_source, everything())
# 
# ## write 
# write.csv(df, "data-processed/polygon_areas.csv", row.names = FALSE)
# 
# df <- read.csv("data-processed/polygon_areas.csv")
# 
# ## join to v3
# v3 <- read.csv("data-processed/v3_new-cvs_preliminary.csv")
# 
# v3 <- left_join(v3, df)
# 
# ## get rid of precip 
# v3 <- filter(v3, !str_detect(cv_type, "map"))
# 
# ## convert to km2
# v3 <- v3 %>%
#   mutate(Area_km2_ID = Area_m2_ID/1000000,
#          Area_km2_range = Area_m2_range/1000000,
#          Area_km2_overlap = Area_m2_overlap/1000000) 
# 
# ## plot:
# ## study level cv vs spp specific cv
# v3 %>%
#   ggplot(aes(x = mean_cv_studylevel, y = mean_cv_sppspecific)) +
#   geom_point(size = 0.1) +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~Eco)
# 
# ## study level cv vs shift 
# v3 %>%
#   ggplot(aes(x = mean_cv_studylevel, y = Rate)) +
#   geom_point(size = 0.1) +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~Eco) +
#   geom_smooth(method = "lm")
# 
# ## spp specific cv vs shift 
# v3 %>%
#   ggplot(aes(x = mean_cv_sppspecific, y = Rate)) +
#   geom_point(size = 0.1) +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~Eco) +
#   geom_smooth(method = "lm")
# 
# ## by edge  
# v3 %>%
#   mutate(Eco = ifelse(Eco %in% c("Mar", "Int"), "Mar", "Ter")) %>%
#   ggplot(aes(x = mean_cv_studylevel, y = Rate, colour = Eco)) +
#   geom_point(size = 0.1) +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~Param) +
#   geom_smooth(method = "lm")
# 
# v3 %>%
#   mutate(Eco = ifelse(Eco %in% c("Mar", "Int"), "Mar", "Ter")) %>%
#   ggplot(aes(x = mean_cv_sppspecific, y = Rate, colour = Eco)) +
#   geom_point(size = 0.1) +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~Param) +
#   geom_smooth(method = "lm")
# 
# ## colour by area of study overlap with range / area of range 
# ## look at centroid shifts - need to have high amount of range covered to be meaningful
# v3 %>%
#   filter(Param == "O") %>%
#   mutate(Eco = ifelse(Eco %in% c("Mar", "Int"), "Mar", "Ter"),
#          ratio_overlap_range = Area_km2_overlap/Area_km2_range) %>%
#   ggplot(aes(x = mean_cv_sppspecific, y = Rate, colour = Eco)) +
#   geom_point(aes(size = ratio_overlap_range)) +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~Eco) +
#   geom_smooth(method = "lm") +
#   scale_size(range = c(0, 2))
# 
# ## 1 = all of range is within study area
# 
# v3 %>%
#   filter(Param == "O") %>%
#   mutate(Eco = ifelse(Eco %in% c("Mar", "Int"), "Mar", "Ter"),
#          ratio_overlap_range = Area_km2_overlap/Area_km2_range) %>%
#   mutate(ratio_gt_75 = ratio_overlap_range >= 0.75) %>%
#   ggplot(aes(x = mean_cv_sppspecific, y = Rate, colour = Eco)) +
#   geom_point(aes(size = ratio_overlap_range)) +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_grid(ratio_gt_75~Eco) +
#   geom_smooth(method = "lm") +
#   scale_size(range = c(0, 2))
# 
# 
# ## plot max/min latitude of study area versus max/min latitude of species' range 
# 
# r <- area_ranges %>%
#   rename("xmin_range" = xmin, "xmax_range" = xmax) %>%
#   select(binomial, range_source, xmin_range, xmax_range)
# 
# cp <- area_cut_pols %>%
#   rename("ID" = Name, "species_studyid" = spcs_st, "range_source" = rng_src, 
#          "xmin_overlap" = xmin, "xmax_overlap" = xmax) %>%
#   select(ID, xmax_overlap, xmin_overlap, species_studyid, range_source) %>%
#   mutate(binomial = str_split_fixed(species_studyid, "\\_", 2)[,1])
# 
# test <- left_join(cp, r)
# 
# test %>%
#   left_join(v3, .) %>%
#   mutate(xmin_diff = xmin_range - xmin_overlap, 
#          xmax_diff = xmax_range - xmax_overlap) %>%
#   mutate(small_xmin_diff = ifelse(abs(xmin_diff) <= 20, "YES", "NO")) %>%
#   ggplot(aes(y = xmin_overlap, x = xmin_range, colour = small_xmin_diff)) +
#   geom_point(size = 0.1) +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~Param)
# 
# test %>%
#   left_join(v3, .) %>%
#   mutate(xmin_diff = xmin_range - xmin_overlap, 
#          xmax_diff = xmax_range - xmax_overlap) %>%
#   mutate(small_xmin_diff = ifelse(abs(xmin_diff) <= 20, "YES", "NO"),
#          small_xmax_diff = ifelse(abs(xmax_diff) <= 20, "YES", "NO")) %>%
#   mutate(Eco = ifelse(Eco %in% c("Mar", "Int"), "Mar", "Ter")) %>%
#   ggplot(aes(x = mean_cv_sppspecific, y = Rate, colour = small_xmax_diff)) +
#   geom_point(size = 0.1) +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(small_xmax_diff~Param) 
# 
# 
# 
# data2 <- list[[2]]
# data <- list[[3]]
# 
# which(data$binomial == "Cleisthenes pinetorum")
# 
# 
