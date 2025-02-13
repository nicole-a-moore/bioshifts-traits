## explore spatial scale 
## look at correlation between mean terrestrial climate velocities calculated at different scales 
library(stringr)
library(tidyverse)

## get list of csv files
files = list.files("data-raw/bioshiftsv3/Velocity_SA")
files = files[which(str_detect(files, "\\.csv"))]

## group by scale
files_25 <-  files[which(str_detect(files, "\\_25km.csv"))]
files_50 <-  files[which(str_detect(files, "\\_50km.csv"))]
files_110 <-  files[which(str_detect(files, "\\_110km.csv"))]

## subset to study areas that have been calculated across all 3 scales 
sa_25 <- paste0(str_split_fixed(files_25, "\\_", 3)[,1], "_", str_split_fixed(files_25, "\\_", 3)[,2])
sa_50 <- paste0(str_split_fixed(files_50, "\\_", 3)[,1], "_", str_split_fixed(files_50, "\\_", 3)[,2])
sa_110 <- paste0(str_split_fixed(files_110, "\\_", 3)[,1], "_", str_split_fixed(files_110, "\\_", 3)[,2])

sa_25 <- sa_25[which(sa_25 %in% sa_50 & sa_25 %in% sa_110)]
sa_50 <- sa_50[which(sa_50 %in% sa_25 & sa_50 %in% sa_110)]
sa_110 <- sa_110[which(sa_110 %in% sa_25 & sa_110 %in% sa_50)]

## read in 
file <- read.csv(paste0("data-raw/bioshiftsv3/Velocity_SA/", sa_25[1], "_25km.csv"))

## pick out columns to keep for terrestrial/marine study areas 
terr_cols <- c("ID", "v.lat.mean.mat", # mean latitudinal velocity of air temperature 
               "v.ele.mean.mat", # mean elevational velocity of air temperature
               "v.lat.sd.mat", # sd of both 
               "v.ele.sd.mat")
mar_cols <- c("ID", "v.lat.mean.sst",
              "v.lat.sd.sst")

cvs <- c()
for(f in 1:length(sa_25)) {
  
  ## read in a file 
  file <- read.csv(paste0("data-raw/bioshiftsv3/Velocity_SA/", sa_25[f], "_25km.csv"))
  
  ## if study is terrestrial, extract air temperature columns for ele and lat 
  if(any(str_detect(colnames(file), "mat"))) {
    
    file <- select(file, c("ID", "v.lat.mean.mat", "v.lat.sd.mat"))
    
    ## rearrange
    file <- gather(file, key = "cv_type", value = "val", 
                   c(v.lat.mean.mat, v.lat.sd.mat)) %>%
      mutate(Type = ifelse(str_detect(cv_type, "ele"), "ELE", 
                           "LAT"),
             Eco = "Ter") %>%
      mutate(measure = ifelse(str_detect(cv_type, "mean"), "mean_cv_studylevel", 
                              "sd_cv_studylevel")) %>%
      select(-cv_type) %>%
      spread(key = measure, value = val) %>%
      mutate(cv_res = "25km")
    
    ## bind 
    cvs <- rbind(cvs, file)
    
  }
  
  print(paste0("On file number ", f))
}

for(f in 1:length(sa_50)) {
  
  ## read in a file 
  file <- read.csv(paste0("data-raw/bioshiftsv3/Velocity_SA/", sa_50[f], "_50km.csv"))
  
  ## if study is terrestrial, extract air temperature columns for ele and lat 
  if(any(str_detect(colnames(file), "mat"))) {
    
    file <- select(file, c("ID", "v.lat.mean.mat", "v.lat.sd.mat"))
    
    ## rearrange
    file <- gather(file, key = "cv_type", value = "val", 
                   c(v.lat.mean.mat, v.lat.sd.mat)) %>%
      mutate(Type = ifelse(str_detect(cv_type, "ele"), "ELE", 
                           "LAT"),
             Eco = "Ter") %>%
      mutate(measure = ifelse(str_detect(cv_type, "mean"), "mean_cv_studylevel", 
                              "sd_cv_studylevel")) %>%
      select(-cv_type) %>%
      spread(key = measure, value = val) %>%
      mutate(cv_res = "50km")
    
    ## bind 
    cvs <- rbind(cvs, file)
    
  }
  
  print(paste0("On file number ", f))
}

for(f in 1:length(sa_110)) {
  
  ## read in a file 
  file <- read.csv(paste0("data-raw/bioshiftsv3/Velocity_SA/", sa_110[f], "_110km.csv"))
  
  ## if study is terrestrial, extract air temperature columns for ele and lat 
  if(any(str_detect(colnames(file), "mat"))) {
    
    file <- select(file, c("ID", "v.lat.mean.mat", "v.lat.sd.mat"))
    
    ## rearrange
    file <- gather(file, key = "cv_type", value = "val", 
                   c(v.lat.mean.mat, v.lat.sd.mat)) %>%
      mutate(Type = ifelse(str_detect(cv_type, "ele"), "ELE", 
                           "LAT"),
             Eco = "Ter") %>%
      mutate(measure = ifelse(str_detect(cv_type, "mean"), "mean_cv_studylevel", 
                              "sd_cv_studylevel")) %>%
      select(-cv_type) %>%
      spread(key = measure, value = val) %>%
      mutate(cv_res = "110km")
    
    ## bind 
    cvs <- rbind(cvs, file)
    
  }
  
  print(paste0("On file number ", f))
}

## plot against each other 
cvs_wide = cvs %>%
  select(- "sd_cv_studylevel") %>%
  spread(key = "cv_res", value = c("mean_cv_studylevel"))

colnames(cvs_wide)[4:6] <- paste0("cv_res_", colnames(cvs_wide[4:6]))

cvs_wide %>%
  ggplot(aes(x = cv_res_25km, y = cv_res_110km)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

cvs_wide %>%
  ggplot(aes(x = cv_res_50km, y = cv_res_110km)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

cvs_wide %>%
  ggplot(aes(x = cv_res_25km, y = cv_res_50km)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)


## look at dispersal data
## read in data 
dd <- read.csv("~/Documents/dispersal-limitation/data-processed/v3_potential-dispersal-rate.csv")

## filter to only latitude 
dd <- filter(dd, Type == "LAT")

## filter to leading edge shifts with positive climate velocity 
dd <- filter(dd, !Param %in% c("O", "TE") & ClimVeloTKmY_spp >= 0)

## get rid of range contractions that are father than 1 sd from the mean shift 
dd <- filter(dd, Rate >= (mean(dd$ShiftKmY) - sd(dd$ShiftKmY)))

hist(dd$ClimVeloTKmY_spp)
hist(dd$ShiftKmY)

dd <- left_join(dd, cvs_wide, by = c("ID", "Eco", "Type"))

dd <- dd %>%
  mutate(LimitingRate = ifelse(DispersalPotentialKmY <= ClimVeloTKmY_spp,
                               DispersalPotentialKmY,
                               ClimVeloTKmY_spp)) %>%
  mutate(LimitingRate_25 = ifelse(DispersalPotentialKmY <= cv_res_25km,
                                   DispersalPotentialKmY,
                                  cv_res_25km)) %>%
  mutate(LimitingRate_50 = ifelse(DispersalPotentialKmY <= cv_res_50km,
                                  DispersalPotentialKmY,
                                  cv_res_50km)) %>%
  mutate(LimitingRate_110 = ifelse(DispersalPotentialKmY <= cv_res_110km,
                               DispersalPotentialKmY,
                               cv_res_110km)) %>%
  mutate(cv_relevant = ifelse(DispersalDistanceKm <= 25, cv_res_25km,
                              ifelse(DispersalDistanceKm <= 50, cv_res_50km,
                                     ifelse(DispersalDistanceKm > 50, cv_res_110km,
                                            NA)))) %>%
  mutate(LimitingRate_scaled = ifelse(DispersalPotentialKmY <= cv_relevant,
                                      DispersalPotentialKmY,
                                      cv_relevant)) %>%
  mutate(what_is_limiting = ifelse(DispersalPotentialKmY == LimitingRate, "Dispersal", "Climate")) %>%
  mutate(what_is_limiting_25 = ifelse(DispersalPotentialKmY == LimitingRate_25, "Dispersal", "Climate")) %>%
  mutate(what_is_limiting_50 = ifelse(DispersalPotentialKmY == LimitingRate_50, "Dispersal", "Climate")) %>%
  mutate(what_is_limiting_110 = ifelse(DispersalPotentialKmY == LimitingRate_110, "Dispersal", "Climate")) %>%
  mutate(what_is_limiting_scaled = ifelse(DispersalPotentialKmY == LimitingRate_scaled, "Dispersal", "Climate")) %>%
  mutate(colour = ifelse(what_is_limiting == "Climate", ClimVeloTKmY_spp, NA)) %>%
  mutate(colour25 = ifelse(what_is_limiting_25 == "Climate", cv_res_25km, NA)) %>%
  mutate(colour50 = ifelse(what_is_limiting_50 == "Climate", cv_res_50km, NA)) %>%
  mutate(colour110 = ifelse(what_is_limiting_110 == "Climate", cv_res_110km, NA)) %>%
  mutate(colourscaled = ifelse(what_is_limiting_scaled == "Climate", LimitingRate_scaled, NA))

dd %>%
  ggplot(aes(x = LimitingRate, y = ShiftKmY, colour = ClimVeloTKmY_spp, shape = group)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(dd, is.na(colour), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(dd, is.na(colour), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change(km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") 


dd %>%
  filter(cv_res_25km >= 0) %>%
  ggplot(aes(x = LimitingRate_25, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(dd, is.na(colour25), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate_25, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(dd, is.na(colour25), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate_25, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change(km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_smooth(method = "lm")


dd %>%
  filter(cv_res_110km >= 0) %>%
  ggplot(aes(x = LimitingRate_110, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(dd, is.na(colour110), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate_110, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(dd, is.na(colour110), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate_110, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change(km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_smooth(method = "lm")

## choose resolution of climate velocity based on dispersal distance 
dd %>%
  filter(LimitingRate_scaled >= 0) %>%
  ggplot(aes(x = LimitingRate_scaled, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(dd, is.na(colourscaled), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate_scaled, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(dd, is.na(colourscaled), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate_scaled, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change(km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_smooth(method = "lm")

dd %>%
  filter(LimitingRate_scaled >= 0) %>%
  ggplot(aes(x = LimitingRate_scaled, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  geom_point(data = filter(dd, is.na(colourscaled), group == "Plants"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 2,
             aes(x = LimitingRate_scaled, y = ShiftKmY, shape = group)) +
  geom_point(data = filter(dd, is.na(colourscaled), group == "Birds"), colour = "black", inherit.aes = FALSE,
             fill = "transparent", pch = 1,
             aes(x = LimitingRate_scaled, y = ShiftKmY, shape = group)) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change(km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_smooth(method = "lm")

dd %>%
  filter(cv_res_110km >= 0) %>%
  filter(cv_res_110km != LimitingRate_110) %>%
  ggplot(aes(x = LimitingRate_110, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change(km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_smooth(method = "lm")

dd %>%
  filter(cv_res_110km >= 0) %>%
  filter(cv_res_110km != LimitingRate_110) %>%
  ggplot(aes(x = cv_res_110km, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change(km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_smooth(method = "lm")

dd %>%
  filter(cv_relevant >= 0) %>%
  filter(cv_relevant != LimitingRate_scaled) %>%
  ggplot(aes(x = LimitingRate_scaled, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change(km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_smooth(method = "lm")

dd %>%
  filter(cv_relevant >= 0) %>%
  filter(cv_relevant != LimitingRate_scaled) %>%
  ggplot(aes(x = cv_relevant, y = ShiftKmY, colour = ClimVeloTKmY_spp)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "grey", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(2 , "lines"), 
        legend.position = "none") +
  labs(x = "Minimum of potential dispersal rate\nand rate of climate change(km/y)",
       y = "Observed range shift rate (km/y)", 
       colour = 'Mean\nclimate\nvelocity\n(km/y)', 
       shape = "") +
  scale_y_continuous(limits = c(-6, 26), expand = c(0,0.5)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5,
                         na.value = "black") +
  geom_smooth(method = "lm")





