## join new species-specific climate velocity data to bioshifts v3
library(tidyverse)
library(stringr)

######################################
####   JOIN CV DATA TO BIOSHIFTS  ####
######################################
###### prep bioshifts ######
## read in bioshifts v3 
v3 = read.csv("data-raw/bioshiftsv3/BIOSHIFTS_v3.csv")
v3$species_name = str_replace_all(v3$sp_name_checked, "\\_", " ")
v3$species_name

## make column that will act as joining key (species name_study id)
v3$species_studyid <- paste0(v3$species_name, "_", v3$ID)

###### prep cv data ######
## collate all cv data
files = list.files("data-processed/new-climate-velocities", full.names = TRUE)
cvs <- c()
for(i in 1:length(files)) {
  cvs <- rbind(cvs, read.csv(files[i]))
}
head(cvs)

## make sure only 1 obs per study x species x range source x cv type 
cvs %>%
  group_by(species_studyid, cv_type) %>%
  tally()

## to avoid creating duplicate observations for species with multiple realized range sources, make sure each species has only 1 row in the cv dataframe 
## create a list of priority for range sources:
range_sources <- c("IUCN", "BOTW", "GARD", "BIEN", "Fishmap", "Aquamaps", "Butterfly Maps", "GIFT", "GBIF occurrence")

cvs %>%
  ggplot(aes(x = range_source)) + geom_histogram(stat = "count")

cvs <- cvs %>%
  group_by(species_studyid, cv_type) %>%
  arrange(match(range_source, range_sources),.by_group = TRUE) %>%
  filter(row_number()==1) %>%
  ungroup()

cvs %>%
  ggplot(aes(x = range_source)) + geom_histogram(stat = "count")

## check that there is 1 row per species per type of climate velocity data (e.g., ele, lat, sst)
length(unique(paste(cvs$species_studyid, cvs$cv_type))) == nrow(cvs)

## keep only important columns 
cvs <- cvs %>%
  select(species_studyid, cv_type, range_source, mean_cv_sppspecific, 
         sd_cv_sppspecific, mean_cv_studylevel, sd_cv_studylevel)

## get rid of latitudinal velocities calculated at 1km  
cvs <- filter(cvs, !cv_type == "mat_gVelLat",
              !cv_type == "map_gVelLat")
  
## create column that specifies whether the cv applies to a lat or ele shift 
cvs$Type = ifelse(cvs$cv_type %in% c("mat_gVelEle", "map_gVelEle"), "ELE", "LAT")

## write
write.csv(cvs, "data-processed/new-cvs_preliminary.csv", row.names = FALSE)

## join to bioshifts 
final <- inner_join(cvs, v3)

length(unique(final$species_name)) # 10688 spp

write.csv(final, "data-processed/v3_new-cvs_preliminary.csv", row.names = FALSE)

final <- read.csv("data-processed/v3_new-cvs_preliminary.csv")

## how different are the estimates that use different range sources?
files = list.files("data-processed/new-climate-velocities", full.names = TRUE)
cvs <- c()
for(i in 1:length(files)) {
  cvs <- rbind(cvs, read.csv(files[i]))
}
head(cvs)
cvs <- unique(cvs)

## start with IUCN versus Aquamaps
cvs %>%
  select(-sd_cv_sppspecific) %>%
  filter(range_source %in% c("IUCN", "Aquamaps")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>%
  spread(key = "range_source", value = "mean_cv_sppspecific") %>%
  ggplot(aes(x = IUCN, y = Aquamaps)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
  
cvs %>%
  select(-sd_cv_sppspecific) %>%
  filter(range_source %in% c("IUCN", "BIEN")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>%
  spread(key = "range_source", value = "mean_cv_sppspecific") %>%
  ggplot(aes(x = IUCN, y = BIEN)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

cvs %>%
  select(-sd_cv_sppspecific) %>%
  filter(range_source %in% c("IUCN", "BOTW")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>%
  spread(key = "range_source", value = "mean_cv_sppspecific") %>% 
  ggplot(aes(x = IUCN, y = BOTW)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

cvs %>%
  select(-sd_cv_sppspecific) %>%
  filter(range_source %in% c("IUCN", "GIFT")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>%  
  spread(key = "range_source", value = "mean_cv_sppspecific") %>% 
  ggplot(aes(x = IUCN, y = GIFT)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
  
cvs %>%
  select(-sd_cv_sppspecific) %>%
  filter(range_source %in% c("IUCN", "GARD")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>% 
  spread(key = "range_source", value = "mean_cv_sppspecific") %>%
  ggplot(aes(x = IUCN, y = GARD)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

cvs %>%
  select(-sd_cv_sppspecific) %>%
  filter(range_source %in% c("IUCN", "Fishmap")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>% 
  spread(key = "range_source", value = "mean_cv_sppspecific") %>%
  ggplot(aes(x = IUCN, y = Fishmap)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) 
  

cvs %>%
  select(-sd_cv_sppspecific) %>%
  filter(range_source %in% c("IUCN", "GBIF occurrence")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>% 
  spread(key = "range_source", value = "mean_cv_sppspecific") %>%
  ggplot(aes(x = IUCN, y = `GBIF occurrence`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) 
## GBIF occurrence and IUCN are the most different 

## is the bias directional?
cvs %>%
  select(-sd_cv_sppspecific) %>%
  filter(range_source %in% c("IUCN", "GBIF occurrence")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>% 
  spread(key = "range_source", value = "mean_cv_sppspecific") %>%
  mutate(diff = IUCN - `GBIF occurrence`) %>%
  ggplot(aes(x = diff)) +
  geom_histogram()
## no
  
final %>%
  ggplot(aes(x = mean_cv_sppspecific, y = mean_cv_studylevel)) +
  geom_point() +
  facet_wrap(~cv_type)
  

final %>%
  ggplot(aes(x = mean_cv_sppspecific, y = Rate)) +
  geom_point() +
  geom_abline(intercept = 0, slope =1) +
  facet_wrap(~Param)

