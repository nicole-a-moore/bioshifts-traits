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
cvs <- unique(cvs)

## NOTE TO SELF: some obs getting duplicated - deal with the source of the problem later 

cvs %>%
  group_by(species_studyid, cv_type) %>%
  tally()
## to avoid creating duplicate observations for species with multiple realized range sources, make sure each species has only 1 row in the cv dataframe 
## for now, just take the range from the range source that comes first in alphabetical order
cvs <- cvs %>%
  group_by(species_studyid, cv_type) %>%
  arrange(range_source, .by_group = TRUE) %>%
  filter(row_number()==1) %>%
  ungroup()

## check that there is 1 row per species per type of climate velocity data (e.g., ele, lat, sst)
length(unique(paste(cvs$species_studyid, cvs$cv_type))) == nrow(cvs)

## keep only important columns 
cvs <- cvs %>%
  select(species_studyid, cv_type, range_source, mean_cv_sppspecific, 
         sd_cv_sppspecific, mean_cv_studylevel, sd_cv_studylevel)

## get rid of sst_gVel 
cvs <- filter(cvs, !cv_type == "sst_gVel")

## create column that specifies whether the cv applies to a lat or ele shift 
cvs$Type = ifelse(cvs$cv_type == "mat_gVelEle", "ELE", "LAT")

## join to bioshifts 
final <- inner_join(cvs, v3)

length(unique(final$species_name)) # 7226 spp

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
  filter(range_source %in% c("IUCN", "Aquamaps")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>%
  spread(key = "range_source", value = "mean_cv_sppspecific") %>%
  ggplot(aes(x = IUCN, y = Aquamaps)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
  
cvs %>%
  filter(range_source %in% c("IUCN", "BIEN")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>%
  spread(key = "range_source", value = "mean_cv_sppspecific") %>%
  ggplot(aes(x = IUCN, y = BIEN)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

cvs %>%
  filter(range_source %in% c("IUCN", "BOTW")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>%
  spread(key = "range_source", value = "mean_cv_sppspecific") %>%
  ggplot(aes(x = IUCN, y = BOTW)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

cvs %>%
  filter(range_source %in% c("IUCN", "GIFT")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>% View
  spread(key = "range_source", value = "mean_cv_sppspecific") %>%
  ggplot(aes(x = IUCN, y = GIFT)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
  
cvs %>%
    filter(range_source %in% c("IUCN", "GARD")) %>%
    group_by(species_studyid, cv_type) %>%
    mutate(n = n()) %>%
    filter(n >= 2) %>% 
    spread(key = "range_source", value = "mean_cv_sppspecific") %>%
    ggplot(aes(x = IUCN, y = GARD)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1)

cvs %>%
  filter(range_source %in% c("IUCN", "Fishmap")) %>%
  group_by(species_studyid, cv_type) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>% 
  spread(key = "range_source", value = "mean_cv_sppspecific") %>%
  ggplot(aes(x = IUCN, y = "Fishmap")) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) 
  

