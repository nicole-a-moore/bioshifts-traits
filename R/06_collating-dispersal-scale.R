## bring together empirical dispersal distance data
library(tidyverse)
library(readxl)
library(readr)
library(taxadb)
library(parallel)
library(pbapply)
library(traitdataform)
library(data.table)

## read in list of all bioshifts species 
sp <- read_csv("data-raw/splist.csv")

## read function to harmonize taxonomy 
source("R/harmonize.R")

cols_to_keep <- c("reported_name","reported_name_fixed", "scientificName", "kingdom", "phylum",
                  "class", "order", "family", "db", "db_code")

## Paradis et al 2002:
## read in data
par = read_csv("data-raw/dispersal/Paradis_et_al_2002.csv") 
colnames(par) <- str_replace_all(colnames(par), "\\ ", "_")
length(unique(par$Species)) #75 spp

par_harm <- harmonize(par$Species)

notfound <- filter(par_harm, is.na(db_code))

## rename columns 
par <- par %>%
  rename("reported_name" = Species) %>%
  mutate(reported_name_fixed = reported_name)

par <- left_join(par, par_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## reorganize
par_dd <- par %>%
  select(all_of(cols_to_keep), 
         AM_breeding, AM_natal) %>%
  gather(key = "Field", value = "DispersalDistance", c(AM_breeding, AM_natal)) %>%
  mutate(Field = ifelse(Field == "AM_breeding", "ArithmeticMeanBreedingDispersal",
                        ifelse(Field == "AM_natal", 
                               "ArithmeticMeanNatalDispersal",
                               NA))) %>%
  mutate(Code = "MeanDispersalDistance") %>%
  filter(!is.na(DispersalDistance)) %>%
  mutate(Sex = NA, Source = NA, Unit = "km",
         ObservationType = "ring recoveries",
         Database = "Paradis et al. 2002") 

## check how many species in bioshifts 
length(which(unique(par_dd$scientificName) %in% unique(sp$scientificName))) ## 74
par_sp <- unique(par_dd$scientificName)[which(unique(par_dd$scientificName) %in% unique(sp$scientificName))]


## data on dispersal scale:
## comte and olden - fish mean/max dispersal based on mark recapture, telemetry, etc.
co <- read_excel("data-raw/dispersal/ComteOlden_FishFisheries_SOM.xlsx")
colnames(co) <- str_replace_all(colnames(co), "\\ ", "_")
length(unique(co$Species)) #142 spp

co_harm <- harmonize(co$Species)

notfound <- filter(co_harm, is.na(db_code))

## rename columns 
co <- co %>%
  rename("reported_name" = Species) %>%
  mutate(reported_name_fixed = reported_name)

co <- left_join(co, co_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## reformat the data
co_dd <- co %>%
  select(all_of(cols_to_keep), 
         `Mean_dispersal_distance_(m)`, 
         `Maximum_dispersal_distance_(m)`, 
         `Median_dispersal_distance_(m)`, 
         Type, Source) %>%
  gather(key = "Field", value = "DispersalDistance", c(`Mean_dispersal_distance_(m)`, 
                                                       `Maximum_dispersal_distance_(m)`, 
                                                       `Median_dispersal_distance_(m)`)) %>%
  mutate(Code = ifelse(Field == "Mean_dispersal_distance_(m)",
                       "MeanDispersalDistance", 
                       ifelse(Field == "Maximum_dispersal_distance_(m)",
                              "MaxDispersalDistance", 
                              ifelse(Field == "Median_dispersal_distance_(m)",
                                     "MedianDispersalDistance", 
                                     NA)))) %>%
  mutate(Database = "Comte & Olden", Unit = "m", Sex = NA) %>%
  filter(!is.na(DispersalDistance)) %>%
  rename("ObservationType" = Type) %>%
  mutate(ObservationType = ifelse(ObservationType == "Radio-tracking", 
                                  "radio-tracking", 
                                  ifelse(ObservationType == "Mark-Recapture", 
                                  "mark-release-recapture",
                                  ObservationType)))

## check how many species in bioshifts 
length(which(unique(co_dd$scientificName) %in% unique(sp$scientificName))) ## 28
co_sp <- unique(co_dd$scientificName)[which(unique(co_dd$scientificName) %in% unique(sp$scientificName))]


## Sutherland 2000 - birds and mammals natal and breeding dispersal
suth_mamm <- read_csv("data-raw/dispersal/Sutherland_2000_mammals.csv")
colnames(suth_mamm) <- str_replace_all(colnames(suth_mamm), "\\ ", "_")
suth_mamm$Species = ifelse(suth_mamm$Species == "", NA, suth_mamm$Species)

## fill
suth_mamm = fill(suth_mamm, Species, .direction = "down")

sm_harm <- harmonize(suth_mamm$Species)

notfound <- filter(sm_harm, is.na(db_code))

## rename columns 
suth_mamm <- suth_mamm %>%
  rename("reported_name" = Species) %>%
  mutate(reported_name_fixed = reported_name)

suth_mamm <- left_join(suth_mamm, sm_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## Sutherland mammals:
suth_mamm_dd <- suth_mamm %>%
  select(all_of(cols_to_keep), 
         `Natal_dispersal_median_distance_(km)`,
         `Natal_dispersal_maximum_distance_(km)`, Obs_type, Source) %>%
  mutate(ObservationType = "natal dispersal") %>%
  gather(key = "Field", value = "DispersalDistance", c(`Natal_dispersal_median_distance_(km)`,
                                                       `Natal_dispersal_maximum_distance_(km)`)) %>%
  mutate(Code = ifelse(Field == "Natal_dispersal_median_distance_(km)",
                       "MedianDispersalDistance", 
                       ifelse(Field == "Natal_dispersal_maximum_distance_(km)",
                              "MaxDispersalDistance", 
                              NA))) %>%
  mutate(Database = "Sutherland (mammals)", Unit = "km") %>%
  filter(!is.na(DispersalDistance), DispersalDistance != "...") %>%
  mutate(Sex = str_split_fixed(DispersalDistance, " ", 2)[,2], 
         DispersalDistance = str_split_fixed(DispersalDistance, " ", 2)[,1]) %>%
  select(-Obs_type) %>%
  mutate()

## check how many species in bioshifts 
length(which(unique(suth_mamm_dd$scientificName) %in% unique(sp$scientificName))) ## 16
suth_mamm_sp <- unique(suth_mamm_dd$scientificName)[which(unique(suth_mamm_dd$scientificName) %in% unique(sp$scientificName))]


suth_bird <- read_csv("data-raw/dispersal/Sutherland_2000_birds.csv")
colnames(suth_bird) <- str_replace_all(colnames(suth_bird), "\\ ", "_")
length(unique(suth_bird$Species)) #78 spp
suth_bird$Species = ifelse(suth_bird$Species == "", NA, suth_bird$Species)

## fill
suth_bird = fill(suth_bird, Species, .direction = "down")

## fix taxonomy
suth_bird$reported_name_fixed = suth_bird$Species
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Accipiter cooperi")] <- "Accipiter cooperii"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Bonasa bonasia")] <- "Tetrastes bonasia"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Dendragapus canadensis")] <- "Falcipennis canadensis"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Actitis macularia")] <- "Actitis macularius"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Speotyto cunicularia")] <- "Athene cunicularia"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Delichon urbica")] <- "Delichon urbicum"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Parus atricapillus")] <- "Poecile atricapillus" 
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Parus palustris")] <- "Poecile palustris"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Sitta europea")] <- "Sitta europaea"
suth_bird$reported_name_fixed[which(suth_bird$reported_name_fixed == "Carduelis chloris")] <- "Chloris chloris"

sb_harm <- harmonize(suth_bird$reported_name_fixed)

## fix taxonomy manually
notfound <- filter(sb_harm, is.na(db_code))

## rename columns 
suth_bird <- suth_bird %>%
  rename("reported_name" = Species) 

suth_bird <- left_join(suth_bird, sb_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## Sutherland birds:
suth_bird_dd <- suth_bird %>%
  select(all_of(cols_to_keep), 
         `Natal_dispersal_median_distance_(km)`,
         `Natal_dispersal_maximum_distance_(km)`, Obs_type,
         Source) %>%
  mutate(ObservationType = "natal dispersal") %>%
  gather(key = "Field", value = "DispersalDistance", c(`Natal_dispersal_median_distance_(km)`,
                                                       `Natal_dispersal_maximum_distance_(km)`)) %>%
  mutate(Code = ifelse(Field == "Natal_dispersal_median_distance_(km)",
                       "MedianDispersalDistance", 
                       ifelse(Field == "Natal_dispersal_maximum_distance_(km)",
                              "MaxDispersalDistance", 
                              NA))) %>%
  mutate(Database = "Sutherland (birds)", Unit = "km") %>%
  filter(!is.na(DispersalDistance), DispersalDistance != "...") %>%
  mutate(Sex = str_split_fixed(DispersalDistance, " ", 2)[,2], 
         DispersalDistance = str_split_fixed(DispersalDistance, " ", 2)[,1]) %>%
  select(-Obs_type)

## check how many species in bioshifts 
length(which(unique(suth_bird_dd$scientificName) %in% unique(sp$scientificName))) ## 72
suth_bird_sp <- unique(suth_bird_dd$scientificName)[which(unique(suth_bird_dd$scientificName) %in% unique(sp$scientificName))]


## Whitmee & Orme 2013 - mammal natal dispersal distances
## strong phylogenetic signal
## dispersal distance per generation can be estimated for mammals with comparatively little data availability!!! 
## Home range area, geographic range size and body mass 
wo <- read_csv("data-raw/dispersal/Whitmee_and_Orme_2013.csv")
colnames(wo) <- str_replace_all(colnames(wo), "\\ ", "_")
length(unique(wo$Species)) #104 spp

wo_harm <- harmonize(wo$Species)

## fix some taxonomy
notfound <- filter(wo_harm, is.na(db_code))
notfound$species[which(notfound$species == "Spermophilus saturatus")] <- "Callospermophilus saturatus"
wo_harm_notfound <- harmonize(notfound$species)

wo$reported_name_fixed = wo$Species
wo$reported_name_fixed[which(wo$reported_name_fixed == "Spermophilus saturatus")] <- "Callospermophilus saturatus"

## bind:
wo_harm <-filter(wo_harm, !is.na(db_code)) %>%
  rbind(., wo_harm_notfound)

## rename columns 
wo <- wo %>%
  rename("reported_name" = Species)

wo <- left_join(wo, wo_harm, by = c("reported_name" = "species")) %>%
  unique()

## Whitmee & Orme
wo_dd <- wo %>%
  select(all_of(cols_to_keep), Value, Units, Sex, Measure, Ref_no) %>%
  mutate(Code = ifelse(Measure == "Mean",
                       "MeanDispersalDistance", 
                       ifelse(Measure == "Maximum",
                              "MaxDispersalDistance", 
                              ifelse(Measure == "Median",
                                     "MedianDispersalDistance", 
                                     NA)))) %>%
  mutate(Unit = ifelse(str_detect(.$Units, "km"), "km",
                       ifelse(str_detect(.$Units, "Metres"), "m", NA))) %>%
  mutate(Database = "Whitmee & Orme 2013", ObservationType = "individual movement distances", Field = "Value") %>%
  rename("DispersalDistance" = Value, "Source" = Ref_no) %>%
  filter(!is.na(DispersalDistance)) %>%
  select(-Measure, -Units) 

## get rid of species in Sutherland to avoid duplicating observations
length(which(unique(wo_dd$scientificName) %in% unique(suth_mamm_dd$scientificName))) ## 49
length(which(unique(wo_dd$scientificName) %in% unique(suth_bird_dd$scientificName))) ## 0

wo_dd <- filter(wo_dd, !scientificName %in% c(suth_mamm_dd$scientificName, suth_bird_dd$scientificName))

## check how many species in bioshifts 
length(which(unique(wo_dd$scientificName) %in% unique(sp$scientificName))) ## 5
wo_sp <- unique(wo_dd$scientificName)[which(unique(wo_dd$scientificName) %in% unique(sp$scientificName))]


## Shanks 2009/2003:
shanks9 <- read_csv("data-raw/dispersal/Shanks_2009.csv")
colnames(shanks9) <- str_replace_all(colnames(shanks9), "\\ ", "_")
shanks9$genus_species <- paste(shanks9$Genus, shanks9$Species, sep = " ")
length(unique(shanks9$genus_species)) #44 spp

s9_harm <- harmonize(shanks9$genus_species)

notfound <- filter(s9_harm, is.na(db_code))

## rename columns 
shanks9 <- shanks9 %>%
  rename("reported_name" = genus_species) %>%
  mutate(reported_name_fixed = reported_name)

shanks9 <- left_join(shanks9, s9_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## Shanks 2009
shanks9_dd <- shanks9 %>%
  select(all_of(cols_to_keep),
         Dispersal_distance, source, References) %>%
  mutate(Unit = str_split_fixed(.$Dispersal_distance, " ", 2)[,2],
         DispersalDistance = str_split_fixed(.$Dispersal_distance, " ", 2)[,1]) %>%
  mutate(Unit = ifelse(Unit == "m1", "m", Unit), 
         DispersalDistance = str_replace_all(DispersalDistance, "\\<", "")) %>%
  filter(Unit %in% c("km", "m")) %>%
  mutate(Code = ifelse(str_detect(.$DispersalDistance, "\\–"),  "DispersalDistanceRange",
                       "DispersalDistance")) %>%
  mutate(Database = "Shanks 2009", Sex = NA, Field = "Dispersal_distance", ObservationType = source) %>%
  select(-Dispersal_distance, -source) %>%
  filter(!is.na(DispersalDistance)) %>%
  rename("Source" = References)

## check how many species in bioshifts 
length(which(unique(shanks9_dd$scientificName) %in% unique(sp$scientificName))) ## 7
shanks9_sp <- unique(shanks9_dd$scientificName)[which(unique(shanks9_dd$scientificName) %in% unique(sp$scientificName))]


shanks3 <- read_csv("data-raw/dispersal/Shanks_et_al_2003.csv")
colnames(shanks3) <- str_replace_all(colnames(shanks3), "\\ ", "_")
length(unique(shanks3$Genus_species)) #32 spp

s3_harm <- harmonize(shanks3$Genus_species)

notfound <- filter(s3_harm, is.na(db_code))

## rename columns 
shanks3 <- shanks3 %>%
  rename("reported_name" = Genus_species) %>%
  mutate(reported_name_fixed = reported_name)

shanks3 <- left_join(shanks3, s3_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## Shanks 2003
shanks3_dd <- shanks3 %>%
  select(all_of(cols_to_keep),
         `Realized_dispersal_distance_(mean)`, source, References) %>%
  mutate(Unit = str_split_fixed(.$`Realized_dispersal_distance_(mean)`, " ", 2)[,2],
         DispersalDistance = str_split_fixed(.$`Realized_dispersal_distance_(mean)`, " ", 2)[,1]) %>%
  mutate(DispersalDistance = str_replace_all(DispersalDistance, "\\<", "")) %>%
  mutate(Code = ifelse(str_detect(.$DispersalDistance, "\\-"),  
                       "DispersalDistanceRange",
                       "DispersalDistance")) %>%
  mutate(Database = "Shanks et al. 2003", Sex = NA, Field = "Realized_dispersal_distance_(mean)",
         ObservationType = source) %>%
  select(-`Realized_dispersal_distance_(mean)`, -source) %>%
  filter(!is.na(DispersalDistance))%>%
  rename("Source" = References) %>%
  unique() 

## check how many species in bioshifts 
length(which(unique(shanks3_dd$scientificName) %in% unique(sp$scientificName))) ## 4
shanks3_sp <- unique(shanks3_dd$scientificName)[which(unique(shanks3_dd$scientificName) %in% unique(sp$scientificName))]


## Santini et al. 2013:
## trapping, radio tracking and parent son distance identified through genetics
## read in data
sant = read_csv("data-raw/dispersal/Santini_et_al_2013.csv")
colnames(sant) <- str_replace_all(colnames(sant), "\\ ", "_")
length(unique(sant$Species)) #174 spp

sant_harm <- harmonize(sant$Species)

notfound <- filter(sant_harm, is.na(db_code))

## rename columns 
sant <- sant %>%
  rename("reported_name" = Species) %>%
  mutate(reported_name_fixed = reported_name)

sant <- left_join(sant, sant_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## Santini
## must clean out sex from dispersal value
## clean out sex by gathering
sant_dd <- sant %>%
  select(all_of(cols_to_keep), Mean, Median, Maximum, n, Type, Source) 

d_max = as.data.frame(str_split_fixed(sant_dd$Maximum, " ", 4)) %>%
  mutate(reported_name = sant_dd$reported_name, 
         Source = sant_dd$Source,
         n = sant_dd$n,
         Type = sant_dd$Type) %>%
  unique() %>%
  mutate(D1 = paste(V1, V2, sep = " "), D2 = paste(V3, V4, sep = " ")) %>%
  gather(key = "est", value = "Maximum", c(D1, D2)) %>%
  select(-c(V1, V2, V3, V4, est)) %>%
  filter(!is.na(Maximum)) %>%
  mutate(Sex = str_split_fixed(Maximum, " ", 2)[,2], 
         Maximum = str_split_fixed(Maximum, " ", 2)[,1]) %>%
  filter(Maximum != "", Maximum != "f", Maximum != "m",
         Maximum != "mf") %>%
  mutate(Maximum = paste(Maximum, Sex, sep = " ")) %>%
  select(-Sex)

d_mean = as.data.frame(str_split_fixed(sant_dd$Mean, " ", 8)) %>%
  mutate(reported_name = sant_dd$reported_name,
         Source = sant_dd$Source,
         n = sant_dd$n,
         Type = sant_dd$Type) %>%
  unique() %>%
  mutate(D1 = paste(V1, V2, sep = " "), D2 = paste(V3, V4, sep = " "),  D3 = paste(V5, V6, sep = " "),
         D4 = paste(V7, V8, sep = " ")) %>%
  gather(key = "est", value = "Mean", c(D1, D2, D3, D4)) %>%
  select(-c(V1, V2, V3, V4, V5, V6, V7, V8,  est)) %>%
  filter(!is.na(Mean)) %>%
  mutate(Sex = str_split_fixed(Mean, " ", 2)[,2], 
         Mean = str_split_fixed(Mean, " ", 2)[,1]) %>%
  filter(Mean != "", Mean != "f", Mean != "m",
         Mean != "mf")%>%
  mutate(Mean = paste(Mean, Sex, sep = " ")) %>%
  select(-Sex)

d_med = as.data.frame(str_split_fixed(sant_dd$Median, " ", 8)) %>%
  mutate(reported_name = sant_dd$reported_name,
         Source = sant_dd$Source,
         n = sant_dd$n,
         Type = sant_dd$Type) %>%
  unique() %>%
  mutate(D1 = paste(V1, V2, sep = " "), D2 = paste(V3, V4, sep = " "),  D3 = paste(V5, V6, sep = " "),
         D4 = paste(V7, V8, sep = " ")) %>%
  gather(key = "est", value = "Median", c(D1, D2, D3, D4)) %>%
  select(-c(V1, V2, V3, V4, V5, V6, V7, V8,  est)) %>%
  filter(!is.na(Median)) %>%
  mutate(Sex = str_split_fixed(Median, " ", 2)[,2], 
         Median = str_split_fixed(Median, " ", 2)[,1]) %>%
  filter(Median != "", Median != "f", Median != "m",
         Median != "mf") %>%
  mutate(Median = paste(Median, Sex, sep = " ")) %>%
  select(-Sex)

sant_dd = sant_dd %>%
  select(-Mean, -Maximum, -Median) %>%
  unique() %>%
  left_join(., d_max) %>%
  left_join(., d_med) %>%
  left_join(., d_mean) %>%
  gather(key = "Field", value = "DispersalDistance", 
         c(Mean, Median, Maximum)) %>%
  unique() %>%
  filter(!is.na(DispersalDistance)) %>%
  mutate(Code = ifelse(Field == "Mean", 
                       "MeanDispersalDistance", 
                       ifelse(Field == "Maximum",
                              "MaxDispersalDistance", 
                              ifelse(Field == "Median",
                                     "MedianDispersalDistance", 
                                     NA)))) %>%
  mutate(Unit = "km", 
         Database = "Santini et al. 2013") %>%
  filter(!is.na(DispersalDistance))%>%
  rename("ObservationType" = Type) %>%
  mutate(Sex = str_split_fixed(DispersalDistance, " ", 2)[,2], 
         DispersalDistance = str_split_fixed(DispersalDistance, " ", 2)[,1]) %>%
  select(-n) %>%
  unique() %>%
  mutate(ObservationType = ifelse(ObservationType == "H", 
                                  "homing study",
                                  ifelse(ObservationType == "D",
                                         "natural dispersal", 
                                         ifelse(ObservationType == "P",
                                                "post release dispersal",
                                                ifelse(ObservationType == "A",
                                                       "unknown",
                                                       NA)))))

## Santini took from Sutherland 
## exclude observations from Sutherland with same Code and DispersalDistance as Santini since has more info
key = suth_mamm_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(data_id = paste(scientificName, Code, DispersalDistance, sep = "_")) %>%
  select(data_id)

sant_dd <- sant_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  filter(!paste(scientificName, Code, DispersalDistance, sep = "_") %in% key$data_id)

## check how many species in bioshifts 
length(which(unique(sant_dd$scientificName) %in% unique(sp$scientificName))) ## 26
sant_sp <- unique(sant_dd$scientificName)[which(unique(sant_dd$scientificName) %in% unique(sp$scientificName))]


## Bowman et al  2002
bowman <- read_csv("data-raw/dispersal/Bowman_et_al_2002.csv")
colnames(bowman) <- str_replace_all(colnames(bowman), "\\ ", "_")
length(unique(bowman$Genus_species)) #21 spp

bowman_harm <- harmonize(bowman$Genus_species)

notfound <- filter(bowman_harm, is.na(db_code))

## rename columns 
bowman <- bowman %>%
  rename("reported_name" = Genus_species) %>%
  mutate(reported_name_fixed = reported_name)

bowman <- left_join(bowman, bowman_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## Bowman
bowman_dd <- bowman %>%
  select(all_of(cols_to_keep),
         `Distance_(km)`, Reference) %>%
  rename("DispersalDistance" = `Distance_(km)`, "Source" = Reference) %>%
  mutate(Sex = NA, ObservationType = "natal dispersal", Unit = "km", Field = "Distance_(km)",
         Code = "DispersalDistance", Database = "Bowman et al. 2002")%>%
  filter(!is.na(DispersalDistance))%>%
  unique()

## get rid of species from Sutherland 
bowman_dd <- filter(bowman_dd, !scientificName %in% suth_mamm_dd$scientificName)

## duplicates between Bowman and Santini
## exclude observations from Bowman with same DispersalDistance as Santini since has more info
## (not Code because it was inferred for Bowman)
key = sant_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(data_id = paste(scientificName, DispersalDistance, sep = "_")) %>%
  select(data_id)

bowman_dd <- bowman_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  filter(!paste(scientificName, DispersalDistance, sep = "_") %in% key$data_id)

## check how many species in bioshifts 
length(which(unique(bowman_dd$scientificName) %in% unique(sp$scientificName))) ## 0
bowman_sp <- unique(bowman_dd$scientificName)[which(unique(bowman_dd$scientificName) %in% unique(sp$scientificName))]


## Smith and Green 2005:
## read in data
sg = read_csv("data-raw/dispersal/Smith_and_Green_2005.csv") 
colnames(sg) <- str_replace_all(colnames(sg), "\\ ", "_")
length(unique(sg$Species)) #94 spp

sg_harm <- harmonize(sg$Species)

notfound <- filter(sg_harm, is.na(db_code))

## rename columns 
sg <- sg %>%
  rename("reported_name" = Species) %>%
  mutate(reported_name_fixed = reported_name)

## join with bioshifts 
sg <- left_join(sg, sg_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## reformat:
sg_dd <- sg %>%
  select(all_of(cols_to_keep), Reference, 
         `Max_distance_recorded_(m)`, Method) %>%
  rename("DispersalDistance" = `Max_distance_recorded_(m)`,
         "Source" = Reference, "ObservationType" = Method) %>%
  mutate(Sex = NA, Unit = "m", Field = "Max_distance_recorded_(m)",
         Database = "Smith & Green 2005", Code = "MaxDispersalDistance") %>%
  filter(!is.na(DispersalDistance)) %>%
  mutate(ObservationType = ifelse(str_detect(ObservationType, "MRR"),
                                  "mark-release-recapture",
                                  ifelse(str_detect(ObservationType, "RATE"),
                                         "yearly rate of movement from introduction",
                                         ifelse(str_detect(ObservationType, "UN"),
                                                "unknown",
                                                ifelse(str_detect(ObservationType, "WETLAND"),
                                                       "distance to nearest wetland",
                                                       ifelse(str_detect(ObservationType, "RAD"),
                                                              "radio-tracking",
                                                              ifelse(str_detect(ObservationType, "MARK"),
                                                                     "mark-recapture",
                                                                     NA))))))) 

## check how many species in bioshifts 
length(which(unique(sg_dd$scientificName) %in% unique(sp$scientificName))) ## 11
sg_sp <- unique(sg_dd$scientificName)[which(unique(sg_dd$scientificName) %in% unique(sp$scientificName))]


## Jenkins et al. 2007
jenkins <- read_csv("data-raw/dispersal/Jenkins_et_al_2007.csv")
colnames(jenkins) <- str_replace_all(colnames(jenkins), "\\ ", "_")
length(unique(jenkins$Scientific_Name)) #792 spp

## get rid of subspecies 
jenkins$genus = str_split_fixed(jenkins$Scientific_Name, " ", 2)[,1]
jenkins$species = str_split_fixed(jenkins$Scientific_Name, " ", 3)[,2]
jenkins$subspecies = str_split_fixed(jenkins$Scientific_Name, " ", 3)[,3]

jenkins = jenkins %>%
  mutate(genus = ifelse(genus == "", NA, genus)) %>%
  mutate(species = ifelse(species == "", NA, species)) %>%
  mutate(subspecies = ifelse(subspecies == "", NA, subspecies)) %>%
  mutate(genus_species = ifelse(is.na(species), genus,
                                paste(genus, species, sep = " ")))

## harmonize taxonomy in each 
jenkins_harm <- harmonize(jenkins$genus_species)

## fix some taxonomy
notfound <- filter(jenkins_harm, is.na(db_code))

## rename columns 
jenkins <- jenkins %>%
  rename("reported_name_fixed" = genus_species,
         "reported_name" = Scientific_Name)

jenkins <- left_join(jenkins, jenkins_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## Jenkins
jenkins_dd <- jenkins %>%
  select(all_of(cols_to_keep), 
         `Max._Indiv._Dispersal_Distance_(m)`) %>%
  rename("DispersalDistance" = `Max._Indiv._Dispersal_Distance_(m)`) %>%
  mutate(Sex = NA, ObservationType = "individual dispersal distance", Unit = "m", 
         Field = "Max._Indiv._Dispersal_Distance_(m)",
         Code = "MaxDispersalDistance", Database = "Jenkins et al. 2007",
         Source = NA) %>%
  filter(!is.na(DispersalDistance))%>%
  unique()

## Jenkins took from Smith and Green, but not everything
## exclude observations from Jenkins with same Code and DispersalDistance as those from Smith and Green
key = sg_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(data_id = paste(scientificName, Code, DispersalDistance, sep = "_")) %>%
  select(data_id)

jenkins_dd <- jenkins_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  filter(!paste(scientificName, Code, DispersalDistance, sep = "_") %in% key$data_id)

## check how many species in bioshifts 
length(which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))) ## 310
jenkins_sp <- unique(jenkins_dd$scientificName)[which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))]


## Flores et al. 2013:
## seed traps, tracked individual seeds, marked and recaptured seeds and estimated dispersal distances based on tracking vectors and calculating gut or fur retention times
flores <- read.delim("data-raw/dispersal/Flores_et_al_2013.txt") %>%
  select(-family, -order)
colnames(flores) <- str_replace_all(colnames(flores), "\\ ", "_")
length(unique(flores$Species)) #56 spp

flores_harm <- harmonize(flores$Species)

notfound <- filter(flores_harm, is.na(db_code))

## rename columns 
flores <- flores %>%
  rename("reported_name" = Species) %>%
  mutate(reported_name_fixed = reported_name)

flores <- left_join(flores, flores_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## Flores
flores_dd <- flores %>%
  select(all_of(cols_to_keep), Mean.dispersal.distance.m., Maximum.dispersal.distance.m.,
         Reference) %>%
  gather(key = "Field", value = "DispersalDistance", 
         c(Mean.dispersal.distance.m., Maximum.dispersal.distance.m.)) %>%
  mutate(Code = ifelse(Field == "Mean.dispersal.distance.m.", 
                       "MeanDispersalDistance", 
                       ifelse(Field == "Maximum.dispersal.distance.m.",
                              "MaxDispersalDistance", 
                              NA))) %>%
  mutate(Sex = NA, ObservationType = "seed/plant dispersal (field)", Unit = "m", 
         Database = "Flores et al. 2013") %>%
  filter(!is.na(DispersalDistance)) %>%
  rename("Source" = Reference) %>%
  mutate(Source = ifelse(str_detect(Source, "Bossard CC. The Role of Habitat Disturbance"),
                         "Bossard CC. The Role of Habitat Disturbance, Seed Predation and Ant Dispersal on Establishment of the Exotic Shrub Cytisus scoparius in California. American Midland Naturalist. 1991;126(1):1",
                         Source)) %>%
  distinct()

## some duplicates from Flores in Jenkins
## exclude observations from Jenkins with same Code and DispersalDistance as those from Flores since has more info
key = flores_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(data_id = paste(scientificName, Code, DispersalDistance, sep = "_")) %>%
  select(data_id)

jenkins_dd <- jenkins_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  filter(!paste(scientificName, Code, DispersalDistance, sep = "_") %in% key$data_id)

## check how many species in bioshifts 
length(which(unique(flores_dd$scientificName) %in% unique(sp$scientificName))) ## 37
flores_sp <- unique(flores_dd$scientificName)[which(unique(flores_dd$scientificName) %in% unique(sp$scientificName))]

## check how many species in bioshifts 
length(which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))) ## 308
jenkins_sp <- unique(jenkins_dd$scientificName)[which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))]


## Tamme et al 2014
tamme <- read_csv("data-raw/dispersal/Tamme_DispersalDistanceData.csv")
colnames(tamme) <- str_replace_all(colnames(tamme), "\\ ", "_")
length(unique(tamme$Species)) #576 spp

## get rid of subspecies 
tamme$genus = str_split_fixed(tamme$Species, " ", 3)[,1]
tamme$species = str_split_fixed(tamme$Species, " ", 3)[,2]
tamme$subspecies = str_split_fixed(tamme$Species, " ", 3)[,3]

tamme = tamme %>%
  mutate(genus = ifelse(genus == "", NA, genus)) %>%
  mutate(species = ifelse(species == "", NA, species)) %>%
  mutate(subspecies = ifelse(subspecies == "", NA, subspecies)) %>%
  mutate(genus_species = ifelse(is.na(species), genus,
                                paste(genus, species, sep = " ")))

## harmonize taxonomy in each 
tamme_harm <- harmonize(tamme$genus_species)

## fix some taxonomy
notfound <- filter(tamme_harm, is.na(db_code))
notfound$species[which(notfound$species == "Picea excelsa")] <- "Picea abies"
tamme_harm_notfound <- harmonize(notfound$species)

tamme$reported_name_fixed = tamme$genus_species
tamme$reported_name_fixed[which(tamme$reported_name_fixed == "Picea excelsa")] <- "Picea abies"
tamme$reported_name_fixed[which(tamme$reported_name_fixed == "Inula conyzae")] <- "Inula conyzae"

## bind:
tamme_harm <-filter(tamme_harm, !is.na(db_code)) %>%
  rbind(., tamme_harm_notfound)

tamme$reported_name_fixed[which(tamme$Species == "Picea excelsa")] <- "Picea abies"
tamme$reported_name_fixed[which(tamme$Species == "Inula conyzae")] <- "Inula conyza"

## rename columns 
tamme <- tamme %>%
  rename("reported_name" = Species)

tamme <- left_join(tamme, tamme_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()


## Tamme
tamme_dd <- tamme %>%
  select(all_of(cols_to_keep), 
         `Maximum_recorded_dispersal_distance_(m)`,
         `Mean_dispersal_distance_(m)`,
         `Median_dispersal_distance_(m)`, 
         `Mode_dispersal_distance_(m)`,
         `90th_percentile_dispersal_distance_(m)`,
         `99th_percentile_dispersal_distance_(m)`,
         Data_type, Reference) %>%
  gather(key = "Field", value = "DispersalDistance", c(`Maximum_recorded_dispersal_distance_(m)`,
                                                       `Mean_dispersal_distance_(m)`,
                                                       `Median_dispersal_distance_(m)`, 
                                                       `Mode_dispersal_distance_(m)`,
                                                       `90th_percentile_dispersal_distance_(m)`,
                                                       `99th_percentile_dispersal_distance_(m)`)) %>%
  mutate(Code = ifelse(str_detect(.$Field, "Mean"),
                       "MeanDispersalDistance", 
                       ifelse(str_detect(.$Field, "Maximum"),
                              "MaxDispersalDistance", 
                              ifelse(str_detect(.$Field, "Median"),
                                     "MedianDispersalDistance", 
                                     ifelse(str_detect(.$Field, "90"),
                                            "90thPercentileDispersalDistance", 
                                            ifelse(str_detect(.$Field, "99"),
                                                   "99thPercentileDispersalDistance",
                                                   ifelse(str_detect(.$Field, "Mode"),
                                                          "ModeDispersalDistance",
                                                          NA))))))) %>%
  mutate(Unit = "m", Database = "Tamme et al. 2014", Sex = NA) %>%
  filter(!is.na(DispersalDistance)) %>%
  rename("ObservationType" = Data_type, "Source" = Reference) %>%
  mutate(ObservationType = ifelse(ObservationType == "field", 
                                  "seed dispersal (field)", 
                                  ifelse(ObservationType == "model", 
                                         "seed dispersal (model)",
                                         NA)))

## Tamme took from Jenkins
## exclude observations from Jenkins with same Code and DispersalDistance as Tamme since has more info
key = tamme_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(data_id = paste(scientificName, Code, DispersalDistance, sep = "_")) %>%
  select(data_id)

jenkins_dd <- jenkins_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  filter(!paste(scientificName, Code, DispersalDistance, sep = "_") %in% key$data_id)

## also exclude ones with different Code from source Augspurger 1986
key = tamme_dd %>%
  filter(Source == "Augspurger 1986") %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(data_id = paste(scientificName, DispersalDistance, sep = "_")) %>%
  select(data_id)

jenkins_dd <- jenkins_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  filter(!paste(scientificName, DispersalDistance, sep = "_") %in% key$data_id)

## check how many species in bioshifts 
length(which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))) ## 224
jenkins_sp <- unique(jenkins_dd$scientificName)[which(unique(jenkins_dd$scientificName) %in% unique(sp$scientificName))]


## looks like Tamme and Flores have some overlapping data but rounded differently 
## if duplicated after rounding and source is from Stamp 1989 or Yumoto 1999, get rid of it in Tamme
key = flores_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(data_id = paste(scientificName, Code, round(DispersalDistance, 2), sep = "_")) %>%
  select(data_id)

tamme_dd <- tamme_dd %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  filter(!(paste(scientificName, Code, round(DispersalDistance, 2), sep = "_") %in% 
           key$data_id & Source %in% c("Stamp 1989", "Yumoto 1999")))

## check how many species in bioshifts 
length(which(unique(tamme_dd$scientificName) %in% unique(sp$scientificName))) ## 346
tamme_sp <- unique(tamme_dd$scientificName)[which(unique(tamme_dd$scientificName) %in% unique(sp$scientificName))]


## TRY
## read in TRY query
try = read_delim("data-raw/primary-trait-data/TRY/23454.txt")
unique(try$TraitName)

## get only dispersal data 
try_dd <- try %>%
  filter(TraitName == "Dispersal distance")
length(unique(try_dd$SpeciesName)) #152 spp

try_dd_harm <- harmonize(try_dd$SpeciesName)

notfound <- filter(try_dd_harm, is.na(db_code))

## rename columns 
try_dd <- try_dd %>%
  rename("reported_name" = SpeciesName) %>%
  mutate(reported_name_fixed = reported_name)

try_dd <- left_join(try_dd, try_dd_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## TRY
try_dd_dd <- try_dd %>%
  select(all_of(cols_to_keep), 
         OrigValueStr,
         OriglName, OrigUnitStr, Reference) %>%
  rename("Source" = Reference, "Field" = OriglName, "DispersalDistance" = OrigValueStr, "Unit" = OrigUnitStr) %>%
  mutate(Code = ifelse(Field %in% c("MeanDispersalDistanceMean", "MeanDispersalDistanceMax", 
                                    "MeanDispersalDistanceMin"), 
                       "MeanDispersalDistance", 
                       ifelse(Field %in%c("MaxDispersalDistanceMax", "Max Seed Dispersal Distance",
                                          "MaxDispersalDistanceMin", "MaxDispersalDistanceMean",
                                          "Seed_dispersal_distance_95",
                                          "Seed_dispersal_distance_5", "Max Seed Dispersal Distance"),
                              "MaxDispersalDistance", 
                              ifelse(Field %in% c("Effective Seed Dispersal Distance", "Seed Dispersal distance"),
                                     "DispersalDistance",
                                     NA)))) %>%
  mutate(Sex = NA, ObservationType = "unknown",  
         Database = "TRY database") %>% 
  filter(!is.na(Unit)) 

## check how many species in bioshifts 
length(which(unique(try_dd_dd$scientificName) %in% unique(sp$scientificName))) ## 55
try_dd_sp <- unique(try_dd_dd$scientificName)[which(unique(try_dd_dd$scientificName) %in% unique(sp$scientificName))]


## Trochet_et_al_2014:
## read in data
troc = read_csv("data-raw/dispersal/Trochet_et_al_2014_movement.csv") 
colnames(troc) <- str_replace_all(colnames(troc), "\\ ", "_")
length(unique(troc$Species)) #86 spp

troc_harm <- harmonize(troc$Species)

notfound <- filter(troc_harm, is.na(db_code))

## rename columns 
troc <- troc %>%
  rename("reported_name" = Species) %>%
  mutate(reported_name_fixed = reported_name)

troc <- left_join(troc, troc_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## reorganize
troc_dd <- troc %>%
  select(all_of(cols_to_keep), 
         Mean_dispersal_distance, Maximum_dispersal_distance, Minimal_dispersal_distance, 
         Mean_migration_distance, Maximum_migration_distance, Minimal_migration_distance) %>%
  gather(key = "Field", value = "DispersalDistance", 
         c(Mean_dispersal_distance, Maximum_dispersal_distance, Minimal_dispersal_distance, 
         Mean_migration_distance, Maximum_migration_distance, Minimal_migration_distance)) %>%
  filter(DispersalDistance != "DD") %>%
  mutate(Code = ifelse(Field %in% c("Mean_dispersal_distance", "Mean_migration_distance"), 
                       "MeanDispersalDistance", 
                       ifelse(Field %in%c("Maximum_dispersal_distance", "Maximum_migration_distance"),
                              "MaxDispersalDistance", 
                              ifelse(Field %in% c("Minimal_dispersal_distance", "Minimal_migration_distance"),
                                     "MinDispersalDistance",
                                     NA)))) %>%
  mutate(Sex = NA, Source = NA, Unit = "m",
         ObservationType = "mark-release-recapture and/or radio-tracking studies",
         Database = "Trochet et al. 2014") 

## check how many species in bioshifts 
length(which(unique(troc_dd$scientificName) %in% unique(sp$scientificName))) ## 17
troc_sp <- unique(troc_dd$scientificName)[which(unique(troc_dd$scientificName) %in% unique(sp$scientificName))]


## Stevens:
## Mammals - came from Whitmee & Orme - exclude 
## Birds - Paradis et al. 2002 - manually extracted
## Amphibians - Smith & Green 2005 - extracted from tables 3 and 4
## Spiders - Bonte et al. 2002,2003, Entling et al. 2011 and Pétillon et al. 2012 - not looked into yet
## Beetles - Turin 1999 - Turin, H. (1999) De Nederlandse Loopkevers. Uitgeverij KNNV & EIS, Nederland, 666p.- not looked into yet
## Butterflies - Stevens et al. 2010b - A meta-analysis of dispersal in butterflies. - not available in paper 

## to do: 
## go through not found species and see why they were not found (clean names more)
## check to make sure no data is replicated (eg. some data is the same between sg and jenkins, tamme and jenkins)
## clean reference column (remove brackets, standardize)


## check how many unique species we have data for:
species_with_dd <- append(co_sp, wo_sp) %>%
  append(., tamme_sp) %>%
  append(., shanks9_sp) %>%
  append(., shanks3_sp) %>%
  append(., suth_mamm_sp) %>%
  append(., suth_bird_sp) %>%
  append(., bowman_sp) %>%
  append(., jenkins_sp) %>%
  append(., flores_sp) %>%
  append(., try_dd_sp) %>%
  append(., sant_sp) %>%
  append(., troc_sp) %>%
  append(., par_sp) %>%
  append(., sg_sp)

length(unique(species_with_dd)) # 651 species 

## now: make subsets of each database with only bioshifts species 
co_sub = filter(co_dd, scientificName %in% co_sp)
wo_sub = filter(wo_dd, scientificName %in% wo_sp)
tamme_sub = filter(tamme_dd, scientificName %in% tamme_sp)
shanks9_sub = filter(shanks9_dd, scientificName %in% shanks9_sp)
shanks3_sub = filter(shanks3_dd, scientificName %in% shanks3_sp)
suth_mamm_sub = filter(suth_mamm_dd, scientificName %in% suth_mamm_sp)
suth_bird_sub = filter(suth_bird_dd, scientificName %in% suth_bird_sp)
bowman_sub = filter(bowman_dd, scientificName %in% bowman_sp)
jenkins_sub = filter(jenkins_dd, scientificName %in% jenkins_sp)
flores_sub = filter(flores_dd, scientificName %in% flores_sp)
try_dd_sub = filter(try_dd_dd, scientificName %in% try_dd_sp)
sant_sub = filter(sant_dd, scientificName %in% sant_sp)
troc_sub = filter(troc_dd, scientificName %in% troc_sp)
par_sub = filter(par_dd, scientificName %in% par_sp)
sg_sub = filter(sg_dd, scientificName %in% sg_sp)

## collate all bioshifts data
dd_collated <- rbind(co_sub, wo_sub) %>%
  rbind(., tamme_sub) %>%
  rbind(., shanks9_sub) %>%
  rbind(., shanks3_sub) %>%
  rbind(., suth_mamm_sub) %>%
  rbind(., suth_bird_sub) %>%
  rbind(., bowman_sub) %>%
  rbind(., jenkins_sub) %>%
  rbind(., flores_sub) %>%
  rbind(., try_dd_sub) %>%
  rbind(., sant_sub) %>%
  rbind(., troc_sub) %>%
  rbind(., par_sub) %>%
  rbind(., sg_sub) %>%
  unique()

length(unique(dd_collated$scientificName)) # 651


## fix classes 
sp <- read_csv("data-raw/splist.csv")
missing_class <- dd_collated %>%
  filter(is.na(class)) 

missing_class_sp <- filter(sp, scientificName %in% missing_class$scientificName) %>%
  select(-v1, -v2, -species)

missing_class <- select(missing_class, -class) %>%
  left_join(., missing_class_sp)

dd_collated <- filter(dd_collated, !scientificName %in% missing_class$scientificName) %>%
  rbind(., missing_class)

dd_collated$class[which(is.na(dd_collated$class))] <- "Thecostraca"

write.csv(dd_collated, "data-processed/dispersal-distance-collated.csv", row.names = FALSE)
dd_collated <- read.csv("data-processed/dispersal-distance-collated.csv")

## collate all data 
dd_collated_ALL <- rbind(co_dd, wo_dd) %>%
  rbind(., tamme_dd) %>%
  rbind(., shanks9_dd) %>%
  rbind(., shanks3_dd) %>%
  rbind(., suth_mamm_dd) %>%
  rbind(., suth_bird_dd) %>%
  rbind(., bowman_dd) %>%
  rbind(., jenkins_dd) %>%
  rbind(., flores_dd) %>%
  rbind(., try_dd_dd) %>%
  rbind(., sant_dd) %>%
  rbind(., troc_dd) %>%
  rbind(., par_dd) %>%
  rbind(., sg_dd) %>%
  unique()

length(unique(dd_collated_ALL$scientificName)) # 1444

write.csv(dd_collated_ALL, "data-processed/dispersal-distance-collated_ALL.csv", row.names = FALSE)
dd_collated_ALL <- read.csv("data-processed/dispersal-distance-collated_ALL.csv")


## make some plots 
dd_collated %>%
  select(scientificName, class, kingdom) %>%
  unique() %>% 
  group_by(class) %>%
  tally() %>%
  left_join(., select(dd_collated, class, kingdom)) %>%
  group_by(kingdom) %>%
  unique() %>%
  ungroup() %>%
  mutate(class = factor(class, levels = as.character(class), ordered = TRUE)) %>%
  ggplot(., aes(x = class, y = n, fill = kingdom)) + geom_col() + 
  coord_flip() +
  theme_minimal() + 
  labs(x = "", y = "Species with known dispersal distance")

left_join(sp, dd_collated) %>%
  mutate(has_dd = ifelse(!is.na(DispersalDistance), 
                          "Yes", "No")) %>%
  select(scientificName, class, has_dd) %>%
  unique() %>%
  ggplot(aes(x = class, fill = has_dd)) + geom_bar() + 
  theme_minimal() +
  coord_flip()

dd_collated %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  ggplot(aes(x = log(DispersalDistance))) + geom_histogram() + 
  facet_wrap(~Code)


### collate dispersal frequency/longevity data
## longevity: 
lifespan <- read.csv("data-raw/CompilationLifeSpan_02242022.csv")
unique(lifespan$Database)

## filter to species with dispersal distance
lifespan <- filter(lifespan, SpeciesChecked %in% dd_collated$scientificName)
length(unique(lifespan$SpeciesChecked)) ## 490 of our species 

## how many are animals?
lifespan = filter(lifespan, phylum == "Chordata")
length(unique(lifespan$SpeciesChecked)) ## 199 of our species 

dd_collated %>%
  select(kingdom, scientificName) %>%
  unique() %>%
  summarise(length(which(kingdom == "Animalia"))) # out of the 234 animals - nice!!! 

unique(lifespan$class)

## look at missing ones
missing_long = dd_collated %>%
  filter(kingdom == "Animalia" & !scientificName %in% lifespan$SpeciesChecked) %>%
  select(scientificName, family, order, class, phylum, kingdom) %>%
  unique()

## see if we can get these from anage after harmonizing
anage <- read.delim("data-raw/primary-trait-data/AnAge/anage_data.txt") %>%
  mutate(SpeciesName = paste(Genus, Species, sep = " "))

## harmonize names 
anage_harm <- harmonize(anage$SpeciesName)

notfound <- filter(anage_harm, is.na(db_code))

## rename columns 
anage <- anage %>%
  rename("reported_name" = SpeciesName) %>%
  mutate(reported_name_fixed = reported_name)

anage <- left_join(anage, anage_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## search it for missing species
anage_bs <- filter(anage, scientificName %in% missing_long$scientificName) %>%
  filter(!is.na(Maximum.longevity..yrs.))
## 14 more species!!!! 

## reformat data and add to lifespan 
colnames(lifespan)
anage_long <- anage_bs %>%
  select(scientificName, reported_name, phylum, class, order, family, Genus, Species,
         Maximum.longevity..yrs., Source) %>%
  mutate(Unit = "yrs", Database = "anAge", Field = "Maximum.longevity..yrs.") %>%
  rename("LifeSpan" = Maximum.longevity..yrs., "SpeciesChecked" = scientificName,
         "Taxa" = reported_name, "genus" = Genus, "species" = Species)

lifespan = lifespan[,which(colnames(lifespan) %in% colnames(anage_long))]
anage_long = anage_long[,which(colnames(anage_long) %in% colnames(lifespan))]

lifespan = rbind(lifespan, anage_long)

## write
write.csv(lifespan, "data-processed/longevity.csv", row.names = FALSE)

## look at missing ones
missing_long = dd_collated %>%
  filter(kingdom == "Animalia" & !scientificName %in% lifespan$SpeciesChecked) %>%
  select(scientificName, family, order, class, phylum, kingdom) %>%
  unique()


## for plants, want to know age at first maturity 
## TRY 
try_am = read_delim("data-raw/primary-trait-data/TRY/23659.txt")
unique(try_am$TraitName)

## age at maturity 
try_am <- try_am %>%
  filter(TraitName == "Plant ontogeny: age of maturity (first flowering)") 

## harmonize names 
try_am_harm <- harmonize(try_am$SpeciesName)

notfound <- filter(try_am_harm, is.na(db_code))

## rename columns 
try_am <- try_am %>%
  rename("reported_name" = SpeciesName) %>%
  mutate(reported_name_fixed = reported_name)

try_am <- left_join(try_am, try_am_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## subset to bioshifts species with dispersal distance
try_am_bs <- filter(try_am, scientificName %in% dd_collated$scientificName)

length(unique(try_am_bs$scientificName)) # 269 species 

dd_collated %>%
  select(kingdom, scientificName) %>%
  unique() %>%
  summarise(length(which(kingdom == "Plantae"))) # out of the 415 plants - nice!!! 

## clean the data 
unique(try_am_bs$OriglName)

try_am_bs = try_am_bs %>%
  filter(!OriglName %in% c("flowering", "Secondary juvenile period")) %>%
  mutate(Unit = ifelse(str_detect(OrigValueStr, "yrs") | OrigUnitStr %in% c("years?", "yr", "years"), 
                        "yrs",
                       ifelse(str_detect(OrigUnitStr, "day"), 
                              "days", 
                              OrigUnitStr)))
try_am_bs$Unit[which(str_detect(try_am_bs$OrigValueStr, "weeks"))] = "weeks"

unique(try_am_bs$OrigUnitStr)
unique(try_am_bs$Unit)

try_am_bs$Unit[which(is.na(try_am_bs$Unit))] <- "yrs"

## want minimum age at maturity to give a maximum dispersal estimate 
try_am_bs = try_am_bs %>%
  mutate(AgeAtMaturity = str_replace_all(OrigValueStr,"\\<", ""),
         AgeAtMaturity = str_replace_all(AgeAtMaturity,"\\>", ""),
         AgeAtMaturity = str_split_fixed(AgeAtMaturity, "\\-", 2)[,1]) %>%
  mutate(AgeAtMaturity = str_replace_all(AgeAtMaturity, "Over ", "")) %>%
  mutate(AgeAtMaturity = str_replace_all(AgeAtMaturity, "Within ", "")) %>%
  mutate(AgeAtMaturity = str_replace_all(AgeAtMaturity, "Between ", "")) %>%
  mutate(AgeAtMaturity = str_split_fixed(AgeAtMaturity, " ", 2)[,1])

unique(try_am_bs$AgeAtMaturity)

try_am_bs <- try_am_bs %>%
  select(all_of(cols_to_keep), 
         Unit, AgeAtMaturity, OriglName, 
         Reference) %>%
  mutate(Database = "TRY") %>%
  filter(!is.na(AgeAtMaturity)) %>%
  rename("Source" = Reference, "Field" = OriglName)

## save 
write.csv(try_am_bs, "data-processed/age-at-maturity-TRY.csv", row.names = FALSE)


## fill in the gaps:
missing_am = dd_collated %>%
  filter(kingdom == "Plantae" & !scientificName %in% try_am_bs$scientificName) %>%
  select(scientificName, family, order, class, phylum, kingdom) %>%
  unique()


## using lifespan/growth form to infer age at maturity for annual plants
## for species that live 1 year, age at maturity = 1 year 
## longevity: 
lifespan <- read.csv("data-raw/CompilationLifeSpan_02242022.csv")
unique(lifespan$Database)

## filter to species with dispersal distance
lifespan <- filter(lifespan, SpeciesChecked %in% dd_collated$scientificName)
length(unique(lifespan$SpeciesChecked)) ## 490 of our species 

## how many are plants?
lifespan = filter(lifespan, phylum %in% c("Streptophyta", "Tracheophyta"))
length(unique(lifespan$SpeciesChecked)) ## 291 of our species 

## harmonize names 
lifespan_harm <- harmonize(lifespan$SpeciesChecked)

notfound <- filter(lifespan_harm, is.na(db_code))

## rename columns 
lifespan <- lifespan %>%
  rename("reported_name" = SpeciesChecked) %>%
  mutate(reported_name_fixed = reported_name) %>%
  select(-family, -genus, -phylum, -class, -order)

lifespan <- left_join(lifespan, lifespan_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## clean lifespan
unique(lifespan$LifeSpan)
unique(lifespan$LifeSpan)[str_detect(unique(lifespan$LifeSpan), "1\\-")]
unique(lifespan$LifeSpan)[str_detect(unique(lifespan$LifeSpan), "1\\,")]

lifespan_am <- lifespan %>%
  mutate(Unit = ifelse(LifeSpan %in% c("annuals", "summer annuals",
                                       "winter annuals"),
                       "years",
                       Unit)) %>%
  mutate(AgeAtMaturity = ifelse(LifeSpan %in% c("annuals", "summer annuals",
                                                "winter annuals"),
                                1,
                                ifelse(str_detect(LifeSpan, "1\\-"), 
                                       1,
                                       ifelse(str_detect(LifeSpan, "1\\,"), 
                                              1,
                                              ifelse(str_detect(LifeSpan, "\\<1"), 
                                                     1,
                                                     ifelse(LifeSpan == 1,
                                                            1, 
                                                            NA)))))) %>%
  filter(AgeAtMaturity == 1) %>% 
  select(all_of(cols_to_keep), -X, -Taxonomic.Group, -Ecosystem.Type, -Group, -ecotype,
         AgeAtMaturity, Unit, Database, Field) %>%
  unique() %>%
  filter(scientificName %in% missing_am$scientificName) %>%
  mutate(Source = "Lifespan compilation")

length(unique(lifespan_am$scientificName)) #8 more species 
  
## BROT
brot = read.csv("data-raw/primary-trait-data/BROT/BROT2_dat.csv")
unique(brot$Trait)
## GrowthForm, DispMode
brot = brot %>%
  filter(Trait %in% c("GrowthForm"))

## harmonize names 
brot_harm <- harmonize(brot$Taxon)

notfound <- filter(brot_harm, is.na(db_code))

## rename columns 
brot <- brot %>%
  rename("reported_name" = Taxon) %>%
  mutate(reported_name_fixed = reported_name)

brot <- left_join(brot, brot_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## subset to bioshifts plant species missing data on age at maturity 
brot_sub <- filter(brot, scientificName %in% missing_am$scientificName) %>%
  filter(!scientificName %in% lifespan_am$scientificName)

## use information on whether each species is annual or perennial 
unique(brot_sub$Data)

brot_am <- brot_sub %>%
  select(all_of(cols_to_keep), Data, SourceID) %>%
  mutate(AgeAtMaturity = ifelse(str_detect(Data, "annual"),
                                "1", 
                                NA)) %>%
    rename("Source" = SourceID) %>%
  mutate(Field = "Growth form",
         Unit = "yrs", Database = "TRY") %>%
  filter(!is.na(AgeAtMaturity)) %>%
  select(-Data)

length(unique(brot_am$scientificName))# 1 more

## LEDA
## Growth form
leda = read.delim("data-raw/primary-trait-data/LEDA/plant_growth_form copy.txt", sep = ";")
unique(leda$plant.growth.form)

# https://hosho.ees.hokudai.ac.jp/tsuyu/top/dct/lf.html
# plant life form classifications:
# Therophytes - annual herb

## harmonize names 
leda_harm <- harmonize(leda$SBS.name)

notfound <- filter(leda_harm, is.na(db_code))

## rename columns 
leda <- leda %>%
  rename("reported_name" = SBS.name) %>%
  mutate(reported_name_fixed = reported_name)

leda <- left_join(leda, leda_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## subset to bioshifts plant species missing data on age at maturity 
leda_sub <- filter(leda, scientificName %in% missing_am$scientificName) %>%
  filter(!scientificName %in% brot_am$scientificName) %>%
  filter(!scientificName %in% lifespan_am$scientificName)

## use information on whether each species is annual or perennial 
unique(leda_sub$plant.growth.form)

leda_am <- leda_sub %>%
  select(all_of(cols_to_keep), plant.growth.form, original.reference) %>%
  mutate(AgeAtMaturity = ifelse(str_detect(plant.growth.form, "Therophyte"),
                                "1", 
                                NA)) %>%
  rename("Source" = original.reference) %>%
  mutate(Field = "plant.growth.form",
         Unit = "yrs", Database = "BROT") %>%
  filter(!is.na(AgeAtMaturity)) %>%
  select(-plant.growth.form) 

length(unique(leda_am$scientificName)) # 0 more 


## TRY
## Growth form 
## read in TRY query
try_gf_full = read_delim("data-raw/primary-trait-data/TRY/23454.txt")
unique(try_gf_full$TraitName)

try_gf_full <- try_gf_full %>%
  filter(TraitName == "Plant growth form")

## do a quick filter of data we don't want 
## otherwise there are too many species to harmonize 
unique(try_gf_full$OrigValueStr)
unique(try_gf_full$OriglName)

try_gf <- try_gf_full %>%
  filter(str_detect(OrigValueStr, "annual") | str_detect(OrigValueStr, "Annual") | 
           str_detect(OrigValueStr, "therophyte"))

## harmonize names 
try_gf_harm <- harmonize(try_gf$SpeciesName)

notfound <- filter(try_gf_harm, is.na(db_code))

## rename columns 
try_gf <- try_gf %>%
  rename("reported_name" = SpeciesName) %>%
  mutate(reported_name_fixed = reported_name)

try_gf <- left_join(try_gf, try_gf_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## subset to bioshifts plant species missing data on age at maturity 
try_gf_sub <- filter(try_gf, scientificName %in% missing_am$scientificName) %>%
  filter(!scientificName %in% brot_am$scientificName) %>%
  filter(!scientificName %in% leda_am$scientificName) %>%
  filter(!scientificName %in% lifespan_am$scientificName)

## use information on whether each species is annual or perennial 
unique(try_gf_sub$OrigValueStr)
unique(try_gf_sub$OriglName)

try_gf_am <- try_gf_sub %>%
  select(all_of(cols_to_keep), 
         OrigValueStr, OriglName, 
         Reference) %>%
  mutate(Database = "TRY", Unit = "yrs") %>%
  mutate(AgeAtMaturity = ifelse(str_detect(OrigValueStr, "Annual"),
                                "1", 
                                NA)) %>%
  filter(!is.na(AgeAtMaturity)) %>%
  rename("Source" = Reference, "Field" = OriglName) %>%
  select(-OrigValueStr) %>%
  unique()

length(unique(try_gf_am$scientificName))# 0 species

## combine
all_am <- rbind(try_am_bs, lifespan_am) %>%
  rbind(., brot_am) %>%
  rbind(., leda_am) %>%
  rbind(., try_gf_am)

length(unique(all_am$scientificName)) ## 276 spp

## write 
write.csv(all_am, "data-processed/age-at-maturity.csv", row.names = FALSE)


