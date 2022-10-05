## sorting the trait databases by taxonomic group 
## and looking at the traits studied within each taxonomic group

## goal: gather all body size data for one taxonomic group in bioshifts 
## and make list of species with missing data for Clare to search for
## also give them a list of synonyms for each species 

## Brunno says plants and birds have most data - maybe focus on those?

library(tidyverse)
library(broom)
library(taxadb)
library(readxl)
library(readr)

theme_set(theme_bw())

## write function to harmonize taxonomy 
harmonize_taxnomy <- function(names) {
  
  ## first, search for GBIF accepted names
  ## GBIF
  gbif = filter_name(names, "gbif") %>%
    filter(taxonomicStatus == "accepted") %>%
    mutate(db = "gbif") %>%
    rename("db_code" = acceptedNameUsageID) %>%
    select(input, scientificName, db, db_code, kingdom, order, class, family, genus)
  
  not_found <- names[which(!names %in% gbif$input)]
  
  ## if not found, search itis
  ## itis
  itis = filter_name(not_found, "itis") %>%
    filter(taxonomicStatus == "accepted") %>%
    mutate(db = "itis") %>%
    rename("db_code" = acceptedNameUsageID) %>%
    select(input, scientificName, db, db_code, kingdom, order, class, family, genus)
  
  not_found <- not_found[which(!not_found %in% itis$input)]
  
  ## if not found, search col
  ## col
  col = filter_name(not_found, "col") %>%
    filter(taxonomicStatus == "accepted") %>%
    mutate(db = "col") %>%
    rename("db_code" = acceptedNameUsageID) %>%
    select(input, scientificName, db, db_code, kingdom, order, class, family, genus)
  
  not_found <- not_found[which(!not_found %in% col$input)]
  
  ## if not found, search ncbi
  ## ncbi
  ncbi = filter_name(not_found, "ncbi") %>%
    filter(taxonomicStatus == "accepted") %>%
    mutate(db = "ncbi") %>%
    rename("db_code" = acceptedNameUsageID) %>%
    select(input, scientificName, db, db_code, kingdom, order, class, family, genus)
  
  not_found <- not_found[which(!not_found %in% ncbi$input)]
  
  names_db <- rbind(itis, ncbi) %>% rbind(., col) %>% rbind(., gbif) %>%
    unique(.) %>%
    rename("genus_species" = input) %>%
    select(-scientificName)
  
  ## return database of names and list of species not found
  
  return(list(names_db, not_found))
}

## read in the list of databases:
dbs <- read_csv("data-raw/TraitDatabases.csv")
colnames(dbs) <- str_replace_all(colnames(dbs), "\\ ", "_")

## look at coverage
dbs %>%
  ggplot(., aes(x = plant_animal_fungi_both)) +
  geom_bar() +
  coord_flip() +
  scale_y_continuous(name = "Number of databases") +
  scale_x_discrete(name = "")

ggsave("figures/database-coverage/plant-animal-both.png", height = 2, width = 3)

dbs %>%
  ggplot(., aes(x = taxonomic_group_general)) +
  geom_bar() +
  coord_flip() +
  scale_y_continuous(name = "Number of databases") +
  scale_x_discrete(name = "")

ggsave("figures/database-coverage/general-group.png", height = 2, width = 6)

dbs %>%
  ggplot(., aes(x = taxonomic_group_specific)) +
  geom_bar() +
  coord_flip() +
  scale_y_continuous(name = "Number of databases") +
  scale_x_discrete(name = "")

ggsave("figures/database-coverage/general-specfic-unrefined.png", height = 4, width = 6)

## now get breakdown by specific taxonomic groups
specific_groups <- dbs %>%
  select(Database_name, taxonomic_group_specific) 

split <- str_split(specific_groups$taxonomic_group_specific, "\\, ")

df <- data.frame(do.call(rbind, split))
df$Database_name <- specific_groups$Database_name

df <- gather(df, key = "old_col", value = "taxonomic_group_specific", c(1:6)) %>%
  select(-old_col) %>%
  unique(.)

df <- left_join(df, select(dbs, -taxonomic_group_specific))

by_group <- df %>%
  mutate(group = ifelse(taxonomic_group_specific %in% c("squamates", "tetrapods", "reptiles", "lizards", "vertebrates",
                                                        "all"),
                                                        "Reptiles", 
                                                        NA)) %>%
  mutate(group = ifelse(taxonomic_group_specific %in% c("amphibians", "vertebrates", "all", "tetrapods"),
                        ifelse(is.na(group), 
                               "Amphibians", 
                               paste(group, "Amphibians", sep = "_")),
                        group))  %>%
  mutate(group = ifelse(taxonomic_group_specific %in% c("spiders", "insects", "copepods", "benthic invertebrates", 
                                                        "lotic invertebrates", "plankton", "beetles", "arthropods",
                                                        "ants", "aquatic invertebrates", "invertebrates", "all"),
                        ifelse(is.na(group), 
                               "Arthropods",
                               paste(group, "Arthropods", sep = "_")),
                        group)) %>%
  mutate(group = ifelse(taxonomic_group_specific %in% c("plants", "phytoplankton", "macroalgae", "plankton", "all"),
                        ifelse(is.na(group), 
                               "Aquatic and terrestrial plants",
                               paste(group, "Aquatic and terrestrial plants", sep = "_")),
                        group)) %>%
  mutate(group = ifelse(taxonomic_group_specific %in% c("fungi", "all"),
                        ifelse(is.na(group), 
                               "Fungi", 
                               paste(group, "Fungi", sep = "_")),
                        group)) %>%
  mutate(group = ifelse(taxonomic_group_specific %in% c("fish", "all", "vertebrates"),
                        ifelse(is.na(group), 
                               "Fish",
                               paste(group, "Fish", sep = "_")),
                        group)) %>%
  mutate(group = ifelse(taxonomic_group_specific %in% c("mammals", "all", "vertebrates","tetrapods"),
                        ifelse(is.na(group), 
                               "Mammals",
                               paste(group, "Mammals", sep = "_")),
                        group)) %>%
  mutate(group = ifelse(taxonomic_group_specific %in% c("birds", "all", "vertebrates"),
                        ifelse(is.na(group), 
                               "Birds", 
                               paste(group, "Birds", sep = "_")),
                        group)) %>%
  mutate(group = ifelse(taxonomic_group_specific %in% c("molluscs", "all", "invertebrates", 
                                                        "benthic invertebrates",
                                                        "lotic invertebrates"),
                        ifelse(is.na(group), 
                               "Molluscs",
                               paste(group, "Molluscs", sep = "_")),
                        group)) %>%
  mutate(group = ifelse(taxonomic_group_specific %in% c("annelids", "all", "invertebrates", "non-chordates", 
                                                        "worms", "corals"),
                        ifelse(is.na(group), 
                               "Annelids, non-chordates, and other worms",
                               paste(group, "Annelids, non-chordates, and other worms", sep = "_")),
                        group)) 

by_group <- by_group %>%
  select(Database_name, group) 

split <- str_split(by_group$group, "\\_")

df_by_group <- data.frame(do.call(rbind, split))
df_by_group$Database_name <- by_group$Database_name

df_by_group <- gather(df_by_group, key = "old_col", value = "group", c(1:10)) %>%
  select(-old_col) %>%
  unique(.)

df_by_group %>%
  ggplot(., aes(x = group)) +
  geom_bar() +
  coord_flip() +
  scale_y_continuous(name = "Number of databases") +
  scale_x_discrete(name = "") 

ggsave("figures/database-coverage/general-group-refined.png", height = 4, width = 6)

## write out table of each group + their possible databases and possible traits 
data_key <- df_by_group %>%
  left_join(., dbs) %>%
  group_by(group) %>%
  summarize(possible_databases = paste(Database_name, collapse = ", "),
            possible_traits = paste(Trait, collapse = ", "))

write.csv(data_key, "data-processed/trait-databases-by-group.csv", row.names = FALSE)


## read in the list of traits in each database:
traits <- read_csv("data-raw/DetailTraits.csv")

## break into one dataframe for each database 
trait_list <- t(traits) %>%
  split(., seq(nrow(.)))

names(trait_list) <- colnames(traits)

## subset to ONLY databases for birds 
names(trait_list)
bird_dbs <- str_split(data_key$possible_databases[which(data_key$group == "Birds")], "\\, ")[[1]]

## make sure all are in the list of traits
length(which(!bird_dbs %in% names(trait_list)))
bird_dbs[which(!bird_dbs %in% names(trait_list))]
## they are

## next: see what traits we have for birds
birds <- trait_list[names(trait_list) %in% bird_dbs]

birds <- data.frame(do.call(rbind, birds))
birds$Database_name = rownames(birds)

## make data base with columns:
## group 
## trait
## Database_name 
## trait_common_name
## trait_description
## useful_y_n
birds <- gather(birds, key = "old_col", value = "trait", c(1:2091)) %>%
  select(-old_col) %>%
  filter(!is.na(trait)) %>%
  unique(.) %>%
  arrange(Database_name) %>%
  mutate(trait_common_name = NA,
         useful_y_n = NA) 

## 298 traits 

## write out and start to decide on common names
write.csv(birds, "data-processed/traits-by-group/birds.csv", row.names = FALSE)


## Notes 

## AVONET
## Sex and Age columns indicate what the age and sex of specimen was 

## BioTIME_BodySize-0.1.1
## mishmash of other datasets - look at source scripts before adding in case it duplicates data

## EltonTraits
## trophic_niche options must be harmonized 
## variables listed not real variable names 
## species level body mass and other body mass available 

## GARD [appendix_1]
## maximum body length coverted to mass using clade-specific alometric equations

## GARD [novosolov_et_al._2017]
## body mass not compiled from the primary literature, but from other data (Feldman, Sabath, Pyron, Mayrose, Meiri, BirdLife International, Jones et al.)
## trophic_level - no scavenger category (check why)

## Heinen et al. 2017
## all species frugivores - decode diet category somehow 

## Lislevand_et_al_2007
## sample size for mass measurements

## Marine Species traits
## length measurement varies - see other column for age/sex of individual
## feeding type needs to be harmonized between trophic_level and trophic_niche
## is life span categorical? also check life stage 

## read in list of bioshifts species 
splist <- read.csv("data-raw/splist.csv")

## read in bird data 
birds_filled <- read.csv("data-processed/traits-by-group/birds_filled.csv") 

## get rid of traits that aren't common 
birds_filled <- birds_filled %>%
  filter(!is.na(trait_common_name))

unique(birds_filled$trait_common_name)

## read in all of the bird trait databases 
unique(birds_filled$Database_name)

#############
## amniota ##
#############
amniota <- read.csv("data-raw/primary-trait-data/amniota/ECOL_96_269/Data_Files/Amniote_Database_Aug_2015.csv")

## change -999 to NA
amniota <- data.frame(lapply(amniota, function(x) {gsub("-999", NA, x)})) %>%
  mutate(genus_species = paste(genus, species, sep = " "))

## taxonomic harmonization: change name to accepted name 
amniota_tax <- harmonize_taxnomy(amniota$genus_species)

amniota_accsp <- amniota_tax[[1]]
amniota_notf <- amniota_tax[[2]]

amniota_tax <- filter(amniota, genus_species %in% amniota_notf) %>%
  mutate(db = NA, db_code = NA, kingdom = "Animalia") %>%
  select(genus_species, db, db_code, kingdom, order, class, family, genus) %>%
  rbind(., amniota_accsp)

## merge with amniota
amniota <- left_join(amniota, amniota_tax)

## subset to bioshifts species 
amniota <- filter(amniota, genus_species %in% splist$species)
## from 21322 species -> 2137 species 

## save subset 
write.csv(amniota, "data-processed/primary-traits-bioshifts/amniota.csv", row.names = FALSE)


#############
##  AnAge  ##
#############
anage <- read_delim("data-raw/primary-trait-data/AnAge/anage_data.txt") %>%
  mutate(genus_species = paste(Genus, Species, sep = " ")) %>%
  rename("kingdom" = Kingdom, "order" = Order, "class" = Class, "family" = Family,
         "genus" = Genus, "species" = Species)

## taxonomic harmonization: change name to accepted name 
anage_tax <- harmonize_taxnomy(anage$genus_species)

anage_accsp <- anage_tax[[1]]
anage_notf <- anage_tax[[2]]

anage_tax <- filter(anage, genus_species %in% anage_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, kingdom, order, class, family, genus) %>%
  rbind(., anage_accsp)

## merge with anage
anage <- left_join(anage, anage_tax) 

## subset to bioshifts species 
anage <- filter(anage, genus_species %in% splist$species)
## from 4219 species -> 1042 species 

## save subset 
write.csv(anage, "data-processed/primary-traits-bioshifts/AnAge.csv", row.names = FALSE)


##############
##  AVONET  ##
##############
avonet <- read_csv("data-raw/primary-trait-data/AVONET/AVONET_Raw.csv") %>%
  mutate(genus_species = Species1, kingdom = "Animalia", class = "Aves",
         genus = str_split_fixed(Species1, " ", 2)[,1]) %>%
  rename("order" = Order1, "family" = Family1) %>%
  select(-Species1)

## taxonomic harmonization: change name to accepted name 
avonet_tax <- harmonize_taxnomy(avonet$genus_species)

avonet_accsp <- avonet_tax[[1]]
avonet_notf <- avonet_tax[[2]]

avonet_tax <- filter(avonet, genus_species %in% avonet_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, kingdom, order, class, family, genus) %>%
  rbind(., avonet_accsp)

## merge with avonet
avonet <- left_join(avonet, avonet_tax) 

## subset to bioshifts species 
avonet <- filter(avonet, genus_species %in% splist$species)
## from 1109 species -> 1912 species 

## save subset 
write.csv(avonet, "data-processed/primary-traits-bioshifts/AVONET.csv", row.names = FALSE)

###################
##  EltonTraits  ##
###################
etraits <- read_delim("data-raw/primary-trait-data/EltonTraits/BirdFuncDat.txt") %>%
  rename("genus_species" = Scientific, "order" = IOCOrder, "family" = BLFamilyLatin) %>%
  mutate(kingdom = "Animalia", class = "Aves",
         genus = str_split_fixed(genus_species, " ", 2)[,1]) %>%
  filter(!is.na(genus_species))

## taxonomic harmonization: change name to accepted name 
etraits_tax <- harmonize_taxnomy(etraits$genus_species)

etraits_accsp <- etraits_tax[[1]]
etraits_notf <- etraits_tax[[2]]

etraits_tax <- filter(etraits, genus_species %in% etraits_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, kingdom, order, class, family, genus) %>%
  rbind(., etraits_accsp)

## merge with etraits
etraits <- left_join(etraits, etraits_tax) 

## subset to bioshifts species 
etraits <- filter(etraits, genus_species %in% splist$species)
## from 9993 species -> 1807 species 

## save subset 
write.csv(etraits, "data-processed/primary-traits-bioshifts/EltonTraits.csv", row.names = FALSE)

#######################
##  GARD appendix 1  ##
#######################
gardapp1 <- read_excel("data-raw/primary-trait-data/GARD/appendix_1_-_data.xlsx") %>%
  rename("genus_species" = species, "class" = Class) %>%
  mutate(genus = str_split_fixed(genus_species, " ", 2)[,1],
         kingdom = "Animalia")

## taxonomic harmonization: change name to accepted name 
gardapp1_tax <- harmonize_taxnomy(gardapp1$genus_species)

gardapp1_accsp <- gardapp1_tax[[1]]
gardapp1_notf <- gardapp1_tax[[2]]

gardapp1_tax <- filter(gardapp1, genus_species %in% gardapp1_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, class, genus, kingdom) %>%
  mutate(order = NA, family = NA) %>%
  rbind(., gardapp1_accsp)

## merge with gardapp1
gardapp1 <- left_join(gardapp1, gardapp1_tax) %>%
  unique(.) ## get rid of duplicated species 

## subset to bioshifts species 
gardapp1 <- filter(gardapp1, genus_species %in% splist$species)
## from 3970 species -> 543 species 

## save subset 
write.csv(gardapp1, "data-processed/primary-traits-bioshifts/GARD_appendix_1.csv", row.names = FALSE)
 
####################################
##  GARD [novosolov_et_al._2017]  ##
####################################
gardnovo <- read_excel("data-raw/primary-trait-data/GARD/novosolov_et_al._2017_appendix_1_population_densities___range_sizes.xlsx") %>%
  rename("genus_species" = Binomial) %>%
  mutate(genus = str_split_fixed(genus_species, " ", 2)[,1],
         kingdom = "Animalia")

## taxonomic harmonization: change name to accepted name 
gardnovo_tax <- harmonize_taxnomy(gardnovo$genus_species)

gardnovo_accsp <- gardnovo_tax[[1]]
gardnovo_notf <- gardnovo_tax[[2]]

gardnovo_tax <- filter(gardnovo, genus_species %in% gardnovo_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, genus, kingdom) %>%
  mutate(order = NA, family = NA, class = NA) %>%
  rbind(., gardnovo_accsp)

## merge with gardnovo
gardnovo <- left_join(gardnovo, gardnovo_tax) 

## subset to bioshifts species 
gardnovo <- filter(gardnovo, genus_species %in% splist$species)
## from 1435 species -> 265 species 

## save subset 
write.csv(gardnovo, "data-processed/primary-traits-bioshifts/GARD_novosolov_et_al_2017.csv", row.names = FALSE)

##########################
##  Heinen et al. 2017  ##
##########################
heinen <- read_delim("data-raw/primary-trait-data/Heinen/doi_10.5061_dryad.s522m__v1/Data_Traits_IslandFrugivores.txt") %>%
  rename("family" = Family, "class" = Class, "order" = Order, "genus_species" = Species) %>%
  mutate(kingdom = "Animalia",
         genus = str_split_fixed(genus_species, " ", 2)[,1]) 

## taxonomic harmonization: change name to accepted name 
heinen_tax <- harmonize_taxnomy(heinen$genus_species)

heinen_accsp <- heinen_tax[[1]]
heinen_notf <- heinen_tax[[2]]

heinen_tax <- filter(heinen, genus_species %in% heinen_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, genus, kingdom, family, order, class) %>%
  rbind(., heinen_accsp)

## merge with heinen
heinen <- left_join(heinen, heinen_tax) 

## subset to bioshifts species 
heinen <- filter(heinen, genus_species %in% splist$species)
## from 387 species -> 46 species 

## save subset 
write.csv(heinen, "data-processed/primary-traits-bioshifts/Heinen_et_al_2017.csv", row.names = FALSE)


############################
##  Lislevand_et_al_2007  ##
############################
lislevand <- read_excel("data-raw/primary-trait-data/Lislevand_et_al_2007/Lislevand2007.xlsx") %>%
  rename("genus_species" = Species_name) %>%
  mutate(kingdom = "Animalia",
         genus = str_split_fixed(genus_species, " ", 2)[,1]) 

## change -999 to NA
lislevand <- data.frame(lapply(lislevand, function(x) {gsub("-999", NA, x)}))
  
## read in family key and create column with family name 
lislevand <-  read_csv("data-raw/primary-trait-data/Lislevand_et_al_2007/family_key.csv") %>%
  mutate(Family = as.character(Family)) %>%
  left_join(lislevand, .) %>%
  select(-Family)

## taxonomic harmonization: change name to accepted name 
lislevand_tax <- harmonize_taxnomy(lislevand$genus_species)

lislevand_accsp <- lislevand_tax[[1]]
lislevand_notf <- lislevand_tax[[2]]

lislevand_tax <- filter(lislevand, genus_species %in% lislevand_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, genus, kingdom, family) %>%
  mutate(order = NA, class = NA) %>%
  rbind(., lislevand_accsp)

## merge with lislevand
lislevand <- left_join(lislevand, lislevand_tax) 

## subset to bioshifts species 
lislevand <- filter(lislevand, genus_species %in% splist$species)
## from 3767 species -> 1344 species 

## save subset 
write.csv(lislevand, "data-processed/primary-traits-bioshifts/lislevand_et_al_2007.csv", row.names = FALSE)

#############################
##  Marine Species traits  ##
#############################
mspt <- read_csv("data-raw/primary-trait-data/MarineSpeciesTraits/Copy of 20220302_listAphiaID_Bioshifts_Lise_Comte.csv") %>%
  rename("genus_species" = scientificName)

## taxonomic harmonization: change name to accepted name 
mspt_tax <- harmonize_taxnomy(mspt$genus_species)

mspt_accsp <- mspt_tax[[1]]
mspt_notf <- mspt_tax[[2]]

mspt_tax <- filter(mspt, genus_species %in% mspt_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, genus, kingdom, family, order, class) %>%
  rbind(., mspt_accsp)

## merge with mspt
mspt <- left_join(mspt, mspt_tax) 

## subset to bioshifts species 
mspt <- filter(mspt, genus_species %in% splist$species)
## from 904 species -> 893 species 

## note: multiple rows per species (!)
## "measurementType" column describes traits 
## spread traits here

mspt = spread(mspt, key = measurementType, value = measurementValue) %>%
  unique(.)

## still needs work later

## save subset 
write.csv(mspt, "data-processed/primary-traits-bioshifts/MarineSpeciesTraits.csv", row.names = FALSE)


#########################
##  Pigot et al. 2020  ##
#########################
pigot <- read_excel("data-raw/primary-trait-data/Pigot/Pigot_41559_2019_1070_MOESM3_ESM.xlsx") %>%
  rename("genus_species" = Binomial) %>%
  mutate(genus = str_split_fixed(genus_species, "\\_", 2)[,1], kingdom = "Animalia", class = "Aves",
         genus_species = str_replace_all(genus_species, "\\_", " ")) 

## taxonomic harmonization: change name to accepted name 
pigot_tax <- harmonize_taxnomy(pigot$genus_species)

pigot_accsp <- pigot_tax[[1]]
pigot_notf <- pigot_tax[[2]]

pigot_tax <- filter(pigot, genus_species %in% pigot_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, genus, kingdom, class) %>%
  mutate(family = NA, order = NA) %>%
  rbind(., pigot_accsp)

## merge with pigot
pigot <- left_join(pigot, pigot_tax) 

## subset to bioshifts species 
pigot <- filter(pigot, genus_species %in% splist$species)
## from 9963 species -> 1807 species 

## save subset 
write.csv(pigot, "data-processed/primary-traits-bioshifts/Pigot_et_al_2020.csv", row.names = FALSE)


##############################
##  Storchova & Horak 2018  ##
##############################
storchova <- read_delim("data-raw/primary-trait-data/Storchova/Life-history characteristics of European birds.txt") %>%
  rename("family" = Family, "order" = Order, "genus_species" = Species) %>%
  mutate(kingdom = "Animalia", class = "Aves",
         genus = str_split_fixed(genus_species, " ", 2)[,1])

## taxonomic harmonization: change name to accepted name 
storchova_tax <- harmonize_taxnomy(storchova$genus_species)

storchova_accsp <- storchova_tax[[1]]
storchova_notf <- storchova_tax[[2]]

storchova_tax <- filter(storchova, genus_species %in% storchova_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, genus, kingdom, class, family, order) %>%
  rbind(., storchova_accsp)

## merge with storchova
storchova <- left_join(storchova, storchova_tax) 

## subset to bioshifts species 
storchova <- filter(storchova, genus_species %in% splist$species)
## from 499 species -> 319 species 

## save subset 
write.csv(storchova, "data-processed/primary-traits-bioshifts/Storchov_and_Horak_2018.csv", row.names = FALSE)


#######################
##  Sutherland 2000  ##
#######################
sutherland <- read_csv("data-raw/primary-trait-data/Sutherland/Sutherland_2000_birds.csv") %>%
  fill(colnames(.), .direction = "down") %>%
  rename("genus_species" = Species) %>%
  mutate(kingdom = "Animalia", class = "Aves",
         genus = str_split_fixed(genus_species, " ", 2)[,1])

## taxonomic harmonization: change name to accepted name 
sutherland_tax <- harmonize_taxnomy(sutherland$genus_species)

sutherland_accsp <- sutherland_tax[[1]]
sutherland_notf <- sutherland_tax[[2]]

sutherland_tax <- filter(sutherland, genus_species %in% sutherland_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, genus, kingdom, class) %>%
  mutate(family = NA, order = NA) %>%
  rbind(., sutherland_accsp)

## merge with sutherland
sutherland <- left_join(sutherland, sutherland_tax) 

## subset to bioshifts species 
sutherland <- filter(sutherland, genus_species %in% splist$species)
## from 61 species -> 58 species 

## save subset 
write.csv(storchova, "data-processed/primary-traits-bioshifts/Sutherland_2000_birds.csv", row.names = FALSE)

#################
##  TraitBank  ##
#################
traitbank <- read_csv("data-raw/primary-trait-data/TraitBank/trait_bank/traits.csv")

## clean name 
names <- str_split_fixed(traitbank$scientific_name, "<i>", 2)[,2]
names <- str_split_fixed(names, "<", 2)[,1]
traitbank$genus_species <- names
traitbank <- filter(traitbank, genus_species != "")
traitbank$genus <- str_split_fixed(traitbank$genus_species, " ", 2)[,1]

## taxonomic harmonization: change name to accepted name 
traitbank_tax <- harmonize_taxnomy(traitbank$genus_species)

traitbank_accsp <- traitbank_tax[[1]]
traitbank_notf <- traitbank_tax[[2]]

traitbank_tax <- filter(traitbank, genus_species %in% traitbank_notf) %>%
  mutate(db = NA, db_code = NA) %>%
  select(genus_species, db, db_code, genus) %>%
  mutate(family = NA, order = NA, class = NA, kingdom = NA) %>%
  rbind(., traitbank_accsp)

## merge with traitbank
traitbank <- left_join(traitbank, traitbank_tax) 

## subset to bioshifts species 
traitbank <- filter(traitbank, genus_species %in% splist$species)
## from 1001512 species -> 13062 species 

## save subset 
write.csv(traitbank, "data-processed/primary-traits-bioshifts/TraitBank.csv", row.names = FALSE)
