## sorting the trait databases by taxonomic group 
## and looking at the traits studied within each taxonomic group

## goal: gather all body size data for one taxonomic group in bioshifts 
## and make list of species with missing data for Clare to search for
## also give them a list of synonyms for each species 

## Brunno says plants and birds have most data - maybe focus on those?

library(tidyverse)
library(broom)

theme_set(theme_bw())

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







