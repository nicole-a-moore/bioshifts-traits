## script to check availability of traits of interest for plants
## renames columns in primary trait databases to common names
library(tidyverse)

## read in the list of traits in each database:
traits <- read_csv("data-raw/DetailTraits.csv")

## break into one dataframe for each database 
trait_list <- t(traits) %>%
  split(., seq(nrow(.)))

names(trait_list) <- colnames(traits)

## subset to ONLY databases for plants 
names(trait_list)
plant_dbs <- str_split(data_key$possible_databases[which(data_key$group == "Aquatic and terrestrial plants")], "\\, ")[[1]]

## make sure all are in the list of traits
length(which(!plant_dbs %in% names(trait_list)))
plant_dbs[which(!plant_dbs %in% names(trait_list))]
## they are

## next: see what traits we have for plants
plants <- trait_list[names(trait_list) %in% plant_dbs]

plants <- data.frame(do.call(rbind, plants))
plants$Database_name = rownames(plants)

## make data base with columns:
## group 
## trait
## Database_name 
## trait_common_name
## trait_description
## useful_y_n
plants <- gather(plants, key = "old_col", value = "trait", c(1:2091)) %>%
  select(-old_col) %>%
  filter(!is.na(trait)) %>%
  unique(.) %>%
  arrange(Database_name) %>%
  mutate(trait_common_name = NA,
         useful_y_n = NA) 

## 2356 traits 

## write out and start to decide on common names
write.csv(plants, "data-processed/traits-by-group/plants.csv", row.names = FALSE)

## Notes: 

##BIOTIME:
## data is from TRY, units are unclear
## right now guessing - deep dive later to make sure units are not mixed in TRY




