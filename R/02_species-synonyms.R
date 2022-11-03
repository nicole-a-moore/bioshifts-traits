## get list of synonyms for all accepted species names in bioshifts V1 and V2
library(taxadb)
library(tidyverse)

## read in harmonized species list 
sp_list <- read.csv("data-raw/splist.csv")

## get vector of accepted species names 
acc_names <- sp_list$species

length(unique(acc_names)) ## some are duplicated - get rid 
acc_names <- unique(acc_names)

## try getting all synonyms for every species
getname = filter_name(acc_names, "itis") %>%
  select(acceptedNameUsageID, input)

itis = taxa_tbl("itis") %>%
  collect() %>%
  filter(acceptedNameUsageID %in% unique(getname$acceptedNameUsageID)[which(!is.na(unique(getname$acceptedNameUsageID)))]) %>%
  select(scientificName, acceptedNameUsageID, phylum, class, order, family, genus) %>% 
  rename("synonym" = scientificName) %>%
  left_join(., getname) %>%
  rename("accepted_name" = "input")

saveRDS(itis, "data-processed/itis-species-synonyms.rds")

getname = filter_name(acc_names, "col") %>%
  select(acceptedNameUsageID, input)

col = taxa_tbl("col") %>%
  collect() %>%
  filter(acceptedNameUsageID %in% unique(getname$acceptedNameUsageID)[which(!is.na(unique(getname$acceptedNameUsageID)))]) %>%
  select(scientificName, acceptedNameUsageID, phylum, class, order, family, genus) %>% 
  rename("synonym" = scientificName)%>%
  left_join(., getname) %>%
  rename("accepted_name" = "input")

saveRDS(col, "data-processed/col1-species-synonyms.rds")

getname = filter_name(acc_names, "gbif") %>%
  select(acceptedNameUsageID, input)

gbif = taxa_tbl("gbif") %>%
  collect() %>%
  filter(acceptedNameUsageID %in% unique(getname$acceptedNameUsageID)[which(!is.na(unique(getname$acceptedNameUsageID)))]) %>%
  select(scientificName, acceptedNameUsageID, phylum, class, order, family, genus) %>% 
  rename("synonym" = scientificName) %>%
  left_join(., getname) %>%
  rename("accepted_name" = "input")

saveRDS(gbif, "data-processed/gbif-species-synonyms.rds")

getname = filter_name(acc_names, "ncbi") %>%
  select(acceptedNameUsageID, input)

ncbi = taxa_tbl("ncbi") %>%
  collect() %>%
  filter(acceptedNameUsageID %in% unique(getname$acceptedNameUsageID)[which(!is.na(unique(getname$acceptedNameUsageID)))]) %>%
  select(scientificName, acceptedNameUsageID, phylum, class, order, family, genus) %>% 
  rename("synonym" = scientificName) %>%
  left_join(., getname) %>%
  rename("accepted_name" = "input")

saveRDS(ncbi, "data-processed/ncbi-species-synonyms.rds")

syn_db <- rbind(itis, ncbi) %>% rbind(., col) %>% rbind(., gbif) %>%
  unique(.)

## save the list
saveRDS(syn_db, "data-processed/species-synonyms.rds")

syn_db <- readRDS("data-processed/species-synonyms.rds")
write.csv(syn_db, "data-processed/species-synonyms.csv", row.names = F)

## clean the data to rid of the ? in some GBIF species synonyms
## replace the ? with the contents of the 'genus' column




