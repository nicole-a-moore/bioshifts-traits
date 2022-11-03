## creating a database for someone to fill in body length data
library(tidyverse)

## read in the database Lise gave me of all species we have body size data for
bsize <- read.csv("data-raw/CompilationBodySizeBioshiftsv1harmonized.csv")

## read in list of all bioshifts species 
sp <- read.csv("data-raw/splist.csv")

## filter out insects and amphibians (Sarah's collaborators are already working on these)
sp <- sp %>%
  filter(!class %in% c("Insecta", "Amphibia")) 

## filter out species we already have body length for 
sp <- sp %>%
  filter(!species %in% bsize$SpeciesChecked)

## get rid of not accepted name column, other columns we don't need
sp <- sp %>%
  select(-scientificName, -v1, -v2, -reported_name, -reported_name_fixed) %>%
  unique()

## make a template for search for body size 
colnames(bsize)

sp[,9:20] <- NA

colnames(sp)[9:20] <- c("BodyLength", "Unit", "MeasurementType", "Metric", "NumberOfIndv", "Sex", "Age", "Lifestage", "Level",
                        "SearchedBy", "Notes", "Citation")

## write out
write.csv(sp, "data-processed/targeted-search-template/BodyLength.csv", row.names = FALSE)

