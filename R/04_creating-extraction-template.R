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



### -----------------
## subsetting to cadillac data 
source("R/taxonomic-harmonization/clean_taxa_functions.R")

## read in cadillac data:
cadillac = read.csv("data-raw/bioshiftsv1/cadillac_data.csv") 

## get list of species in cadillac data 
cadillac$reported_name = Clean_Names(gsub("_", " ", cadillac$publi_name), return_gen_sps = F)

length(unique(cadillac$reported_name)) ## 781 species in cadillac 

## get rid of insects, amphibians 
cadillac <- filter(cadillac, !class %in% c("Insecta", "Amphibia"))

length(unique(cadillac$reported_name)) ## 708 species are not amphibians or insects

## filter out species we already have body length for 
sp <-  read.csv("data-raw/splist.csv") %>%
  select(reported_name, scientificName)

cadillac = left_join(cadillac, sp)

cadillac <- filter(cadillac, !scientificName %in% bsize$SpeciesChecked)

length(unique(cadillac$scientificName)) ## 157 cadillac species are missing body size  


## mark them in Clare's database 
## read in Clare's trait data
## downloaded on Nov 18 from
# https://docs.google.com/spreadsheets/d/1TGZmFGgAzilT1WXmUxL2XHLKSvj7kdG2/edit?rtpof=true&sd=true
trait = read.csv("data-raw/BodyLength_Clare.csv")

sp <- read.csv("data-raw/splist.csv") %>%
  filter(!class %in% c("Insecta", "Amphibia")) 

trait <- left_join(trait, sp, by = c("Species" = "species", "db_code", "db"))

cad_sp <- filter(trait, scientificName %in% cadillac$scientificName)

length(unique(cad_sp$scientificName))

which(!cadillac$scientificName %in% trait$scientificName)

## add column identifying cadillac-ness 

trait <- trait %>%
  mutate(cadillac = ifelse(scientificName %in% cad_sp$scientificName, "y", "n"))

length(which(trait$cadillac == "y"))

## get rid of duplicates
length(unique(trait$scientificName))

trait = trait %>%
  select(-Species, -kingdom, -class, -phylum, -order, -family, -db_code, -db,
         -reported_name, -reported_name_fixed, -v1, -v2) %>%
  select(scientificName, cadillac, everything()) %>%
  unique()

length(unique(test$scientificName))

## write:
write.csv(trait, "data-processed/files for Clare/BodyLength_Clare_cadillac.csv", row.names = F)

