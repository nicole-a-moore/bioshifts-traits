## quick analysis before meeting
library(tidyverse)
library(readxl)
library(readr)
library(taxadb)

## read in list of all bioshifts species 
sp <- read_csv("data-raw/splist.csv")

## function to harmonize taxonomy 

harmonize <- function(sp_names){
  
  library(taxadb)
  library(bdc)
  source("R/clean_taxa_functions.R")
  
  sp_names = Clean_Names(sp_names,
                         return_gen_sps = TRUE)
  
  
  tofind <- data.frame(matrix(nrow = length(sp_names), ncol = 8))
  names(tofind) = c("scientificName", "kingdom", "phylum", "class", "order", "family", "db", "db_code")
  
  tofind <- data.frame(species = sp_names, tofind)
  
  tofind <- tofind %>%
    mutate(across(everything(), as.character))
  
  togo <- tofind[which(is.na(tofind$scientificName)),]
  
  mycols <- c("speciesKey", "kingdom","phylum","class","order","family", "species")
  
  ## GBIF
  # retrieve sp names
  cl <- makeCluster(detectCores()-2)
  clusterExport(cl, c("togo","standardize_taxa"))
  
  gbif_names <- pblapply(togo$species, function(x){
    try(standardize_taxa(data.frame(verbatimScientificName = x), 
                         fuzzy = FALSE,
                         silent = TRUE))
  }, cl = cl)
  
  stopCluster(cl)
  
  rem <- sapply(gbif_names, class)
  rem <- which(rem=="try-error")
  
  if(length(rem) != 0) {
    gbif_names <- gbif_names[-rem]
  }
  
  gbif_names <- rbindlist(gbif_names)
  gbif_names <- gbif_names[-which(is.na(gbif_names$scientificName)),]
  gbif_names <- gbif_names[which(gbif_names$taxonRank=='species'),]
  
  # remove duplicates
  if(any(duplicated(gbif_names$verbatimScientificName))){
    gbif_names <- gbif_names[-which(duplicated(gbif_names$verbatimScientificName)),]
  }
  
  any(!gbif_names$original_search %in% splist$species)
  
  cat("--- Summary ---\n",
      "N taxa:",nrow(togo),"\n",
      "N taxa found:",nrow(gbif_names), "\n",
      "N taxa not found:", nrow(togo)-nrow(gbif_names))
  
  gbif_names <- gbif_names[,c("verbatimScientificName","scientificName","kingdom","phylum","class","order","family","taxonID")]
  names(gbif_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
  gbif_names$db <- "gbif"
  gbif_names$db_code <- gsub("http://www.gbif.org/species/","",gbif_names$db_code)
  gbif_names$db_code <- paste("GBIF:",gbif_names$db_code,sep = "")
  
  # Feed
  tofind <- tofind %>% 
    rows_patch(gbif_names, 
               by = "species")
  
  gc()
  
  ## ITIS
  togo <- tofind[which(is.na(tofind$scientificName)),]
  
  # retrieve sp names
  itis_names <- Find_Names(spnames = togo$species,
                           db = "itis",
                           suggest_names = FALSE, # using exact names
                           return_accepted_only = TRUE) 
  
  # remove duplicates
  if(any(duplicated(itis_names$original_search))){
    itis_names <- itis_names[-which(duplicated(itis_names$original_search)),]
  }
  
  any(!itis_names$original_search %in% splist$species)
  
  itis_names <- itis_names[,c("original_search","scientificName","kingdom","phylum","class","order","family","acceptedNameUsageID")]
  names(itis_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
  itis_names$db <- "itis"
  
  # Feed
  tofind <- tofind %>% 
    rows_patch(itis_names, 
               by = "species")
  
  gc()
  
  
  ## NCBI
  togo <- tofind[which(is.na(tofind$scientificName)),]
  
  # retrieve sp names
  if(nrow(togo) == 0) {
    return(tofind)
  }
  else {
    ncbi_names <- Find_Names(spnames = togo$species,
                             db = "ncbi",
                             suggest_names = FALSE, # using exact names
                             return_accepted_only = TRUE) 
    
    # remove duplicates
    if(any(duplicated(ncbi_names$original_search))){
      ncbi_names <- ncbi_names[-which(duplicated(ncbi_names$original_search)),]
    }
    
    any(!ncbi_names$original_search %in% splist$species)
    
    ncbi_names <- ncbi_names[,c("original_search","scientificName","kingdom","phylum","class","order","family","acceptedNameUsageID")]
    names(ncbi_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
    ncbi_names$db <- "ncbi"
    
    if(nrow(ncbi_names) != 0) {
      # some names from ncbi come with duplicated genus
      # remove duplicated names
      new <- sapply(ncbi_names$scientificName, function(x){
        tmp <- strsplit(x," ")[[1]]
        if(any(duplicated(tmp))){
          paste(tmp[-1], collapse = " ")
        } else {
          x
        }
      })
      
      ncbi_names$scientificName <- new
    }
    
    # Feed
    tofind <- tofind %>% 
      rows_patch(ncbi_names, 
                 by = "species")
    
    gc()
    
    
    ## IUCN
    togo <- tofind[which(is.na(tofind$scientificName)),]
    
    # retrieve sp names
    ## PS: Cant filter on iucn method
    iucn_names<-bdc_query_names_taxadb(togo$species, 
                                       db = "iucn",
                                       suggest_names = FALSE) # using exact names
    iucn_names$scientificName <- paste(iucn_names$genus,iucn_names$specificEpithet) 
    
    # remove duplicates
    if(any(duplicated(iucn_names$original_search))){
      iucn_names <- iucn_names[-which(duplicated(iucn_names$original_search)),]
    }
    
    any(!iucn_names$original_search %in% splist$species)
    
    # select only accepted names
    iucn_names <- iucn_names[which(iucn_names$taxonomicStatus == "accepted"),]
    
    cat("--- Summary ---\n",
        "N taxa:",nrow(togo),"\n",
        "N taxa found:",length(which(iucn_names$taxonomicStatus == "accepted")), "\n",
        "N taxa not found:", nrow(togo)-length(which(iucn_names$taxonomicStatus == "accepted")))
    
    iucn_names <- iucn_names[,c("original_search","scientificName","kingdom","phylum","class","order","family","acceptedNameUsageID")]
    names(iucn_names) <- c("species","scientificName","kingdom","phylum","class","order","family","db_code")
    iucn_names$db <- "iucn"
    
    # Feed
    tofind <- tofind %>% 
      rows_patch(iucn_names, 
                 by = "species")
    
    gc()
    
    ## fix issues with Scientific names 
    new <- sapply(tofind$scientificName, function(x){
      tmp <- strsplit(x," ")[[1]]
      if(any(duplicated(tmp))){
        paste(tmp[-1], collapse = " ")
      } else {
        x
      }
    })
    
    tofind$scientificName <- new
    
    ## return
    return(tofind)
  }
}

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

## check how many species in bioshifts 
length(which(unique(co$scientificName) %in% unique(sp$scientificName))) ## 28
co_sp <- unique(co$scientificName)[which(unique(co$scientificName) %in% unique(sp$scientificName))]


## Sutherland 2000 - birds and mammals natal and breeding dispersal
suth_mamm <- read_csv("data-raw/dispersal/Sutherland_2000_mammals.csv")
colnames(suth_mamm) <- str_replace_all(colnames(suth_mamm), "\\ ", "_")
suth_mamm$Species = ifelse(suth_mamm$Species == "", NA, suth_mamm$Species)

sm_harm <- harmonize(suth_mamm$Species)

notfound <- filter(sm_harm, is.na(db_code))

## rename columns 
suth_mamm <- suth_mamm %>%
  rename("reported_name" = Species) %>%
  mutate(reported_name_fixed = reported_name)

suth_mamm <- left_join(suth_mamm, sm_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## check how many species in bioshifts 
length(which(unique(suth_mamm$scientificName) %in% unique(sp$scientificName))) ## 16
suth_mamm_sp <- unique(suth_mamm$scientificName)[which(unique(suth_mamm$scientificName) %in% unique(sp$scientificName))]


suth_bird <- read_csv("data-raw/dispersal/Sutherland_2000_birds.csv")
colnames(suth_bird) <- str_replace_all(colnames(suth_bird), "\\ ", "_")
length(unique(suth_bird$Species)) #78 spp
suth_bird$Species = ifelse(suth_bird$Species == "", NA, suth_bird$Species)

sb_harm <- harmonize(suth_bird$Species)

## fix taxonomy manually
notfound <- filter(sb_harm, is.na(db_code))
notfound$species[which(notfound$species == "Accipiter cooperi")] <- "Accipiter cooperii"
notfound$species[which(notfound$species == "Bonasa bonasia")] <- "Tetrastes bonasia"
notfound$species[which(notfound$species == "Dendragapus canadensis")] <- "Falcipennis canadensis"
notfound$species[which(notfound$species == "Actitis macularia")] <- "Actitis macularius"
notfound$species[which(notfound$species == "Speotyto cunicularia")] <- "Athene cunicularia"
notfound$species[which(notfound$species == "Delichon urbica")] <- "Delichon urbicum"
notfound$species[which(notfound$species == "Parus atricapillus")] <- "Poecile atricapillus"
notfound$species[which(notfound$species == "Parus palustris")] <- "Poecile palustris"
notfound$species[which(notfound$species == "Sitta europea")] <- "Sitta europaea"
notfound$species[which(notfound$species == "Carduelis chloris")] <- "Chloris chloris"
sb_harm_notfound <- harmonize(notfound$species)

## bind:
sb_harm <-filter(sb_harm, !is.na(db_code)) %>%
  rbind(., sb_harm_notfound)

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

## rename columns 
suth_bird <- suth_bird %>%
  rename("reported_name" = Species) 

suth_bird <- left_join(suth_bird, sb_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## check how many species in bioshifts 
length(which(unique(suth_bird$scientificName) %in% unique(sp$scientificName))) ## 66
suth_bird_sp <- unique(suth_bird$scientificName)[which(unique(suth_bird$scientificName) %in% unique(sp$scientificName))]


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

## check how many species in bioshifts 
length(which(unique(wo$scientificName) %in% unique(sp$scientificName))) ## 18
wo_sp <- unique(wo$scientificName)[which(unique(wo$scientificName) %in% unique(sp$scientificName))]


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

## check how many species in bioshifts 
length(which(unique(tamme$scientificName) %in% unique(sp$scientificName))) ## 349
tamme_sp <- unique(tamme$scientificName)[which(unique(tamme$scientificName) %in% unique(sp$scientificName))]


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

## check how many species in bioshifts 
length(which(unique(shanks9$scientificName) %in% unique(sp$scientificName))) ## 8
shanks9_sp <- unique(shanks9$scientificName)[which(unique(shanks9$scientificName) %in% unique(sp$scientificName))]


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

## check how many species in bioshifts 
length(which(unique(shanks3$scientificName) %in% unique(sp$scientificName))) ## 4
shanks3_sp <- unique(shanks3$scientificName)[which(unique(shanks3$scientificName) %in% unique(sp$scientificName))]


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

## check how many species in bioshifts 
length(which(unique(bowman$scientificName) %in% unique(sp$scientificName))) ## 6
bowman_sp <- unique(bowman$scientificName)[which(unique(bowman$scientificName) %in% unique(sp$scientificName))]


## Jenkins et al. 2007
jenkins <- read_csv("data-raw/dispersal/Jenkins_et_al_2007.csv")
colnames(jenkins) <- str_replace_all(colnames(jenkins), "\\ ", "_")
length(unique(jenkins$Scientific_Name)) #793 spp

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

## check how many species in bioshifts 
length(which(unique(jenkins$scientificName) %in% unique(sp$scientificName))) ## 319
jenkins_sp <- unique(jenkins$scientificName)[which(unique(jenkins$scientificName) %in% unique(sp$scientificName))]


## check how many unique species we have data for:
species_with_dd <- append(co_sp, wo_sp) %>%
  append(., tamme_sp) %>%
  append(., shanks9_sp) %>%
  append(., shanks3_sp) %>%
  append(., suth_mamm_sp) %>%
  append(., suth_bird_sp) %>%
  append(., bowman_sp) %>%
  append(., jenkins_sp)

length(unique(species_with_dd)) # 586 species 

## now: gather the data
co_sub = filter(co, scientificName %in% co_sp)
wo_sub = filter(wo, scientificName %in% wo_sp)
tamme_sub = filter(tamme, scientificName %in% tamme_sp)
shanks9_sub = filter(shanks9, scientificName %in% shanks9_sp)
shanks3_sub = filter(shanks3, scientificName %in% shanks3_sp)
suth_mamm_sub = filter(suth_mamm, scientificName %in% suth_mamm_sp)
suth_bird_sub = filter(suth_bird, scientificName %in% suth_bird_sp)
bowman_sub = filter(bowman, scientificName %in% bowman_sp)
jenkins_sub = filter(jenkins, scientificName %in% jenkins_sp)

cols_to_keep <- c("reported_name","reported_name_fixed", "scientificName", "kingdom", "phylum",
                  "class", "order", "family", "db", "db_code")

## Comte and Olden:
co_dd <- co_sub %>%
  select(all_of(cols_to_keep), 
         `Mean_dispersal_distance_(m)`, 
         `Maximum_dispersal_distance_(m)`, 
         `Median_dispersal_distance_(m)`, 
         Type) %>%
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
  rename("ObservationType" = Type)

## Sutherland mammals:
suth_mamm_dd <- suth_mamm_sub %>%
  select(all_of(cols_to_keep), 
         `Natal_dispersal_median_distance_(km)`,
         `Natal_dispersal_maximum_distance_(km)`, Obs_type) %>%
  mutate(ObservationType = ifelse(Obs_type == "S", "single observation",
                                  ifelse(Obs_type == "D",  "distance-density distribution", 
                                         NA))) %>%
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
  select(-Obs_type)

## Sutherland mammals:
suth_bird_dd <- suth_bird_sub %>%
  select(all_of(cols_to_keep), 
         `Natal_dispersal_median_distance_(km)`,
         `Natal_dispersal_maximum_distance_(km)`, Obs_type) %>%
  mutate(ObservationType = ifelse(Obs_type == "S", "single observation",
                                  ifelse(Obs_type == "D",  "distance-density distribution", 
                                         NA))) %>%
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
  select(-Obs_type)

## Whitmee & Orme
wo_dd <- wo_sub %>%
  select(all_of(cols_to_keep), Value, Units, Sex, Measure) %>%
  mutate(Code = ifelse(Measure == "Mean",
                       "MeanDispersalDistance", 
                       ifelse(Measure == "Maximum",
                              "MaxDispersalDistance", 
                              ifelse(Measure == "Median",
                                     "MedianDispersalDistance", 
                                     NA)))) %>%
  mutate(Unit = ifelse(str_detect(.$Units, "km"), "km",
                        ifelse(str_detect(.$Units, "Metres"), "m", NA))) %>%
  mutate(Database = "Whitmee & Orme", ObservationType = NA, Field = "Value") %>%
  rename("DispersalDistance" = Value) %>%
  filter(!is.na(DispersalDistance)) %>%
  select(-Measure, -Units)

## Tamme
tamme_dd <- tamme_sub %>%
  select(all_of(cols_to_keep), 
         `Maximum_recorded_dispersal_distance_(m)`,
         `Mean_dispersal_distance_(m)`,
         `Median_dispersal_distance_(m)`, 
         `Mode_dispersal_distance_(m)`,
         `90th_percentile_dispersal_distance_(m)`,
         `99th_percentile_dispersal_distance_(m)`,
         Data_type) %>%
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
                                                   NA)))))) %>%
  mutate(Unit = "m", Database = "Tamme", Sex = NA) %>%
  filter(!is.na(DispersalDistance)) %>%
  rename("ObservationType" = Data_type)

## Shanks 2009
shanks9_dd <- shanks9_sub %>%
  select(all_of(cols_to_keep),
         Dispersal_distance) %>%
  mutate(Unit = str_split_fixed(.$Dispersal_distance, " ", 2)[,2],
         DispersalDistance = str_split_fixed(.$Dispersal_distance, " ", 2)[,1]) %>%
  mutate(Unit = ifelse(Unit == "m1", "m", Unit), 
         DispersalDistance = str_replace_all(DispersalDistance, "\\<", "")) %>%
  filter(Unit %in% c("km", "m")) %>%
  mutate(Code = ifelse(str_detect(.$DispersalDistance, "\\–"),  "DispersalDistanceRange",
                                  "DispersalDistance")) %>%
  mutate(Database = "Shanks 2009", Sex = NA, Field = "Dispersal_distance", ObservationType = NA) %>%
  select(-Dispersal_distance) %>%
  filter(!is.na(DispersalDistance))

## Shanks 2003
shanks3_dd <- shanks3_sub %>%
  select(all_of(cols_to_keep),
         `Realized_dispersal_distance_(mean)`) %>%
  mutate(Unit = str_split_fixed(.$`Realized_dispersal_distance_(mean)`, " ", 2)[,2],
         DispersalDistance = str_split_fixed(.$`Realized_dispersal_distance_(mean)`, " ", 2)[,1]) %>%
  mutate(DispersalDistance = str_replace_all(DispersalDistance, "\\<", "")) %>%
  mutate(Code = ifelse(str_detect(.$DispersalDistance, "\\-"),  
                       "DispersalDistanceRange",
                       "DispersalDistance")) %>%
  mutate(Database = "Shanks 2003", Sex = NA, Field = "Realized_dispersal_distance_(mean)",
         ObservationType = NA) %>%
  select(-`Realized_dispersal_distance_(mean)`) %>%
  filter(!is.na(DispersalDistance))%>%
  unique()

## Bowman
bowman_dd <- bowman_sub %>%
  select(all_of(cols_to_keep),
         `Distance_(km)`) %>%
  rename("DispersalDistance" = `Distance_(km)`) %>%
  mutate(Sex = NA, ObservationType = NA, Unit = "km", Field = "Distance_(km)",
         Code = "DispersalDistance", Database = "Bowman")%>%
  filter(!is.na(DispersalDistance))%>%
  unique()
  

## Jenkins
jenkins_dd <- jenkins_sub %>%
  select(all_of(cols_to_keep), 
         `Max._Indiv._Dispersal_Distance_(m)`) %>%
  rename("DispersalDistance" = `Max._Indiv._Dispersal_Distance_(m)`) %>%
  mutate(Sex = NA, ObservationType = NA, Unit = "m", 
         Field = "Max._Indiv._Dispersal_Distance_(m)",
         Code = "MaxDispersalDistance", Database = "Jenkins") %>%
  filter(!is.na(DispersalDistance))%>%
  unique()


dd_collated <- rbind(co_dd, wo_dd) %>%
  rbind(., tamme_dd) %>%
  rbind(., shanks9_dd) %>%
  rbind(., shanks3_dd) %>%
  rbind(., suth_mamm_dd) %>%
  rbind(., suth_bird_dd) %>%
  rbind(., bowman_dd) %>%
  rbind(., jenkins_dd) %>%
  unique()


## dispersal traits
## D3 = plant dispersal traits 
## DISPERSE = dispersal traits for European aquatic macroinverts





## how many species do we have body size for?
bsize <- read_csv("data-raw/CompilationBodySizeBioshiftsv1harmonized.csv")

lifespan <- read_csv("data-raw/CompilationLifeSpan_02242022.csv")

plants <- sp %>% filter(kingdom == "Plantae") %>%
  select(-v1, -v2, -db, -db_code, -reported_name, -reported_name_fixed) %>%
  unique()

write.csv(plants, "plant_spp.csv", row.names = FALSE)



