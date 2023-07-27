## using trait-dispersal distance relationships to predict dispersal distance for plant and bird species in Bioshifts
library(tidyverse)
library(nlme)

## read function to harmonize taxonomy 
source("R/harmonize.R")
library(taxadb)
library(parallel)
library(pbapply)
library(traitdataform)
library(data.table)

select = dplyr::select

###################
####   plants  ####
###################
## use models from Tamme et al.:
## package dispeRsal
## https://figshare.com/collections/Predicting_species_maximum_dispersal_distances_from_simple_plant_traits/3306522

## try for best model (R2 = 0.60)
## dispersal syndrome, growth form, and terminal velocity 
## for plant species with dispersal syndromes that are not ant/animal

## otherwise use next best model (R2 = 0.53) 
## dispersal syndrome, growth form, species taxonomy data
## for plant species with dispersal syndromes are not ant/animal and for other species without known terminal velocity

## predicted values in general overestimate for spp with short-distance dispersal and underestimate for spp with long-distance dispersal

#########################
####  gather traits  ####
#########################
## read in bioshifts
v1 <- read.csv("data-processed/corrected-bioshifts_fixed.csv")

## filter to plants 
v1 <- v1 %>% filter(Kingdom == "Plantae") ## 91649 obs

## start with LEDA
## growth form (GF):
ledagf = read.delim("data-raw/primary-trait-data/LEDA/plant_growth_form copy.txt", sep = ";")

tamme <- read.csv("data-raw/dispersal/Tamme_DispersalDistanceData.csv")

## combine to get a growth form key
key_gf <- inner_join(ledagf, tamme, by = c("SBS.name" = "Species")) %>%
  select(plant.growth.form, Growth_form) %>%
  distinct() 

herb_key <- key_gf$plant.growth.form[which(key_gf$Growth_form == "herb")]
shrub_key <- key_gf$plant.growth.form[which(key_gf$Growth_form == "shrub")]

both <- herb_key[which(herb_key %in% shrub_key)] 
## these could be considered herb or shrub - remove for now
herb_key <- herb_key[which(!herb_key %in% both)]
shrub_key <- shrub_key[which(!shrub_key %in% both)]

# https://hosho.ees.hokudai.ac.jp/tsuyu/top/dct/lf.html
# plant life form classifications:

# Phanerophytes, epiphyte, sclerophyte - tree
# Chamaephyte - shrub

tree_key <- c("Phanerophytes", "Epiphyte", "Sclerophyte")
herb_key <- append(herb_key, "Chamaephyte")

unique(ledagf$plant.growth.form)
## exclude parasites since likely not able to accurately predict dispersal distance

## clean data so everything is tree shrub or herb
ledagf_final = ledagf %>%
  filter(!plant.growth.form %in% c("Vascular semi-parasite", "Vascular parasite")) %>% 
  mutate(GF = ifelse(plant.growth.form %in% herb_key, "herb", 
                     ifelse(plant.growth.form %in% tree_key, "tree",  
                            ifelse(plant.growth.form %in% shrub_key, "shrub", NA))))%>%
  select(SBS.name, GF) %>%
  rename("Taxon" = SBS.name) %>%
  filter(!is.na(GF)) %>%
  unique()

## get rid of ones with two growth form classifications for now
ledagf_final <- filter(ledagf_final, ! Taxon %in% ledagf_final$Taxon[which(duplicated(ledagf_final$Taxon))])

## filter to bioshifts species
ledagf_sub <- filter(ledagf_final, Taxon %in% v1$scientificName)
length(unique(ledagf_sub$Taxon)) # 1264 spp

## Terminal velocity
ledatv = read.delim("data-raw/primary-trait-data/LEDA/TV_2016 copy.txt", sep = ";")

## take maximum terminal velocity (as Tamme did) 
ledatv_final = ledatv %>%
  filter(!general.method == "unknown") %>%
  group_by(SBS.name) %>%
  mutate(TV = max(single.value..m.s.)) %>%
  ungroup() %>%
  rename("Taxon" = SBS.name) %>%
  select(Taxon, TV) %>%
  unique() 

## filter to bioshifts species
ledatv_sub <- filter(ledatv_final, Taxon %in% v1$scientificName)  
length(unique(ledatv_sub$Taxon)) # 797 spp

## Dispersal type
ledadt = read.delim("data-raw/primary-trait-data/LEDA/dispersal_type copy.txt", sep = ";")

ledadt = filter(ledadt, dispersal.type != "") # get rid of empty data

unique(ledadt$gen..dispersal.type)
## get rid of hemerochory, nautochory
ledadt <- filter(ledadt, !gen..dispersal.type %in% c('hemerochor', 'nautochor'))
unique(ledadt$gen..dispersal.type)

## filter to bioshifts species
ledadt_sub <- filter(ledadt, SBS.name %in% v1$scientificName)

## clean
ledadt_sub <- ledadt_sub %>% 
  mutate(DS = ifelse(dispersal.type == "dysochor" & str_detect(.$dispersal.vector, "ant"), 
                     "ant",
                     ifelse(gen..dispersal.type == "zoochor",
                            "animal",
                            ifelse(gen..dispersal.type == "autochor", "ballistic",
                                   ifelse(gen..dispersal.type == "meteorochor", "wind", NA))))) %>%
  rename("Taxon" = SBS.name) %>%
  select(Taxon, DS) %>%
  filter(!is.na(DS)) %>%
  unique() 

length(unique(ledadt_sub$Taxon)) # 1524 spp
## note: will need to go back and see which type of wind dispersal features each species has 

## combine all leda subsets: 
leda_sub = full_join(ledatv_sub, ledadt_sub) %>%
  full_join(., ledagf_sub) %>%
  mutate(source = "LEDA")


## EuDiS 
eudis <- read.csv("data-raw/dispersal/EuDiS.csv")

## dispersal syndrome (DS)
eudis_ds <- eudis %>%
  mutate(animal = ifelse(Endozoochorous == 1 | Epizoochorous == 1, 
                         "animal",
                         "No"),
         ant = ifelse(Myrmecochorous == 1, 
                      "ant",
                      "No"),
         wind.special = ifelse(Anemochorous == 1, 
                      "wind.special",
                      "No"),
         ballistic = ifelse(Ballochorous == 1, 
                               "ballistic",
                               "No")) %>%
  gather(key =  "DS_type", value = "DS", c(ant, animal, wind.special, ballistic)) %>%
  select(Species_Flora_Europaea, DS) %>%
  filter(DS != "No") %>%
  rename("Taxon" = Species_Flora_Europaea) %>%
  unique()

## filter to bioshifts spp
eudis_ds_sub <- filter(eudis_ds, Taxon %in% v1$scientificName)

length(unique(eudis_ds_sub$Taxon)) # 819 spp

## growth form (GF)
eudis_gf <- eudis %>%
  filter(Tree_or_Shrub == 0) %>%
  mutate(GF = "herb") %>%
  select(Species_Flora_Europaea, GF) %>%
  rename("Taxon" = Species_Flora_Europaea) %>%
  unique()

## filter to bioshifts spp
eudis_gf_sub <- filter(eudis_gf, Taxon %in% v1$scientificName)

length(unique(eudis_gf_sub$Taxon)) # 1615 spp


## combine all subsets: 
eudis_sub = full_join(eudis_gf, eudis_ds_sub) %>%
  mutate(source = "EuDiS", TV = NA)


## BIEN
library(BIEN)

## download all trait data for each bioshift species
# bien = BIEN_trait_species(unique(v1$scientificName))

# write.csv(bien, "data-processed/predicting-dispersal-distance/BIEN_bioshifts.csv", row.names = FALSE)
bien = read.csv("data-processed/predicting-dispersal-distance/BIEN_bioshifts.csv")
unique(bien$trait_name)

## growth form (GF)
bien_gf <- bien %>%
  filter(trait_name == "whole plant growth form") %>%
  mutate(trait_value = str_remove_all(trait_value, "\\*")) %>%
  mutate(trait_value = str_remove_all(trait_value, "\\-")) %>%
  mutate(GF = ifelse(trait_value %in% c("tree", "Tree", "Small_Tree", "tree ", "giant tree"), "tree", 
                     ifelse(trait_value %in% c("Herb", "herb", "Scandent_Herb","climbing herb", "twining herb",
                                               "herbaceous climber","Grass", "grass", "Trailing_Herb",
                                               "forb", "fern"), "herb", 
                            ifelse(trait_value %in% c("Shrub", "subshrub", "shrub", "shrublet"), "shrub", NA))))

## filter to bioshifts spp
biengf_sub <- filter(bien_gf, scrubbed_species_binomial %in% v1$scientificName)
length(unique(biengf_sub$scrubbed_species_binomial)) # 4909 spp

is_na <- filter(biengf_sub, is.na(GF))
not_na <- filter(biengf_sub, !is.na(GF))

need_resolving <- filter(is_na, !scrubbed_species_binomial %in% not_na$scrubbed_species_binomial)
unique(need_resolving$trait_value)
## okay none really fit naturally into the categories, so exclude them

biengf_sub <- filter(biengf_sub, !is.na(GF)) %>%
  select(scrubbed_species_binomial, GF) %>%
  filter(!is.na(GF)) %>%
  unique() %>%
  rename("Taxon" = scrubbed_species_binomial) 

length(unique(biengf_sub$Taxon)) # 4817 spp

## reconcile species with multiple growth forms 
bien_multi = biengf_sub %>%
  filter(Taxon %in% .$Taxon[which(duplicated(.$Taxon))])
## leave these for later 

biengf_sub <- filter(biengf_sub, !duplicated(Taxon))
length(unique(biengf_sub$Taxon)) # 4817 spp

bien_sub <- biengf_sub %>%
  mutate(TV = NA, DS = NA, source = "BIEN")

## BROT
brot = read.csv("data-raw/primary-trait-data/BROT/BROT2_dat.csv")
unique(brot$Trait)
## GrowthForm, DispMode
brot = brot %>%
  filter(Trait %in% c("GrowthForm", "DispMode", "SeedMass"))

## GrowthForm
gf = brot %>%
  filter(Trait == "GrowthForm") 
unique(gf$Data)

## clean data so everything is tree shrub or herb
gf = gf %>%
  mutate(GF = ifelse(Data %in% c("large shrub", "subshrub", "shrub"), "shrub", 
                     ifelse(Data %in% c("perennial forb", "annual forb", "variable forb",
                                        "perennial graminoid", "variable graminoid",
                                        "annual graminoid", "geophyte"), "herb", 
                            ifelse(Data == "tree", "tree", NA)))) %>%
  filter(!is.na(GF))

## DispMode
dm = brot %>%
  filter(Trait == "DispMode") 
unique(dm$Data)

# - G: autochory, by Gravity (= unassisted dispersal).
# - W: anemochory, by Wind (with wind dispersal adaptations).
# - H: Hydrochory, by water.
# - B: Ballistichory, by launching (= ballochory).
# - M: Myrmecochory, by ants.
# - N: eNdozoochory, internal animal transport.
# - P: ePizoochory, external animal transport (= exozoochory).
# - O: hOarding, scatter and hoarding diaspores by animals (others than ants).
# - Z: Zoochory, dispersal mediated by animals (unknown transport system).

dm = dm %>%
  ## some species have more than one mode; split into multiple columns
  mutate(Data1 = str_split_fixed(dm$Data, "", 4)[,1],
         Data2 = str_split_fixed(dm$Data, "", 4)[,2],
         Data3 = str_split_fixed(dm$Data, "", 4)[,3],
         Data4 = str_split_fixed(dm$Data, "", 4)[,4]) %>%
  gather(key = "type", value = "Data", c(Data1, Data2, Data3, Data4)) %>%
  select(-type) %>%
  mutate(DS = ifelse(Data == "W", "wind.special", 
                     ifelse(Data == "G", "wind.none", 
                            ifelse(Data == "B", "ballistic", 
                                   ifelse(Data == "M", "ant",
                                          ifelse(Data %in% c("N", "P", "O", "Z"),"animal",
                                                 NA)))))) %>%
  filter(!is.na(DS)) 

## make data so that is has 1 row per species per unique trait combination
brot_final <- dm %>%
  select(Taxon, DS) %>%
  full_join(select(gf, c(Taxon, GF))) 

## filter to bioshifts species
brot_sub <- filter(brot_final, Taxon %in% v1$scientificName)%>%
  mutate(TV = NA, source = "BROT")
length(unique(brot_sub$Taxon)) # 722 spp


## D3
d3 = read.delim("data-raw/primary-trait-data/D3/1-s2.0-S1433831913000218-mmc1.txt")
colnames(d3)
# vterm = terminal velocity
# "citation_prop_ane" # wind = but don't know whether gravity or special
# "citation_prop_dyso" # animal
# "citation_prop_endo" # animal
# "citation_prop_epi" # animal

d3_final <- d3 %>%
  mutate(DS = ifelse(citation_prop_dyso > 0, "animal",
                     ifelse(citation_prop_endo > 0, "animal", 
                            ifelse(citation_prop_epi > 0, "animal",
                                   NA)))) %>%
  mutate(TV = vterm) %>%
  mutate(Taxon = paste(str_split_fixed(name, " ", 2)[,1], str_split_fixed(name, " ", 3)[,2], sep = " ")) %>%
  select(Taxon, DS, TV) %>%
  filter(!is.na(DS) | !is.na(TV)) %>%
  mutate(GF = NA)

## filter to bioshifts species
d3_sub <- filter(d3_final, Taxon %in% v1$scientificName) %>%
  mutate(source = "D3")
length(unique(d3_sub$Taxon)) ## 1207


## combine all
all <- rbind(brot_sub, bien_sub) %>%
  rbind(., leda_sub) %>%
  rbind(., d3_sub) %>%
  rbind(., eudis_sub)

## resolve multiple growth forms per species 
gf_all <- select(all, Taxon, GF, source) %>%
  distinct() %>%
  filter(!is.na(GF)) 
length(unique(gf_all$Taxon)) # 12187

tally <- gf_all %>%
  group_by(Taxon) %>%
  tally(length(unique(GF))) %>%
  ungroup() %>%
  group_by(n)%>% 
  tally()
## most species have 1 growth form reported 

## for ones with more than one growth form, prioritze leda, then EuDIS, then BROT, then BIEN
gf_2sources_all <- gf_all %>%
  filter(!is.na(GF)) %>%
  group_by(Taxon) %>%
  filter(Taxon %in% gf_all$Taxon[duplicated(gf_all$Taxon)] & length(unique(GF)) != 1) %>%
  ungroup() 
length(unique(gf_2sources_all$Taxon)) ## 114

gf_2sources_leda <- filter(gf_2sources_all, source == "LEDA")
length(unique(gf_2sources_leda$Taxon)) ## 66

gf_2sources <- filter(gf_2sources_all, !Taxon %in% gf_2sources_leda$Taxon)

gf_2sources_eudis <- filter(gf_2sources, source == "EuDiS")
length(unique(gf_2sources_eudis$Taxon)) ## 27

gf_2sources <- filter(gf_2sources, !Taxon %in% gf_2sources_eudis$Taxon)

gf_2sources_brot <- filter(gf_2sources, source == "BROT")
length(unique(gf_2sources_brot$Taxon)) ## 21

## get rid of repeats
gf_all_norep <- gf_all %>%
  group_by(Taxon) %>%
  filter(!duplicated(GF)) %>%
  ungroup()

gf_all <- filter(gf_all_norep, !Taxon %in% gf_2sources_all$Taxon) %>%
  rbind(gf_2sources_brot, .) %>%
  rbind(., gf_2sources_leda) %>%
  rbind(., gf_2sources_eudis) %>%
  rename("source_GF" = source)
length(unique(gf_all$Taxon)) ## 12187


## resolve multiple terminal velocities per species 
tv_all <- select(all, Taxon, TV, source) %>%
  distinct()%>%
  filter(!is.na(TV)) 
length(unique(tv_all$Taxon)) # 1244

tally <- tv_all %>%
  filter(!is.na(TV)) %>%
  group_by(Taxon) %>%
  tally(length(unique(TV))) %>%
  ungroup() %>%
  group_by(n)%>% 
  tally()

## for ones with more than one growth form, keep maximum
tv_2sources_all <- tv_all %>%
  filter(!is.na(TV)) %>%
  group_by(Taxon) %>%
  filter(Taxon %in% tv_all$Taxon[duplicated(tv_all$Taxon)] & length(unique(TV)) != 1) %>%
  ungroup() 
length(unique(tv_2sources_all$Taxon)) ## 531

tv_2sources <- tv_2sources_all %>%
  group_by(Taxon) %>%
  filter(TV == max(TV)) %>%
  ungroup()
length(unique(tv_2sources$Taxon)) ## 531

tv_all <- filter(tv_all, !Taxon %in% tv_2sources_all$Taxon) %>%
  rbind(tv_2sources, .) %>%
  rename("source_TV" = source)
length(unique(tv_all$Taxon)) ## 1244

## for ones with more than one dispersal mode, let them keep it 
## but resolve wind dispersal type: 
ds_all <- select(all, Taxon, DS, source) %>%
  distinct() %>%
  filter(!is.na(DS)) 
length(unique(ds_all$Taxon)) #1798

tally <- ds_all %>%
  group_by(Taxon) %>%
  tally(length(unique(DS))) %>%
  ungroup() %>%
  group_by(n)%>% 
  tally()

## get rid of repeats
ds_all_norep <- ds_all %>%
  filter(!is.na(DS)) %>%
  group_by(Taxon) %>%
  filter(!duplicated(DS)) %>%
  ungroup()

wind <- filter(ds_all, str_detect(DS, "wind"))

tally <- wind %>%
  filter(!is.na(DS)) %>%
  group_by(Taxon) %>%
  tally(length(unique(DS))) %>%
  ungroup() %>%
  group_by(n)%>% 
  tally()

wind <- filter(wind, DS != "wind")

ds_all <- filter(ds_all_norep, !Taxon %in% wind$Taxon) %>%
  rbind(wind, .) %>%
  rename("source_DS" = source)
length(unique(ds_all$Taxon)) ## 1798


## combine them all back together:
## keep all the unique dispersal syndromes 
## keep only one terminal velocity 
## keep only one growth form

all_clean <- left_join(ds_all, gf_all) %>%
  left_join(tv_all)

all_clean %>%
  filter(Taxon %in% .$Taxon[duplicated(Taxon)]) %>%View

## see how much coverage we have 
## how many bioshifts plants are we completely missing info for?
length(which(!unique(v1$scientificName) %in% all_clean$Taxon)) ## 4942 spp
length(unique(v1$scientificName)) ## out of 6740

## how many do we have all 3 traits for?
all_clean %>%
  filter(!is.na(GF), !is.na(TV), !is.na(DS)) %>%
  select(Taxon) %>% nrow(.)
## 1771 spp

all_clean %>%
  select(Taxon) %>% nrow(.)
## out of 2559


# check
all_clean %>%
  select(Taxon, GF, source_GF) %>%
  distinct() %>%
  group_by(Taxon) %>%
  tally(length(unique(GF))) %>%
  ungroup() %>%
  group_by(n)%>% 
  tally()

all_clean %>%
  select(Taxon, TV, source_TV) %>%
  distinct() %>%
  group_by(Taxon) %>%
  tally(length(unique(TV))) %>%
  ungroup() %>%
  group_by(n)%>% 
  tally()

all_clean %>%
  select(Taxon, DS, source_DS) %>%
  distinct() %>%
  group_by(Taxon) %>%
  tally(length(unique(DS))) %>%
  ungroup() %>%
  group_by(n)%>% 
  tally()


## to do later:
## ~~~~~~~~~~~~
## figure out how to add species in BIEN with multiple GFs
## add other traits to use more predictive models when possible?
## resort to TRY for missing values?
## figure out what type of wind dispersal species have
length(which(all_clean$DS == "wind"))
## for now, exclude
all_clean <- all_clean %>%
  filter(DS != "wind")

## yay! now have a database we can use to make predictions 
############################
####  make predictions  ####
############################
## load required packages 
library(nlme)

## load the functions
load("R/dispeRsal.rda")

## modify the function to:
## 1. avoid broken link
## 2. allow multiple predictions per species (for different dispersal syndromes)
dispeRsal <- function(data.predict, model, CI = F, random = T, tax = "family",
                      write.result=F){
  require(nlme)
  require(AICcmodavg)
  thisVersion <- 0.2
  #currentVersion <- scan("http://www.botany.ut.ee/dispersal/version.txt")
  TPL_nomatch <- NULL
  
  if (is.numeric(model)){
    traits <- switch(model, c("DS", "GF","TV"), c("DS", "GF", "SM","RH"), c("DS", "GF", "RH"),
                     c("DS", "GF", "SM"), c("DS","GF"))
  } else{
    traits <- model
  }
  
  if(any(is.na(match(traits, names(data.predict))))) {
    stop("variable names do not match required names")
  }
  
  if(any(is.na(match(levels(data.predict$GF), levels(model.data$GF))))) {
    stop("wrong category name in column GF")
  }
  
  if(any(is.na(match(levels(data.predict$DS), levels(model.data$DS))))) {
    stop("wrong category name in column DS")
  }
  
  model.data.traits <- model.data[,2:6]
  model.data.rest <- model.data[,c(1,7:9)]
  ind <- names(model.data) %in% traits
  md <- na.omit(data.frame(model.data.rest, model.data[ind]))
  pd <- data.predict[,na.omit(match(names(md),
                                    names(data.predict)))]
  if(random == TRUE){
    pd.species.tpl <- TPLMod(as.character(pd$Species))
    pd$FamilyTPL <- pd.species.tpl$Family
    TPL_nomatch <- pd[which(model.data$FamilyTPL=="" |
                              is.na(model.data$FamilyTPL)),]$Species
    if (any(pd$FamilyTPL=="" | is.na(pd$FamilyTPL)) ==F) {
      pd <- pd
    } else{
      pd <- pd[-which(pd$FamilyTPL=="" | is.na(pd$FamilyTPL)),]
    }
    ln <- length(pd$FamilyTPL)
    Order <- vector("character", ln)
    for (i in 1:ln) {
      Order[i] <- as.character(OrderFamilies$Order[match(pd$FamilyTPL[i],
                                                         OrderFamilies$Family)])
    }
    pd$OrderTPL <- Order
  }
  
  pd <- na.omit(pd)
  
  if(dim(pd)[1] == 0) {
    stop(paste("not enough data to run specified model with traits", paste(traits, collapse = " + ")))
  }
  
  formu <- as.formula(paste("logMax ~", paste(traits, collapse="+")))
  m <- lm(formu,  data=md)
  
  if(CI == T) {
    log10MDD <- predict(m, na.omit(pd), interval = "confidence")
    colnames(log10MDD) <- c("log10MDD", "log10MDD_lwrCL", "log10MDD_uppCL")
  } else{
    log10MDD <- predict(m, na.omit(pd))
  }
  
  Species <- as.character(na.omit(pd)$Species)
  assign("formu", formu,  envir = .GlobalEnv)
  
  if(random == T) {
    
    if (tax == "family") {
      rs <- reStruct(object = ~ 1 | OrderTPL/FamilyTPL, pdClass="pdDiag")
      level.max <- 2
      m <- lme(formu, data = md, random = rs)
      m$call$fixed <- as.call(formu)
      
      if (dim(na.omit(pd))[1] == 1) {
        log10MDD_Family <- data.frame(predict(m, na.omit(pd),
                                              level = 0:level.max))[3, 3]
        log10MDD_Order <- data.frame(predict(m, na.omit(pd),
                                             level = 0:level.max))[2, 3]
      } else{
        log10MDD_Family <- data.frame(predict(m, na.omit(pd),
                                              level = 0:level.max))[,-c(1,2)][,3]
        log10MDD_Order <- data.frame(predict(m, na.omit(pd),
                                             level = 0:level.max))[,-c(1,2)][,2]
      }
      
      if (CI == TRUE){
        se <- predictSE.lme(m, na.omit(pd))$se.fit
        log10MDD_Family <- as.data.frame(log10MDD_Family)
        log10MDD_Family$log10MDD_Family_lwrCL <-
          log10MDD_Family$log10MDD_Family - 2 * se
        log10MDD_Family$log10MDD_Family_uppCL <-
          log10MDD_Family$log10MDD_Family + 2 * se
        log10MDD_Order <- as.data.frame(log10MDD_Order)
        log10MDD_Order$log10MDD_Order_lwrCL <-
          log10MDD_Order$log10MDD_Order - 2 * se
        log10MDD_Order$log10MDD_Order_uppCL <-
          log10MDD_Order$log10MDD_Order + 2 * se
      }
      
      if (any(is.na(log10MDD_Family))) {
        p = data.frame(Species, log10MDD_Family, log10MDD_Order)
      } else{
        p = data.frame(Species, log10MDD_Family)
      }
      
      if (any(is.na(log10MDD_Order))) {
        p = data.frame(Species, log10MDD_Family, log10MDD_Order, log10MDD)
      }
    }
    
    if (tax == "order") {
      rs <- reStruct(object = ~ 1 | OrderTPL, pdClass="pdDiag")
      level.max <- 1
      m <- lme(formu, data = md, random = rs)
      m$call$fixed <- as.call(formu)
      
      if (dim(na.omit(pd))[1] == 1) {
        log10MDD_Order <- data.frame(predict(m, na.omit(pd),
                                             level = 0:level.max))[2,2]
      } else{
        log10MDD_Order <- data.frame(predict(m, na.omit(pd),
                                             level = 0:level.max))[,-c(1,2)]
      }
      
      if (CI == TRUE){
        se <- predictSE.lme(m, na.omit(pd))$se.fit
        log10MDD_Order <- as.data.frame(log10MDD_Order)
        log10MDD_Order$log10MDD_Order_lwrCL <-
          log10MDD_Order$log10MDD_Order - 2 * se
        log10MDD_Order$log10MDD_Order_uppCL <-
          log10MDD_Order$log10MDD_Order + 2 * se
      }
      
      if (any(is.na(log10MDD_Order))) {
        p = data.frame(Species, log10MDD_Order, log10MDD)
      } else{
        p = data.frame(Species, log10MDD_Order)
      }
    }
    
  } else{
    p <-  data.frame(Species, log10MDD)
  }
  
  if(any(traits == "DS")) {
    DS <- pd[1:length(p$Species),]$DS
  } else {
    DS <- NULL
  }
  SpecDS.p <- paste(p$Species, DS)
  SpecDS.md <- paste(md$Species, md$DS)
  md.s <- md[which(match(SpecDS.md, SpecDS.p)!="NA"),]
  p.s <- p[which(match(SpecDS.p, SpecDS.md)!="NA"),]
  r <-  rownames(p[which(match(SpecDS.p, SpecDS.md)!="NA"),])
  if(length(r) > 0){
    log10MDD_measured <- as.vector(rep(NA, dim(p)[1]), mode="numeric")
    p <- cbind(p, log10MDD_measured)
    p[as.character(r),]$log10MDD_measured <- md.s[match(p.s$Species, md.s$Species),]$logMax
  } 
  
  if(random == T) {
    Order <- pd[na.omit(match(p$Species, pd$Species)),]$OrderTPL
    Family <- pd[na.omit(match(p$Species, pd$Species)),]$FamilyTPL
    if(is.null(DS)) {
      p <- cbind(p, Order, Family)
      cs <- dim(p)[2]
      p <- p[,c(1, cs-1, cs, c(2:(cs-2)))]
      p <- p[order(p$Family, p$Species),]    
    } else {
      p <- cbind(p, DS, Order, Family) 
      cs <- dim(p)[2]
      p <- p[,c(1, cs-1, cs, cs-2, c(2:(cs-3)))]
      p <- p[order(p$Family, p$Species),]
    }
    
  } else{
    p <- cbind(p, DS)
    cs <- dim(p)[2]
    p <- p[,c(1, cs, c(2:(cs-1)))]
  }
  
  remove(formu,  envir = .GlobalEnv)
  
  # if (thisVersion != currentVersion) {
  #   ver = cat( "\n\n\nATTENTION\n",
  #              "You are using dispeRsal version ", thisVersion , ".\n",
  #              "You can download the more up to date\n",
  #              "version ", currentVersion , " at www.botany.ut/dispersal\n\n\n", sep="")
  # } else{
  #   ver = cat("You are using dispeRsal version ", thisVersion, ".\n\n", sep="")
  # }
  out <- list(p, TPL_nomatch)
  names(out) <- c("predictions", "unmatched_species")
  if(write.result==T) {
    write.table(out[[1]] , "predictedDD.txt")
    write.table(out[[2]] , "unmatched.txt")
  }
  out
}

## prep the data 
data <- all_clean 
data <- rename(data, "Species" = Taxon)

## required format of data:
## each row = 1 species 

## trait info:
## 1. dispersal syndrome DS
## animal, ant, ballistic, wind.none, wind.special
## 2. growth form GF
## tree, shrub, herb 
## 3. releasing height RH
## log10 transformed in m
## 4. seed mass M
## log10 transformed in mg 
## 5. terminal velocity TV
## log10 transformed in m/s

## transform terminal velocity
data$TV <- log10(data$TV)

## split into data for each different model 
## 1: dispersal distance = DS + GF + TV
data1 <- filter(data, !is.na(DS) & !is.na(GF) & !is.na(TV))

## 5: dispersal distance = DS + GF
data5 <- filter(data, !is.na(DS) & !is.na(GF) & is.na(TV))

## and predict!!
preds1 <- dispeRsal(data1, 
          model = 1, 
          CI = TRUE, 
          random = TRUE, 
          tax = "family", 
          write.result = FALSE)
## check for unmatched species
length(preds1[[2]]) ## 0
preds1 <- preds1[[1]]

preds1_join <- left_join(preds1, data1)

preds1_join <- gather(preds1_join, 
                     key = "CI_type", value = "CI", 
                     c("log10MDD_Family_lwrCL", "log10MDD_Family_uppCL",
                       "log10MDD_Order_lwrCL", "log10MDD_Order_uppCL",
                       "log10MDD_lwrCL", "log10MDD_uppCL")) %>%
  gather(key = "Pred_type", value = "Pred", 
         c("log10MDD_Family", "log10MDD_Order", "log10MDD")) %>%
  filter(!is.na(Pred)) %>%
  mutate(Pred_nonlog = 10^Pred) %>%
  mutate(log10MDD_measured_nonlog = 10^log10MDD_measured) %>%
  mutate(fixed_effects =  "DS + GF + TV")

## other model
preds5 <- dispeRsal(data5, 
                    model = 5, 
                    CI = TRUE, 
                    random = TRUE, 
                    tax = "family", 
                    write.result = FALSE)
## check for unmatched species
length(preds5[[2]]) ## 0
preds5 <- preds5[[1]]

preds5_join <- left_join(preds5, data5)

preds5_join <- gather(preds5_join, 
                      key = "CI_type", value = "CI", 
                      c("log10MDD_Family_lwrCL", "log10MDD_Family_uppCL",
                        "log10MDD_Order_lwrCL", "log10MDD_Order_uppCL",
                        "log10MDD_lwrCL", "log10MDD_uppCL")) %>%
  gather(key = "Pred_type", value = "Pred", 
         c("log10MDD_Family", "log10MDD_Order", "log10MDD")) %>%
  filter(!is.na(Pred)) %>%
  mutate(Pred_nonlog = 10^Pred) %>%
  mutate(log10MDD_measured_nonlog = 10^log10MDD_measured) %>%
  mutate(fixed_effects =  "DS + GF")

## plot predictions
preds1_join %>%
  ggplot(aes(x = Pred_nonlog, y = log10MDD_measured_nonlog, colour = DS,
             shape = GF, size = TV)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_point() +
  theme_bw() +
  geom_abline(intercept = 0, slope= 1) +
  facet_grid(~Pred_type) +
  labs(y = "Predicted dispersal distance (m)",x = "Measured dispersal distance (m)", 
       shape = "Growth form", size = "Log10(terminal velocity)", colour = "Dispersal syndrome")

preds1_join %>%
  ggplot(aes(x = Pred_nonlog, fill = Pred_type)) +
  geom_histogram() +
  scale_x_log10() +
  theme_bw() +
  facet_grid(~Pred_type) +
  labs(x = "Predicted dispersal distance (m)", y = "Frequency") +
  theme(legend.position = "none")

preds5_join %>%
  ggplot(aes(x = Pred_nonlog, y = log10MDD_measured_nonlog, colour = DS,
             shape = GF)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_point() +
  theme_bw() +
  geom_abline(intercept = 0, slope= 1) +
  facet_grid(~Pred_type) +
  labs(y = "Predicted dispersal distance (m)",x = "Measured dispersal distance (m)", 
       shape = "Growth form", colour = "Dispersal syndrome")

preds5_join %>%
  ggplot(aes(x = Pred_nonlog, fill = Pred_type)) +
  geom_histogram() +
  scale_x_log10() +
  theme_bw() +
  facet_grid(~Pred_type) +
  labs(x = "Predicted dispersal distance (m)", y = "Frequency") +
  theme(legend.position = "none")


## combine datasets 
predictions <- rbind(preds1_join, preds5_join)

predictions %>% ggplot(aes(x = Pred_nonlog, fill = fixed_effects)) +
  geom_histogram() +
  scale_x_log10() +
  theme_bw() +
  facet_grid(~Pred_type) +
  labs(x = "Predicted dispersal distance (m)", y = "Frequency") +
  theme(legend.position = "none")

preds1_join <- left_join(preds1, data1) %>%
  mutate(fixed_effects =  "DS + GF + TV")
preds5_join <- left_join(preds5, data5)%>%
  mutate(fixed_effects =  "DS + GF")
predictions <- rbind(preds1_join, preds5_join)

## write out 
write.csv(predictions, "data-processed/predicting-dispersal-distance/predictions_plants.csv", row.names = F)

##############################################################
####  calculate predicted potential range expansion rate  ####
##############################################################
## for species with predicted dispersal distance, get lifespan / age at maturity 
#---------------------
# TRY
#---------------------
try_am <- read.csv("data-processed/TRY_age-at-maturity_harmonized.csv")

cols_to_keep <- c("reported_name","reported_name_fixed", "scientificName", "kingdom", "phylum",
                  "class", "order", "family", "db", "db_code")

## subset to bioshifts species with dispersal distance
try_am_pred <- filter(try_am, scientificName %in% predictions$Species)
length(unique(try_am_pred$scientificName)) # 830 species 
length(unique(predictions$Species)) # out of the 1689 spp 

## clean the data 
unique(try_am_pred$OriglName)


try_am_pred = try_am_pred %>%
  filter(!OriglName %in% c("flowering", "Secondary juvenile period")) %>%
  mutate(Unit = ifelse(str_detect(OrigValueStr, "yrs") | OrigUnitStr %in% c("years?", "yr", "years", "year"), 
                       "yrs",
                       ifelse(str_detect(OrigUnitStr, "day"), 
                              "days", 
                              OrigUnitStr)))
try_am_pred$Unit[which(str_detect(try_am_pred$OrigValueStr, "weeks"))] = "weeks"

unique(try_am_pred$OrigUnitStr)
unique(try_am_pred$OrigValueStr)
unique(try_am_pred$Unit)

try_am_pred$Unit[which(is.na(try_am_pred$Unit))] <- "yrs"

## want minimum age at maturity to give a maximum dispersal estimate 
try_am_pred = try_am_pred %>%
  mutate(AgeAtMaturity = str_replace_all(OrigValueStr,"\\<", ""),
         AgeAtMaturity = str_replace_all(AgeAtMaturity,"\\>", ""),
         AgeAtMaturity = str_split_fixed(AgeAtMaturity, "\\-", 2)[,1]) %>%
  mutate(AgeAtMaturity = str_replace_all(AgeAtMaturity, "Over ", "")) %>%
  mutate(AgeAtMaturity = str_replace_all(AgeAtMaturity, "Within ", "")) %>%
  mutate(AgeAtMaturity = str_replace_all(AgeAtMaturity, "Between ", "")) %>%
  mutate(AgeAtMaturity = str_split_fixed(AgeAtMaturity, " ", 2)[,1]) %>%
  mutate(AgeAtMaturity = ifelse(is.na(OrigValueStr), NA, AgeAtMaturity))

unique(try_am_pred$AgeAtMaturity)

try_am_pred <- try_am_pred %>%
  select(all_of(cols_to_keep), 
         Unit, AgeAtMaturity, OriglName, 
         Reference) %>%
  mutate(Database = "TRY", Code = "AgeAtMaturity") %>%
  filter(!is.na(AgeAtMaturity)) %>%
  rename("Source" = Reference, "Field" = OriglName)

## save 
write.csv(try_am_pred, "data-processed/age-at-maturity-TRY-predictions.csv", row.names = FALSE)

## fill in the gaps:
missing_am = predictions %>%
  filter(!Species %in% try_am_pred$scientificName) %>%
  select(Species) %>%
  unique()

## using lifespan/growth form to infer age at maturity for annual plants
#---------------------
# lifespan compilation
#---------------------
## for species that live 1 year, age at maturity = 1 year 
## longevity: 
lifespan <- read.csv("data-raw/CompilationLifeSpan_02242022.csv")
unique(lifespan$Database)

## filter to species with dispersal distance
lifespan <- filter(lifespan, SpeciesChecked %in% predictions$Species)
length(unique(lifespan$SpeciesChecked)) ## 902 of our species 

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
                       Unit),
         Code = "MaturityFromLifespan") %>%
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
         AgeAtMaturity, Unit, Database, Field, Code) %>%
  unique() %>%
  filter(scientificName %in% missing_am$Species) %>%
  mutate(Source = "Lifespan compilation")

length(unique(lifespan_am$scientificName)) #58 more species woohoo!!

missing_am <- filter(missing_am, !Species %in% lifespan_am$scientificName)

#---------------------
# BROT
#---------------------
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
brot_sub <- filter(brot, scientificName %in% missing_am$Species) %>%
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
         Unit = "yrs", Database = "TRY", Code = "MaturityFromGrowthForm") %>%
  filter(!is.na(AgeAtMaturity)) %>%
  select(-Data)

length(unique(brot_am$scientificName))# 20 more

missing_am <- filter(missing_am, !Species %in% brot_am$scientificName)

## combine it all
all_am <- rbind(try_am_pred, lifespan_am) %>%
  rbind(., brot_am)

length(unique(all_am$scientificName)) ## 896 spp

## write 
write.csv(all_am, "data-processed/age-at-maturity-predictions.csv", row.names = FALSE)










###################
####   birds   ####
###################
## best model uses Kipp's distance, bill depth, tail graduation
# Kipp’s: 0.406
# Bill depth: -0.499
# Tail graduation: 0.004

y = 0.406*kd - 0.499*bd + 0.004*tg

## second best model uses Kipp's distance and bill depth
# Kipp’s: 0.368
# Bill depth: -0.402

y = 0.368*kd - 0.402*bd 

## natal dispersal distance 
## ln[ln(x+1) + 1] transformation to normalize geometric means of natal dispersal distances
## Kipp's distance and bill depth ln transformed 

## read in bioshifts
v1 <- read.csv("data-processed/corrected-bioshifts_fixed.csv")

## filter to passerines
v1 <- filter(v1, Order == "Passeriformes")
length(unique(v1$scientificName)) ## 831 spp

#########################
####  gather traits  ####
#########################
#---------------------
# AVOTNET
#---------------------
## Kipps.Distance
## Beak.Depth
## tail graduation not a measure 

avo <- read.csv("data-raw/primary-trait-data/AVONET/AVONET_raw.csv")

avo$genus_species <- avo$Species1

## harmonize names 
avo_harm <- harmonize(as.character(avo$genus_species))

notfound <- filter(avo_harm, is.na(db_code))

## rename columns 
avo <- avo %>%
  rename("reported_name" = genus_species) %>%
  mutate(reported_name_fixed = reported_name)

avo <- left_join(avo, avo_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## subset to bioshifts species 
avo_bs <- filter(avo, scientificName %in% v1$scientificName)
length(unique(avo_bs$Species1)) # 824 spp

avo_traits = avo_bs %>%
  mutate(Sex = ifelse(Female == 0, "M", 
                      ifelse(Male == 0, "F", 
                             "M/F"))) %>%
  select(all_of(cols_to_keep), Kipps.Distance, Beak.Depth) %>%
  rename("KippsDistance" = Kipps.Distance, 
         "BeakDepth" = Beak.Depth) %>%
  mutate(Database = "AVONET", 
         Unit = "mm", 
         Source = NA) 


############################
####  make predictions  ####
############################
## ln transform variables 
avo_traits$KippsDistanceLn = log(avo_traits$KippsDistance)
avo_traits$BeakDepthLn = log(avo_traits$BeakDepth)

## predict:
predictions <- 0.368*avo_traits$KippsDistanceLn - 0.402*avo_traits$BeakDepthLn
hist(predictions)

## back transform predictions 
## ln[ln(predictons+1) + 1]
back_trans <- exp(1)^(exp(1)^predictions - 1) - 1
hist(back_trans)
## now units are km

## problem: this model is for geometric mean natal dispersal distance
## we want arithmetic mean since closer to max
## solution: create our own model using their database 
## plus the database of Chu et al.

## recreate their models 
ecomorph <- read.csv("data-raw/dispersal/Dawideit_et_al_ecomorph.csv")
ecomorph[ecomorph == "-"] <- NA

paradis <- read.csv("data-raw/dispersal/Paradis_et_al_2002.csv")

## fix names 
paradis[paradis == "Delichon urbica"] <- "Delichon urbicum"
paradis[paradis == "Parus montanus"] <- "Poecile montanus"
paradis[paradis == "Parus ater"] <- "Periparus ater"
paradis[paradis == "Parus caeruleus"] <- "Cyanistes caeruleus"

join <- inner_join(ecomorph, paradis)

## normalize traits 
join$Kipp.s.distanceLn <- log(join$Kipp.s.distance)
join$Bill.depthLn <- log(join$Bill.depth)

## transform geometric mean of dispersal distances 
## ln[ln(x + 1) + 1]
join$GM_natal_transformed <- log(log(join$GM_natal + 1)+1)

hist(join$GM_natal)
hist(join$GM_natal_transformed)

## fit the model
mod <- lm(GM_natal_transformed ~ Kipp.s.distanceLn + Bill.depthLn,
          data = join)
summary(mod)

## get the phylogenetic tree 
library(caper) 
data(BritishBirds)

BritishBirds.tree
join$Species <- str_replace_all(join$Species, " ",  ".")

## make sure all species are there 
length(which(join$Species %in% BritishBirds.tree$tip.label))
## yay

# corr = corPagel(value = 1, 
#                 phy = BritishBirds.tree, 
#                 form = ~Species,
#                 fixed = FALSE)
# 
# corr = Initialize(corr, data = join) # initialize object
# 
# ## extract the correlation matrix and inspect 
# matrix <- matrix(corMatrix(corr),
#                  nrow = length(unique(join$Species)))
# #View(matrix)
# 
# gls = gls(GM_natal_transformed ~ Kipp.s.distanceLn + Bill.depthLn,
#                
#                correlation = corr,
#                
#                data = join)
# summary(gls)
# 
#         
# row.names(matrix) <- BritishBirds.tree$tip.label[BritishBirds.tree$tip.label %in% join$Species]
# row.names(join) <- join$Species


## try phylogenetic comparative model using caper
## make comparative object
join_comp <- comparative.data(BritishBirds.tree, join, Species, vcv = TRUE)

## fit the GLS using the phylogenetic comparative method
mod_pgls <- pgls(GM_natal_transformed ~ Kipp.s.distanceLn + Bill.depthLn,
                 lambda = "ML",
                 data = join_comp)

summary(mod_pgls)

## better adjusted R2 than model without phylogeny









## okay, whatever
## let's see if I can build my own predictive model for arithmetic mean natal dispersal distance
## add data from Chu so it is even more predictive 

## think: do I want to control for phylogeny? means I will need to have a phylogeny for species I want to predict for 
## try with and without 


## transform arithmetic mean of dispersal distances 
## ln[ln(x + 1) + 1]
join$AM_natal_transformed <- log(log(join$AM_natal + 1)+1)

hist(join$AM_natal)
hist(join$AM_natal_transformed)

mod <- lm(AM_natal_transformed ~ Kipp.s.distanceLn + Bill.depthLn,
          data = join)
summary(mod)
## same sign of relationships, better R2

## fit the GLS using the phylogenetic comparative method
## make comparative object
join_comp <- comparative.data(BritishBirds.tree, join, Species, vcv = TRUE)

mod_pgls <- pgls(AM_natal_transformed ~ Kipp.s.distanceLn + Bill.depthLn,
                 lambda = "ML",
                 data = join_comp)

summary(mod_pgls)
## even better R2


## use this model
## ln transform variables 
avo_traits$KippsDistanceLn = log(avo_traits$KippsDistance)
avo_traits$BeakDepthLn = log(avo_traits$BeakDepth)

## predict using simple linear model:
predictions <- 0.3042*avo_traits$KippsDistanceLn - 0.3576*avo_traits$BeakDepthLn
hist(predictions)

## back transform predictions 
## ln[ln(predictons+1) + 1]
back_trans <- exp(1)^(exp(1)^predictions - 1) - 1
hist(back_trans)
## now units are km



## wait a second 
## many of the bird dispersal observations are NOT MEAN natal dispersal (ex. ones from sutherland are max / median, some are breeding dispersal)
## need to filter real data to only plants + geometric/arithmetic mean natal dispersal and see if pattern holds 
## if not, no hope in predicting 


## make many models using Paradis ?
## ex. one for natal, one for breeding 

## need to see if results are sensitive to choosing only observations that are MAX / only observations that are MEAN


