## using other plant traits to infer dispersal scale
library(tidyverse)
library(plyr)

## dispeRsal
## https://figshare.com/collections/Predicting_species_maximum_dispersal_distances_from_simple_plant_traits/3306522

## format of data
## each row = 1 species 
## at minimum, each species must have dispersal syndrome

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

## let's gather these traits!!!

## read in bioshifts species list
sp <- read_csv("data-raw/splist.csv")

## filter to plants 
sp <- sp %>% filter(kingdom == "Plantae")

## BIEN
library(BIEN)

## download all trait data for each bioshift species
bien = BIEN_trait_species(unique(sp$scientificName))

# write.csv(bien, "data-processed/BIEN_bioshifts.csv", row.names = FALSE)
bien = read.csv("data-processed/BIEN_bioshifts.csv")
unique(bien$trait_name)

## GF
bien_gf <- bien %>%
  filter(trait_name == "whole plant growth form") %>%
  mutate(trait_value = str_remove_all(trait_value, "\\*")) %>%
  mutate(trait_value = str_remove_all(trait_value, "\\-")) 

## reconcile species with multiple growth forms 
bien_multi = bien_gf %>%
  filter(scrubbed_species_binomial %in% 
           bien_gf$scrubbed_species_binomial[which(duplicated(bien_gf$scrubbed_species_binomial))]) %>%
  group_by(scrubbed_species_binomial) %>%
  filter(str_detect(method, "Species with multiple growth form classifications")) %>%
  filter(!duplicated(scrubbed_species_binomial)) # remove ones that are still duplicated

bien_gf = bien_gf %>%
  filter(!scrubbed_species_binomial %in% 
           bien_gf$scrubbed_species_binomial[which(duplicated(bien_gf$scrubbed_species_binomial))]) %>%
  rbind(., bien_multi) %>%
  mutate(GF = ifelse(trait_value %in% c("tree", "Tree", "Small_Tree", "tree ", "giant tree"), "tree", 
                     ifelse(trait_value %in% c("Herb", "herb", "Graminoid", "geophyte", "forb", "Scandent_Herb",
                                      "climbing herb", "twining herb", "herbaceous climber",
                                      "Grass", "grass"), "herb", 
                            ifelse(trait_value %in% c("Shrub", "subshrub", "shrub", "shrublet"), "shrub", NA)))) %>%
  select(scrubbed_species_binomial, GF) %>%
  filter(!is.na(GF)) %>%
  unique() %>%
  rename("Taxon" = scrubbed_species_binomial)

unique(bien_gf$GF)

## "seed mass" 
bien_sm <- bien %>%
  filter(trait_name == "seed mass") 

## make sure we do not lost data in numeric conversion
val = bien_sm$trait_value
new_val = as.numeric(as.character(bien_sm$trait_value))
val[which(is.na(new_val))]

bien_sm <- bien_sm %>%
  mutate(SM = as.numeric(as.character(trait_value))) %>%
  rename("Taxon" = scrubbed_species_binomial) %>%
  select(Taxon, SM) %>%
  filter(!is.na(SM)) %>%
  unique() 

bien_final = full_join(bien_gf, bien_sm)

## filter to bioshifts species
bien_sub <- filter(bien_final, Taxon %in% sp$scientificName) %>%
  mutate(source = "BIEN")


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

sm = brot %>%
  filter(Trait == "SeedMass") %>%
  mutate(SM = Data) %>%
  select(Taxon, SM) %>%
  unique()
  
  
## make data so that is has 1 row per species per unique trait combination
brot_final <- dm %>%
  select(Taxon, DS) %>%
  full_join(select(gf, c(Taxon, GF))) %>%
  full_join(sm)

## filter to bioshifts species
brot_sub <- filter(brot_final, Taxon %in% sp$scientificName)%>%
  mutate(source = "BROT")

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
  mutate(SM = dia_mass) %>%
  mutate(Taxon = paste(str_split_fixed(name, " ", 2)[,1], str_split_fixed(name, " ", 3)[,2], sep = " ")) %>%
  select(Taxon, DS, TV, SM) %>%
  filter(!is.na(DS) | !is.na(TV) | !is.na(SM))

## filter to bioshifts species
d3_sub <- filter(d3_final, Taxon %in% sp$scientificName) %>%
  mutate(source = "D3")


## LEDA
## Growth form
ledagf = read.delim("data-raw/primary-trait-data/LEDA/plant_growth_form copy.txt", sep = ";")

# https://hosho.ees.hokudai.ac.jp/tsuyu/top/dct/lf.html
# plant life form classifications:

# Phanerophytes, epiphyte, sclerophyte - tree
# Macrophanerophytes - tall tree
# Mesophanerophytes - short tree or shrub 
# Nanophanerophytes - shrub
# Chamaephytes [Ch]  - shrub
# Hemicryptophytes, Geophytes - perennial herb
# Helophytes - underwater perennial herb
# Therophytes - annual herb
# forb, graminoid, geophyte - herb
# Liana - classify as tree

unique(ledagf$plant.growth.form)
## exclude parasites since likely not able to accurately predict dispersal distance

## clean data so everything is tree shrub or herb
ledagf_final = ledagf %>%
  mutate(GF = ifelse(plant.growth.form %in% c("Chamaephyte"), "shrub", 
                     ifelse(plant.growth.form %in% c("Geophyte", "Hemicryptophyte", 
                                                     "Therophyte"), "herb", 
                            ifelse(plant.growth.form %in% c("Phanerophyte"), "tree", NA)))) %>%
  select(SBS.name, GF) %>%
  rename("Taxon" = SBS.name) %>%
  filter(!is.na(GF)) %>%
  unique()
  
length(unique(ledagf_final$ledagf_final))

## get rid of ones with two growth form classifications for now
ledagf_final <- filter(ledagf_final, ! Taxon %in% ledagf_final$Taxon[which(duplicated(ledagf_final$Taxon))])

## filter to bioshifts species
ledagf_sub <- filter(ledagf_final, Taxon %in% sp$scientificName)

## Terminal velocity
ledatv = read.delim("data-raw/primary-trait-data/LEDA/TV_2016 copy.txt", sep = ";")

## get rid of mean observations (sample size != 1) 
ledatv_final = ledatv %>%
  filter(!general.method == "unknown") %>%
  filter(sample.size == 1) %>%
  mutate(TV = single.value..m.s.) %>%
  rename("Taxon" = SBS.name) %>%
  select(Taxon, TV) %>%
  unique() 
 
## filter to bioshifts species
ledatv_sub <- filter(ledatv_final, Taxon %in% sp$scientificName) 

## Dispersal type
ledadt = read.delim("data-raw/primary-trait-data/LEDA/dispersal_type copy.txt", sep = ";")

ledadt = filter(ledadt, dispersal.type != "") # get rid of empty data

unique(ledadt$dispersal.type)

# https://d-nb.info/1122297556/34
# https://d-nb.info/975149946/34

## animal = c("epizoochor","endozoochor", "dysochor", "zoochor")
## wind.none = 
## wind.special = 
## ant = 
## ballistic = c("autochor", "ballochor", "blastochor", "herpochor")

## wind unknown = "boleochor", "meteorochor"

# "hemerochor" = "speirochor" = "agochor" = human
# "ombrochor" = "bythisochor" = "nautochor" = water
# "chamaechor" = tumble on ground in wind

ledadt_final <- ledadt %>% 
  mutate(DS = ifelse(dispersal.type == "dysochor" & str_detect(.$dispersal.vector, "ant"), 
                     "ant",
                     ifelse(dispersal.type %in% c("epizoochor","endozoochor", "dysochor", "zoochor"),
                            "animal",
                            ifelse(dispersal.type %in% c("autochor", "ballochor", "blastochor", "herpochor"),
                                   "ballistic", 
                                   NA)))) %>%
  rename("Taxon" = SBS.name) %>%
  select(Taxon, DS) %>%
  filter(!is.na(DS)) %>%
  unique() 
                            
length(unique(ledadt_final$Taxon))         

## filter to bioshifts species
ledadt_sub <- filter(ledadt_final, Taxon %in% sp$scientificName)

## Release height
ledarh = read.delim("data-raw/primary-trait-data/LEDA/releasing_height copy.txt", sep = ";")

## calculate a mean rh per species 
## get rid of ones with unknown methodology 
## single value = max - min for a study
## weight mean by sample size
ledarh_final = ledarh %>%
  filter(!general.method == "unknown") %>%
  mutate(RHsample.size = ifelse(is.na(sample.size), 1, sample.size)) %>%
  mutate(RH = single.value..m., sample.size) %>%
  rename("Taxon" = SBS.name) %>%
  select(Taxon, RH, RHsample.size) %>%
  filter(!is.na(RH)) %>%
  unique() 

## filter to bioshifts species
ledarh_sub <- filter(ledarh_final, Taxon %in% sp$scientificName)

## seed mass
ledasm = read.delim("data-raw/primary-trait-data/LEDA/seed_mass copy.txt", sep = ";")
## leave for now, needs cleaning


## combine all leda subsets 
leda_sub = full_join(ledarh_sub, ledatv_sub) %>%
  full_join(., ledadt_sub) %>%
  full_join(., ledagf_sub) %>%
  mutate(source = "LEDA")

## TRY
## read in TRY query
try = read_delim("data-raw/primary-trait-data/TRY/23454.txt")
unique(try$TraitName)

## subset to our bioshift species bc this data is HUGE
try <- filter(try, SpeciesName %in% sp$scientificName)

## GF
try_gf <- try %>%
  filter(TraitName == "Plant growth form")

length(unique(try_gf$SpeciesName)) ## assuming we can trust/clean every data, potential to get GF for 5800 bioshift species

unique(try_gf$OrigValueStr) ## oh god
unique(try_gf$OriglName) ## oh god x2

## get rid of original names that sound like not what we want 
try_gf <- try_gf %>%
  filter(!OriglName %in% c("shoot growth form", "ClimbingMode", "Multistemness/ Growth form ", 
  "Multistemness", "Succulence", "Succulence of Leaves or Stem", "Funct. group", "Leaf succulence",
  "Stem succulent", "Plant size", "succulent species", "Succulence index", "CONSENSUS", "Parasite",
  "carnivory", "aquatic", "climber", "Parasitic", "Transect"))

unique(try_gf$OriglName)

## low-hanging fruit
try_gf <- try_gf %>%
  filter(!str_detect(OrigValueStr, "aquatic"), !str_detect(OrigValueStr, "Aquatic")) %>% ## get rid of aquatic
  filter(!str_detect(OrigValueStr, "moss"), !str_detect(OrigValueStr, "Moss"), !str_detect(OrigValueStr, "MOSS")) %>% ## get rid of moss
  filter(!str_detect(OrigValueStr, "PARASITE"), !str_detect(OrigValueStr, "parasite"),
         !str_detect(OrigValueStr, "Parasite"),!str_detect(OrigValueStr, "carniv"),
         !str_detect(OrigValueStr, "Succulent"),  !str_detect(OrigValueStr, "non-succulent"),
         !str_detect(OrigValueStr, "Free-standing"), !str_detect(OrigValueStr, "No"),
         !str_detect(OrigValueStr, "succulent")) %>% ## get rid of parasites, succulent
  filter(OrigValueStr != "?") %>% ## get rid of ?
  mutate(GF = OrigValueStr) %>%
  mutate(GF = ifelse(GF %in% c("shrub", "Shrub", "Shurb", "subshrub", "Shrub, Subshrub",
                                         "Subshrub, Shrub", "Shrub/Subshrub", "Shrub (S)", 
                               "Subshurb", "sub-shrub", "Sub-Shrub (Chamaephyte)", "shrubs", "shrub or chamaephyt", 
                               "S", "SS", "sh"),
                     "shrub", GF)) %>%
  mutate(GF = ifelse(GF %in% c("tree", "Tree","Tree (T)", "t", "T", "Tree (deciduous)", "Tree (evergreen)", 
                               "TREE", "trees", "trees/Tree", "trees/tree/Tree", "trees/tree", "Tree/Treelet", 
                               "Small_Tree", "small tree", "tree (woody >4m)", "tree/woody", "epiphyte", "epiphytes", 
                               "Epiphyte", "Epiphyte/Tree", "T  resp.  T"),
                     "tree", GF)) %>%
  mutate(GF = ifelse(GF %in% c("herb", "Herb", "herbs", "Herbaceous", "Herbaceous Monocot", "Herbaceous Dicot", 
                               "Herbaceous Forb", 
                               "H", "herb.", "herbaceous legume", "h", "herbaceous", "herb/non-woody", 
                               "herbaceous/Terrestrial Herb", "Trailing_Herb", "Forb/herb", "perennial", "annual",
                               "Terrestrial Herb", "herbaceous plant", "herbaceous perennial", "herbaceous annual-biennial",
                               "perennial legume", "herbaceous dicotyl", "herbaceous monocotyl", "perennial herb", 
                               "forb (herbaceous, with or without woody base)","perennial leguminous herb", 
                               "herbaceous/?", "geophyte", "Geophyte", "graminoid", "perennial graminoid", 
                               "graminoid/non-woody", "Perennial graminoid", "Annual graminoid", "Graminoid", "Graminoids",
                               "therophyte", "forb",  "leguminous forb", "Annual forb", "Perennial forb", "perennial forb", 
                               "annual forb", "Forb", "Forbs", "Herbaceous dicots", "f", "Annual Herb", "F"),
                     "herb", GF))

## see how many per species 
try_gf %>%
  filter(GF %in% c("tree", "shrub", "herb")) %>%
  group_by(SpeciesName, GF) %>%
  tally() %>%
  arrange(n) %>% View
  
## good enough for now - do a better job of cleaning later 
## let's say if more than 10 records call something a tree/herb/shrub, it's safe to call it that 
classif = try_gf %>%
  filter(GF %in% c("tree", "shrub", "herb")) %>%
  group_by(SpeciesName, GF) %>%
  tally() %>%
  filter(n > 10) %>%
  arrange(n)

length(which(duplicated(classif$SpeciesName))) ## 71 have 2 classifications 

## exclude them for now 
classif <- filter(classif, !SpeciesName %in% classif$SpeciesName[which(duplicated(classif$SpeciesName))])

length(unique(classif$SpeciesName)) ## growth form for 1897 bioshift species for now 

## for now, leave this alone
## move on to dispersal syndrome


## DS
try_ds <- try %>%
  filter(TraitName == "Dispersal syndrome")

length(unique(try_ds$SpeciesName)) ## potential for 2247 bioshifts species 

unique(try_ds$OrigValueStr)


## TV
try_tv <- try %>%
  filter(TraitName == "Seed terminal velocity")

length(unique(try_tv$SpeciesName)) ## potential for 2247 bioshifts species 

unique(try_tv$OrigValueStr)

## RH
try_rh <- try %>%
  filter(TraitName == "Seed releasing height")

length(unique(try_rh$SpeciesName)) ## potential for 2247 bioshifts species 

unique(try_rh$OrigValueStr)


## Tamme 
tamme <- read_csv("data-raw/dispersal/Tamme_DispersalDistanceData.csv")
colnames(tamme) <- str_replace_all(colnames(tamme), "\\ ", "_")

## reorganize
tamme_final <- tamme %>%
  rename("Taxon" = Species, "DS" = Dispersal_syndrome,
         "GF" = Growth_form, "TV" = `Seed_terminal_velocity_(m/s)`,
         "RH" = `Seed_release_height_(m)`,
         "SM" = `Seed_weight_(mg)`) %>%
  select(Taxon, DS, GF, RH, TV, SM) %>%
  mutate(DS = ifelse(DS == "wind (special)", "wind.special",
                     ifelse(DS == "wind (none)", "wind.none", 
                            DS)))

## filter to bioshifts species
tamme_sub <- filter(tamme_final, Taxon %in% sp$scientificName) %>%
  mutate(source = "Tamme")

## combine them all!!!!!
all = rbind.fill(tamme_sub, leda_sub) %>%
  rbind.fill(., brot_sub) %>%
  rbind.fill(., d3_sub) %>%
  rbind.fill(., bien_sub)

## need to resolve multiple values per species 
all %>%
  group_by(Taxon) %>%
  tally()






## load dispeRsal objects into workspace 
load("R/dispeRsal.rda")

## fix issue with function 
dispeRsal = function(data.predict, model, CI = F, random = T, tax = "family",
                     write.result=F){
  require(nlme)
  require(AICcmodavg)
  thisVersion <- 0.2
  currentVersion <- 0.2
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
    DS <- pd[na.omit(pmatch(p$Species, pd$Species)),]$DS
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
  
  if (thisVersion != currentVersion) {
    ver = cat( "\n\n\nATTENTION\n",
               "You are using dispeRsal version ", thisVersion , ".\n",
               "You can download the more up to date\n",
               "version ", currentVersion , " at www.botany.ut/dispersal\n\n\n", sep="")
  } else{
    ver = cat("You are using dispeRsal version ", thisVersion, ".\n\n", sep="")
  }
  out <- list(p, TPL_nomatch)
  names(out) <- c("predictions", "unmatched_species")
  if(write.result==T) {
    write.table(out[[1]] , "predictedDD.txt")
    write.table(out[[2]] , "unmatched.txt")
  }
  out
}

## read in practice data
own.data

predictions = dispeRsal(data.predict = own.data, 
              model = 1, 
              CI = TRUE, 
              random = TRUE, 
              tax = "family", 
              write.result = FALSE)

predictions_simple = dispeRsal(data.predict = own.data, 
                        model = 1, 
                        CI = TRUE, 
                        random = FALSE, 
                        tax = "family", 
                        write.result = FALSE)

df = predictions[[1]]
df_simple = predictions_simple[[1]]




