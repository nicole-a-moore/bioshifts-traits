## using other plant traits to infer dispersal scale
library(tidyverse)

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
bien_final <- bien %>%
  filter(trait_name == "whole plant growth form") %>%
  mutate(trait_value = str_remove_all(trait_value, "\\*")) %>%
  mutate(trait_value = str_remove_all(trait_value, "\\-")) 

## reconcile species with multiple growth forms 
bien_multi = bien_final %>%
  filter(scrubbed_species_binomial %in% 
           bien_final$scrubbed_species_binomial[which(duplicated(bien_final$scrubbed_species_binomial))]) %>%
  group_by(scrubbed_species_binomial) %>%
  filter(str_detect(method, "Species with multiple growth form classifications")) %>%
  filter(!duplicated(scrubbed_species_binomial)) # remove ones that are still duplicated

bien_final = bien_final %>%
  filter(!scrubbed_species_binomial %in% 
           bien_final$scrubbed_species_binomial[which(duplicated(bien_final$scrubbed_species_binomial))]) %>%
  rbind(., bien_multi) %>%
  mutate(GF = trait_value,
         GF = ifelse(GF %in% c("tree", "Tree", "Small_Tree", "tree ", "giant tree"), "tree", 
                     ifelse(GF %in% c("Herb", "herb", "Graminoid", "geophyte", "forb", "Scandent_Herb",
                                      "climbing herb", "twining herb", "herbaceous climber",
                                      "Grass", "grass"), "herb", 
                            ifelse(GF %in% c("Shrub", "subshrub", "shrub", "shrublet"), "shrub", NA)))) %>%
  select(scrubbed_species_binomial, GF) %>%
  filter(!is.na(GF)) %>%
  unique() %>%
  rename("Taxon" = scrubbed_species_binomial)

unique(bien_final$GF)

## filter to bioshifts species
bien_sub <- filter(bien_final, Taxon %in% sp$scientificName)


## BROT
brot = read.csv("data-raw/primary-trait-data/BROT/BROT2_dat.csv")
unique(brot$Trait)
## GrowthForm, DispMode
brot = brot %>%
  filter(Trait %in% c("GrowthForm", "DispMode"))

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
  inner_join(select(gf, c(Taxon, GF)))

## filter to bioshifts species
brot_sub <- filter(brot_final, Taxon %in% sp$scientificName)

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
  mutate(TV = log10(vterm)) %>%
  mutate(Taxon = paste(str_split_fixed(name, " ", 2)[,1], str_split_fixed(name, " ", 3)[,2], sep = " ")) %>%
  select(Taxon, DS, TV) %>%
  filter(!is.na(DS) | !is.na(TV))

## filter to bioshifts species
d3_sub <- filter(d3_final, Taxon %in% sp$scientificName)


## LEDA

## PalmTraits

## TRY

## Tamme 






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




