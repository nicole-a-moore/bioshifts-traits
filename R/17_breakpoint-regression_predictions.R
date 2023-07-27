## applying breakpoint regression using predicted dispersal distances 
library(tidyverse)
library(segmented)

select <- dplyr::select

## read function to harmonize taxonomy 
source("R/taxonomic-harmonization/clean_taxa_functions.R")
source("R/harmonize.R")


###########################################################
####   range shift and dispersal distance preparation  ####
###########################################################
#----------------------
## read in prediced dispersal scale data 
dscale = read.csv("data-processed/predicting-dispersal-distance/predictions_plants.csv") 

## choose one dispersal distance prediction per species x dispersal syndrome:
dscale <- dscale %>%
  mutate(DispersalDistance = ifelse(fixed_effects == "DS + GF + TV", log10MDD,
                                    ifelse(fixed_effects ==  "DS + GF" & !is.na(log10MDD_Family),
                                           log10MDD_Family,
                                           ifelse(fixed_effects ==  "DS + GF" & !is.na(log10MDD_Order),
                                                  log10MDD_Order, 
                                                  log10MDD)))) %>%
  mutate(DispersalDistanceModel = ifelse(fixed_effects == "DS + GF + TV", "DS + GF + TV",
                                    ifelse(fixed_effects ==  "DS + GF" & !is.na(log10MDD_Family),
                                           "DS + GF + (1|Order/Family)",
                                           ifelse(fixed_effects ==  "DS + GF" & !is.na(log10MDD_Order),
                                                  "DS + GF + (1|Order)",
                                                  "DS + GF")))) %>%
  select(Species, DS, DispersalDistance, DispersalDistanceModel) %>%
  distinct()


## convert to Km
dscale <- dscale %>%
  mutate(DispersalDistance = 10^DispersalDistance) %>% ## unlog
  mutate(DispersalDistanceKm = DispersalDistance/1000) %>%
  filter(!is.na(DispersalDistanceKm)) 

#----------------------
## read in list of all bioshifts species 
sp <- read_csv("data-raw/splist.csv")

spv1 <- filter(sp, v1 == 1) %>%
  select(scientificName, species) %>%
  unique()

#----------------------
## read in bioshifts v1 
v1 = read.csv("data-processed/corrected-bioshifts.csv")

## do some name fixes 
v1$Species[which(!v1$Species %in% spv1$species)]

v1$Species = Clean_Names(v1$Species, return_gen_sps = F)

v1$Species[which(v1$Species == "Quercus x")] = "Quercus" 
v1$Species[which(v1$Species == "Mentha x")] = "Mentha" 
v1$Species[which(v1$Species == "Circaea x intermedia")] = "Circaea intermedia" 

v1$Species = Clean_Names(v1$Species, return_gen_sps = F)

v1$Species[which(!v1$Species %in% spv1$species)] ## all are there 

## yay! all are there
which(!v1$Species %in% spv1$species)

v1 = left_join(v1, spv1, by = c("Species" = "species"))

## save 
write.csv(v1, "data-processed/corrected-bioshifts_fixed.csv", row.names = FALSE)

#transform study area
v1$Area <- log(v1$Area)

## subset to species with dispersal scale prediction 
v1 <- filter(v1, scientificName %in% dscale$Species)
length(unique(v1$Species)) #1707 species 

#----------------------
## add dispersal scale 
## get rid of old taxonomy columns from v1 (they aren't right)
v1 <- select(v1, -c("Kingdom", "Phylum", "Class", "Order", "Family"))
v1_saved = v1

## rename some columns 
v1 = left_join(v1, dscale, by = c("scientificName" = "Species")) 

## check on merge
length(which(is.na(v1$DispersalDistanceKm))) #0 missing dispersal scale
length(unique(v1$Species)) #still 1707 species


## save dataset
# write.csv(v1, "data-processed/corrected-bioshiftsv1_dispersal-distance-predicted.csv", row.names = FALSE)
v1 = read.csv("data-processed/corrected-bioshiftsv1_dispersal-distance-predicted.csv")

###################################################
####   calculating potential dispersal rate    ####
###################################################
### join age at maturity data with dispersal data 
am <- read.csv("data-processed/age-at-maturity-predictions.csv")

## if multiple estimates of age at maturity per species, keep the lowest 
am_join <- am %>%
  group_by(scientificName) %>%
  mutate(AgeAtMaturity = as.numeric(as.character(AgeAtMaturity))) %>%
  mutate(AgeAtMaturityDays = ifelse(Unit == "yrs", 
                                    AgeAtMaturity*365,
                                    ifelse(Unit == "weeks",
                                           AgeAtMaturity*7,
                                           AgeAtMaturity))) %>% # convert all to days 
  mutate(AgeAtMaturityDays = min(AgeAtMaturityDays)) %>% # select minimum per species 
  ungroup() %>%
  select(scientificName, AgeAtMaturityDays) %>%
  unique() %>%
  mutate(YearOfMaturity = ceiling(AgeAtMaturityDays/365)) ## make new column for a value that's rounded to the nearest year 

## join to dispersal data:
v1 <- left_join(v1, am_join, by = "scientificName")

## get rid of duplicated data (some are same species but have different scientific names)
v1 <- select(v1, -scientificName) %>%
  distinct()

length(unique(v1$Species)) ## still have all the species
length(unique(v1$Species[which(is.na(v1$AgeAtMaturityDays))])) 
## 802 / 1707 species do not have age at maturity data 

unique(v1$Species[which(is.na(v1$AgeAtMaturityDays))])

## filter to only species with age at maturity 
v1 <- filter(v1, !is.na(AgeAtMaturityDays))

## calculate dispersal potential for species with age at maturity 
v1 = v1 %>%
  mutate(DispersalPotentialKmY = ifelse(!is.na(YearOfMaturity), 
                                        DispersalDistanceKm/YearOfMaturity,
                                        NA)) %>%
  mutate(DispersalPotentialmY = ifelse(!is.na(YearOfMaturity), 
                                       (DispersalDistanceKm*1000)/YearOfMaturity,
                                       NA)) %>%
  mutate(DispersalDistancem = DispersalDistanceKm*1000)


############################################
####     subset data, convert units     ####
############################################
## note: duplicate observations for species with multiple dispersal syndromes 

mod_data <- v1 %>%
  ## filter out observations where climate velocity is negative at leading edge/optimum (expect contraction)
  filter(LatVeloT >= 0 | EleVeloT >= 0) %>%
  ## get rid of negative shifts 
  filter(ShiftR > 0) %>%
  ## get rid of trailing edge 
  filter(., Position != "Trailing edge") %>%
  ## make sure none have empty dispersal potential 
  filter(!is.na(DispersalPotentialKmY)) %>%
  ## make one column for potential dispersal rate, climate velo for easier plotting of lat x elev data together
  mutate(annual_dispersal_pot = ifelse(Gradient == "Elevation",
                                       DispersalPotentialmY,
                                       ifelse(Gradient == "Latitudinal",
                                              DispersalPotentialKmY,
                                              NA))) 

## convert units of latitudinal and elevation shifts & dispersal & climate velocity to km/y:
mod_data$ShiftKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$ShiftR / 1000, mod_data$ShiftR)
mod_data$ClimVeloTKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$EleVeloT / 1000, mod_data$LatVeloT)
mod_data$AnnualDispPotKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$annual_dispersal_pot / 1000, mod_data$annual_dispersal_pot)


mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)

mod_data <- select(mod_data, -c(SLDiff, CorrShift, PredSLShift, study_level_shift)) %>%
  distinct()

## make plots to summarize subset 
nrow(mod_data) # 3276 range shifts (but this includes duplicate dispersal syndromes)
length(unique(mod_data$Species)) # 803 species 

mod_data %>%
  filter(AnnualDispPotKmY < ClimVeloTKmY) %>%
  filter(ShiftR > 0.0001) %>%
  filter(AnnualDispPotKmY < 50) %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) +
  theme_bw() +
  geom_point() +
  labs(x = "Potential dispersal rate (km/y)", y = "Observed range expansion rate (km/y)",
       colour = "") +
  scale_y_continuous(limits = c(0, 10)) +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  theme(panel.grid = element_blank()) +
  scale_x_log10() 

mod_data %>%
 # filter(AnnualDispPotKmY > ClimVeloTKmY) %>%
  filter(ShiftKmY > 0.0001) %>%
 # filter(AnnualDispPotKmY < 50) %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) +
  theme_bw() +
  geom_point() +
  labs(x = "Potential dispersal rate (km/y)", y = "Observed range expansion rate (km/y)",
       colour = "") +
  scale_y_continuous(limits = c(0, 10)) +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  theme(panel.grid = element_blank()) +
  scale_x_log10() 




########################################################
####      fitting global breakpoint regression      ####
########################################################
## filter to shifts > 0.0001
mod_data <- mod_data %>%
  filter(ShiftR > 0.0001) 

## split between elev and lat
lat = filter(mod_data, Gradient == "Latitudinal")
ele = filter(mod_data, Gradient == "Elevation")









