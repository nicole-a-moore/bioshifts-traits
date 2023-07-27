## applying breakpoint regression 
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
## read in dispersal scale data 
dscale = read.csv("data-processed/dispersal-distance-collated.csv") 

unique(dscale$Unit)

val = dscale$DispersalDistance
new_val = as.numeric(as.character(dscale$DispersalDistance))
val[which(is.na(new_val))]
## only one that isn't a number is a range 
## select max for now

## convert all to Km
dscale <- dscale %>%
  mutate(DispersalDistance = ifelse(DispersalDistance == "10-40", "40", DispersalDistance)) %>%
  mutate(DispersalDistance = as.numeric(as.character(DispersalDistance))) %>%
  mutate(DispersalDistanceKm = ifelse(Unit == "m",
                                      DispersalDistance/1000, 
                                      DispersalDistance)) %>%
  filter(!is.na(DispersalDistanceKm)) 

## plot distribution
dscale %>%
  ggplot(aes(x = log(DispersalDistanceKm), fill = class)) + geom_histogram() +
  theme_bw()

dscale %>%
  ggplot(aes(x = Code, fill = class)) + geom_bar() +
  theme_bw() + coord_flip()


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

## subset to species with dispersal scale 
v1 <- filter(v1, scientificName %in% dscale$scientificName)
length(unique(v1$Species)) #696 species 


#----------------------
### relate velocity of shift to dispersal scale ###
colnames(v1)
# colnames(v2)


#----------------------
## add dispersal scale 
## get rid of old taxonomy columns from v1 (they aren't right)
v1 <- select(v1, -c("Kingdom", "Phylum", "Class", "Order", "Family"))
v1_saved = v1

## get rid of columns that will cause duplication in dispersal scale database 
dscale <- select(dscale, -c("reported_name", "reported_name_fixed", "db", "db_code")) %>%
  unique()

## rename some columns 
dscale <- rename(dscale, "DispersalSource"= Source, "DispersalUnit"= Unit)

v1 = left_join(v1, dscale, by = c("scientificName" = "scientificName")) 

## check on merge
length(which(is.na(v1$DispersalDistanceKm))) #0 missing dispersal scale
length(unique(v1$Species)) #still 696 species

## get rid of the single tunicate (probably passively disperses)
v1 <- filter(v1, scientificName != "Salix alba")
length(unique(v1$Species)) ## now have 695


## see which taxa have movement studies 
v1 %>%
  filter(ObservationTypeGeneral == "movement study") %>%
  ggplot(aes(x = DispersalDistanceKm, fill = class)) + geom_histogram() + scale_x_log10()
## and which do not
v1 %>%
  filter(ObservationTypeGeneral != "movement study") %>%
  ggplot(aes(x = DispersalDistanceKm, fill = class)) + geom_histogram() + scale_x_log10()

## get rid of movement studies for now:
v1 <- filter(v1, ObservationTypeGeneral != "movement study")
length(unique(v1$Species)) ## now have 609

## save dataset
# write.csv(v1, "data-processed/corrected-bioshiftsv1_max-dispersal-distance.csv", row.names = FALSE)
v1 <- read.csv("data-processed/corrected-bioshiftsv1_max-dispersal-distance.csv")

###################################################
####   calculating potential dispersal rate    ####
###################################################
### join age at maturity data with dispersal data 
am <- read.csv("data-processed/age-at-maturity.csv")

am %>% 
  ggplot(aes(x = log(AgeAtMaturity), fill = class)) + geom_histogram() +
  theme_bw()

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
## 160 / 609 species do not have age at maturity data 

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

## how variable are dispersal potential estimates within species?
sd <- v1 %>%
  group_by(Species) %>%
  summarise(sd_dp = sd(DispersalPotentialKmY), class = unique(class)) %>%
  filter(!is.na(sd_dp))

sd %>%
  ggplot(aes(x = sd_dp, fill = class)) + geom_histogram() + scale_x_log10() 
## birds have highest standard deviation in dispersal potential within species 
max(sd$sd_dp)

## calculate dispersal potential (within species)
## this time: leave multiple for each species for each dispersal estimate 
# v1 = v1 %>%
#   group_by(Species) %>%
#   mutate(MeanDispersalPotentialKmY = mean(DispersalPotentialKmY)) %>%
#   mutate(MeanDispersalPotentialmY = mean(DispersalPotentialmY)) %>%
#   mutate(MaxDispersalPotentialKmY = max(DispersalPotentialKmY)) %>%
#   mutate(MaxDispersalPotentialmY = max(DispersalPotentialmY)) 

## make dataframe that has one row per shift 
v1 <- v1 %>%
  select(#-ObservationTypeSpecific, -ObservationTypeGeneral, 
         -DispersalDistanceKm,
         -DispersalDistancem, 
         #-Code, -Field, 
         -Sex, -DispersalUnit, 
         #-Database, 
         #-DispersalSource,
         # -DispersalPotentialKmY, -DispersalPotentialmY, 
         #-DispersalDistance
         ) %>%
  distinct()

## classify species in taxonomic groups
v1 = v1 %>%
  mutate(Group = ifelse(class %in% c("Pinopsida", "Magnoliopsida", "Liliopsida"),
                        "Plants",
                        ifelse(class == "Aves",
                               "Birds",
                               ifelse(class == "Insecta",
                                      "Insects",
                                      ifelse(class == "Squamata",
                                             "Squamates",
                                             ifelse(class == "Amphibia",
                                                    "Amphibians",
                                                    ifelse(class %in% c("Actinopterygii", "Elasmobranchii"),
                                                           "Fish",
                                                           ifelse(class == "Mammalia",
                                                                  "Mammals", NA))))))))

## plot
v1 %>%
  ggplot(aes(x = Group)) + geom_bar()


############################################
####     subset data, convert units     ####
############################################

######
## CHOOSE DISPERSAL ESTIMATE HERE 
#####
saved_v1 <- v1
v1 = saved_v1
## see what type of estimate was typically the max
test = v1 %>%
  group_by(Species) %>%
  mutate(MaxDispersalPotentialKmY = max(DispersalPotentialKmY),
         MaxDispersalPotentialmY = max(DispersalPotentialmY)) %>%
  ungroup() %>%
  mutate(CodeOfMax = MaxDispersalPotentialKmY == DispersalPotentialKmY) %>%
  filter(CodeOfMax == TRUE) %>%
  filter(Group == "Birds")

test %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, fill = Field)) +
  geom_histogram() +
  scale_x_log10() +
  labs(fill = "Observation time that is the maximum within species:")

## okay, so ones labelled as maximum natal are similar to those labelled as geo mean natal:
test %>%
  filter(Field %in% c("Natal_dispersal_maximum_distance_(km)", "Geometric_mean_natal_dispersal_distance_(km)")) %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, fill = Field)) +
  geom_histogram() +
  scale_x_log10()

## and arithmetic mean is not that different than geometric mean, has a wider distribution
test %>%
  filter(Field %in% c("ArithmeticMeanNatalDispersal", "Geometric_mean_natal_dispersal_distance_(km)")) %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, fill = Field)) +
  geom_histogram() +
  scale_x_log10()

## so, if I use only natal dispersal the pattern should be the same

## recode code category 
## group xx percentile into "max"
v1$Code = ifelse(v1$Code %in% c("90thPercentileDispersalDistance",
                                "99thPercentileDispersalDistance"), 
                 "MaxDispersalDistance", v1$Code)

## split into different subsets of data
natal <- filter(v1, ObservationTypeGeneral %in% c("natal dispersal"))
breeding <- filter(v1, ObservationTypeGeneral %in% c("breeding dispersal"))
plants <- filter(v1, ObservationTypeGeneral %in% c("seed dispersal"))

## get rid of breeding 
v1 <- v1 %>%
  group_by(Species) %>%
  mutate(MaxDispersalPotentialKmY = max(DispersalPotentialKmY),
         MaxDispersalPotentialmY = max(DispersalPotentialmY)) %>%
  ungroup() %>%
  mutate(CodeOfMax = MaxDispersalPotentialKmY == DispersalPotentialKmY) %>%
  filter(CodeOfMax == TRUE) %>%
  filter(Field != "ArithmeticMeanBreedingDispersal") %>%
  select(-CodeOfMax)

## filter to geometric mean 
unique(natal$Field)
geo_natal <- filter(natal, Field %in% c("GeometricMeanNatalDispersal", 
                                        #"Natal_dispersal_maximum_distance_(km)", 
                                        "Geometric_mean_natal_dispersal_distance_(km)"))
ari_natal <- filter(natal, Field %in% c("ArithmeticMeanNatalDispersal"))

## add plant data back
ari_natal <- rbind(ari_natal, plants)
geo_natal <- rbind(geo_natal, plants)

## select max observation per species 
ari_natal <- ari_natal %>%
  group_by(Species) %>%
  mutate(MaxDispersalPotentialKmY = max(DispersalPotentialKmY),
         MaxDispersalPotentialmY = max(DispersalPotentialmY)) %>%
  ungroup() %>%
  mutate(CodeOfMax = MaxDispersalPotentialKmY == DispersalPotentialKmY) %>%
  filter(CodeOfMax == TRUE) %>%
  select(-CodeOfMax)

geo_natal <- geo_natal %>%
  group_by(Species) %>%
  mutate(MaxDispersalPotentialKmY = max(DispersalPotentialKmY),
         MaxDispersalPotentialmY = max(DispersalPotentialmY)) %>%
  ungroup() %>%
  mutate(CodeOfMax = MaxDispersalPotentialKmY == DispersalPotentialKmY) %>%
  filter(CodeOfMax == TRUE) %>%
  select(-CodeOfMax)

## look at distributions
geo_natal %>%
  rbind(., ari_natal) %>%
  filter(Group == "Birds") %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, fill = Field)) +
  geom_histogram() +
  scale_x_log10()


v1 %>%
  filter(Group == "Birds") %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, fill = Field)) +
  geom_histogram() +
  scale_x_log10()

## use max dispersal potential per species
geo_natal$dispersal_potential_kmY = geo_natal$MaxDispersalPotentialKmY
geo_natal$dispersal_potential_mY = geo_natal$MaxDispersalPotentialmY
ari_natal$dispersal_potential_kmY = ari_natal$MaxDispersalPotentialKmY
ari_natal$dispersal_potential_mY = ari_natal$MaxDispersalPotentialmY
v1$dispersal_potential_kmY = v1$MaxDispersalPotentialKmY
v1$dispersal_potential_mY = v1$MaxDispersalPotentialmY

mod_data_geo <- geo_natal %>%
  ## filter out observations where climate velocity is negative at leading edge/optimum (expect contraction)
  filter(LatVeloT >= 0 | EleVeloT >= 0) %>%
  ## get rid of negative shifts 
  filter(ShiftR > 0) %>%
  ## get rid of trailing edge 
  filter(., Position != "Trailing edge") %>%
  ## make sure none have empty dispersal potential 
  filter(!is.na(dispersal_potential_kmY)) %>%
  ## make one column for potential dispersal rate, climate velo for easier plotting of lat x elev data together
  mutate(annual_dispersal_pot = ifelse(Gradient == "Elevation",
                                       dispersal_potential_mY,
                                       ifelse(Gradient == "Latitudinal",
                                              dispersal_potential_kmY,
                                              NA))) %>%
  ## reorder factors 
  mutate(Group = factor(Group, ordered = TRUE, levels = c("Birds", "Plants", "Mammals",
                                                          "Fish", "Amphibians", "Squamates"))) 

mod_data_ari <- ari_natal %>%
  ## filter out observations where climate velocity is negative at leading edge/optimum (expect contraction)
  filter(LatVeloT >= 0 | EleVeloT >= 0) %>%
  ## get rid of negative shifts 
  filter(ShiftR > 0) %>%
  ## get rid of trailing edge 
  filter(., Position != "Trailing edge") %>%
  ## make sure none have empty dispersal potential 
  filter(!is.na(dispersal_potential_kmY)) %>%
  ## make one column for potential dispersal rate, climate velo for easier plotting of lat x elev data together
  mutate(annual_dispersal_pot = ifelse(Gradient == "Elevation",
                                       dispersal_potential_mY,
                                       ifelse(Gradient == "Latitudinal",
                                              dispersal_potential_kmY,
                                              NA))) %>%
  ## reorder factors 
  mutate(Group = factor(Group, ordered = TRUE, levels = c("Birds", "Plants", "Mammals",
                                                          "Fish", "Amphibians", "Squamates"))) 

mod_data <- v1 %>%
  ## filter out observations where climate velocity is negative at leading edge/optimum (expect contraction)
  filter(LatVeloT >= 0 | EleVeloT >= 0) %>%
  ## get rid of negative shifts 
  filter(ShiftR > 0) %>%
  ## get rid of trailing edge 
  filter(., Position != "Trailing edge") %>%
  ## make sure none have empty dispersal potential 
  filter(!is.na(dispersal_potential_kmY)) %>%
  ## make one column for potential dispersal rate, climate velo for easier plotting of lat x elev data together
  mutate(annual_dispersal_pot = ifelse(Gradient == "Elevation",
                                       dispersal_potential_mY,
                                       ifelse(Gradient == "Latitudinal",
                                              dispersal_potential_kmY,
                                              NA))) %>%
  ## reorder factors 
  mutate(Group = factor(Group, ordered = TRUE, levels = c("Birds", "Plants", "Mammals",
                                                          "Fish", "Amphibians", "Squamates"))) 


## convert units of latitudinal and elevation shifts & dispersal & climate velocity to km/y:
mod_data_ari$ShiftKmY <- ifelse(mod_data_ari$Gradient == "Elevation", mod_data_ari$ShiftR / 1000, mod_data_ari$ShiftR)
mod_data_ari$ClimVeloTKmY <- ifelse(mod_data_ari$Gradient == "Elevation", mod_data_ari$EleVeloT / 1000, mod_data_ari$LatVeloT)
mod_data_ari$AnnualDispPotKmY <- ifelse(mod_data_ari$Gradient == "Elevation", mod_data_ari$annual_dispersal_pot / 1000, mod_data_ari$annual_dispersal_pot)

mod_data_geo$ShiftKmY <- ifelse(mod_data_geo$Gradient == "Elevation", mod_data_geo$ShiftR / 1000, mod_data_geo$ShiftR)
mod_data_geo$ClimVeloTKmY <- ifelse(mod_data_geo$Gradient == "Elevation", mod_data_geo$EleVeloT / 1000, mod_data_geo$LatVeloT)
mod_data_geo$AnnualDispPotKmY <- ifelse(mod_data_geo$Gradient == "Elevation", mod_data_geo$annual_dispersal_pot / 1000, mod_data_geo$annual_dispersal_pot)

mod_data$ShiftKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$ShiftR / 1000, mod_data$ShiftR)
mod_data$ClimVeloTKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$EleVeloT / 1000, mod_data$LatVeloT)
mod_data$AnnualDispPotKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$annual_dispersal_pot / 1000, mod_data$annual_dispersal_pot)

mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)

mod_data_ari <- select(mod_data_ari, -c(SLDiff, CorrShift, PredSLShift, study_level_shift)) %>%
  distinct()
mod_data_geo <- select(mod_data_geo, -c(SLDiff, CorrShift, PredSLShift, study_level_shift)) %>%
  distinct()
mod_data <- select(mod_data, -c(SLDiff, CorrShift, PredSLShift, study_level_shift)) %>%
  distinct()

## make plots to summarize subset 
nrow(mod_data_ari) # 1565 range shifts
length(unique(mod_data_ari$Species)) # 318 species 
nrow(mod_data_geo) # 1817 range shifts
length(unique(mod_data_geo$Species)) # 379 species 
nrow(mod_data) # 1920 range shifts
length(unique(mod_data$Species)) # 404 species 

########################################################
####      fitting global breakpoint regression      ####
########################################################
## filter to shifts > 0.0001
mod_data_ari <- mod_data_ari %>%
  filter(ShiftR > 0.0001) 
mod_data_geo <- mod_data_geo %>%
  filter(ShiftR > 0.0001) 
mod_data <- mod_data %>%
  filter(ShiftR > 0.0001) 

## split between elev and lat
lat_ari = filter(mod_data_ari, Gradient == "Latitudinal")
ele_ari = filter(mod_data_ari, Gradient == "Elevation")
lat_geo = filter(mod_data_geo, Gradient == "Latitudinal")
ele_geo = filter(mod_data_geo, Gradient == "Elevation")
lat_v1 = filter(mod_data, Gradient == "Latitudinal")
ele_v1 = filter(mod_data, Gradient == "Elevation")

### ARITHMETIC
## fit normal regression 
mod_lat_ari <- lm(ShiftKmY ~ AnnualDispPotKmY, data = lat_ari)
summary(mod_lat_ari)
# Residual standard error: 2.303 on 1040 degrees of freedom
# Multiple R-squared:  0.1057,	Adjusted R-squared:  0.1049 

## fit breakpoint regression
mod_bp_lat_ari <- segmented(mod_lat_ari, 
                        seg.Z = ~ AnnualDispPotKmY) ## do not set any starting value for break point

summary(mod_bp_lat_ari)
# Residual standard error: 2.127 on 1038 degrees of freedom
# Multiple R-Squared: 0.2385,  Adjusted R-squared: 0.2363 

## get estimated breakpoint and its standard error
mod_bp_lat_ari$psi ## 4.2

## get the slopes
slope(mod_bp_lat_ari) ## 0.612, -0.002

## plot the predictions
pred_ari <- data.frame(AnnualDispPotKmY = lat_ari$AnnualDispPotKmY)
pred_ari <- data.frame(predict(mod_bp_lat_ari, newdata = pred_ari, interval = "confidence"))
ci_ari <- confint(mod_bp_lat_ari)

bp_y <- predict(mod_bp_lat_ari, newdata = data.frame(AnnualDispPotKmY = mod_bp_lat_ari$psi[,2]))

df_ari <- data.frame(AnnualDispPotKmY = lat_ari$AnnualDispPotKmY,
                 ShiftKmY = lat_ari$ShiftKmY,
                 pred_shift = pred_ari$fit,
                 min = min(lat_ari$ClimVeloTKmY), 
                 max = max(lat_ari$ClimVeloTKmY),
                 mean = mean(lat_ari$ClimVeloTKmY), 
                 breakpoint_x = mod_bp_lat_ari$psi[,2],
                 breakpoint_y = bp_y,
                 ci_lower = pred_ari$lwr,
                 ci_upper = pred_ari$upr,
                 ci_bp_lower = ci_ari[,2],
                 ci_bp_upper = ci_ari[,3])

## plot the fitted model
df_ari %>%
  ggplot(aes(x = AnnualDispPotKmY, y = pred_shift)) +
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey70", alpha = 0.7) +
  geom_point(data = lat_ari, aes(colour = ClimVeloTKmY, y = ShiftKmY), alpha = 0.5) + 
  geom_line() +
  theme_bw() +
  stat_function(colour = "black", # add 1:1 line
                fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_pointrange(aes(xmin = ci_bp_lower, xmax = ci_bp_upper, x = breakpoint_x, y = breakpoint_y),
            size = 0.3, colour = "black") +
  geom_vline(aes(xintercept = mean))  # plot theoretical point 
  # geom_rect(aes(xmin = min, xmax = max, ymin = 0, ymax = 41),
  #           alpha = 0.002)

## arithmetic mean break point BIGGER than predicted 

### GEOMETRIC
## fit normal regression 
mod_lat_geo <- lm(ShiftKmY ~ AnnualDispPotKmY, data = lat_geo)
summary(mod_lat_geo)
# Residual standard error: 2.885 on 1254 degrees of freedom
# Multiple R-squared:  0.02171,	Adjusted R-squared:  0.02093 

## fit breakpoint regression
mod_bp_lat_geo <- segmented(mod_lat_geo, 
                            seg.Z = ~ AnnualDispPotKmY) ## do not set any starting value for break point

summary(mod_bp_lat_geo)
# Residual standard error: 2.65 on 1252 degrees of freedom
# Multiple R-Squared: 0.1759,  Adjusted R-squared: 0.1739 

## get estimated breakpoint and its standard error
mod_bp_lat_geo$psi ## 0.75

## get the slopes
slope(mod_bp_lat_geo) ## 3.336, 0.004

## plot the predictions
pred_geo <- data.frame(AnnualDispPotKmY = lat_geo$AnnualDispPotKmY)
pred_geo <- data.frame(predict(mod_bp_lat_geo, newdata = pred_geo, interval = "confidence"))
ci_geo <- confint(mod_bp_lat_geo)

bp_y <- predict(mod_bp_lat_geo, newdata = data.frame(AnnualDispPotKmY = mod_bp_lat_geo$psi[,2]))

df_geo <- data.frame(AnnualDispPotKmY = lat_geo$AnnualDispPotKmY,
                     ShiftKmY = lat_geo$ShiftKmY,
                     pred_shift = pred_geo$fit,
                     min = min(lat_geo$ClimVeloTKmY), 
                     max = max(lat_geo$ClimVeloTKmY),
                     mean = mean(lat_geo$ClimVeloTKmY), 
                     breakpoint_x = mod_bp_lat_geo$psi[,2],
                     breakpoint_y = bp_y,
                     ci_lower = pred_geo$lwr,
                     ci_upper = pred_geo$upr,
                     ci_bp_lower = ci_geo[,2],
                     ci_bp_upper = ci_geo[,3])

## plot the fitted model
df_geo %>%
  ggplot(aes(x = AnnualDispPotKmY, y = pred_shift)) +
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey70", alpha = 0.7) +
  geom_point(data = lat_geo, aes(colour = ClimVeloTKmY, y = ShiftKmY), alpha = 0.5) + 
  geom_line() +
  theme_bw() +
  stat_function(colour = "black", # add 1:1 line
                fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_pointrange(aes(xmin = ci_bp_lower, xmax = ci_bp_upper, x = breakpoint_x, y = breakpoint_y),
                  size = 0.3, colour = "black") +
  geom_vline(aes(xintercept = mean))  # plot theoretical point 
# geom_rect(aes(xmin = min, xmax = max, ymin = 0, ymax = 41),
#           alpha = 0.002)

## geometric mean breakpoint smaller than predicted 

### NORMAL BUT WITHOUT BREEDING
## fit normal regression 
mod_lat_v1 <- lm(ShiftKmY ~ AnnualDispPotKmY, data = lat_v1)
summary(mod_lat_v1)
# Residual standard error: 2.885 on 1254 degrees of freedom
# Multiple R-squared:  0.02171,	Adjusted R-squared:  0.02093 
## with breeding dispersal:
# Residual standard error: 3.221 on 1364 degrees of freedom
# Multiple R-squared:  0.01355,	Adjusted R-squared:  0.01283 

## better R2 without breeding dispersal

## fit breakpoint regression
mod_bp_lat_v1 <- segmented(mod_lat_v1, 
                            seg.Z = ~ AnnualDispPotKmY) ## do not set any starting value for break point

summary(mod_bp_lat_v1)
# Residual standard error: 2.946 on 1332 degrees of freedom
# Multiple R-Squared: 0.1806,  Adjusted R-squared: 0.1787 
## with breeding dispersal:
# Residual standard error: 2.946 on 1362 degrees of freedom
# Multiple R-Squared: 0.1761,  Adjusted R-squared: 0.1743 

## similar fit without breeding dispersal

## get estimated breakpoint and its standard error
mod_bp_lat_v1$psi ## 5.03
## with breeding dispersal: 4.96

## get the slopes
slope(mod_bp_lat_v1) ## 0.59, -0.0002
## with breeding dispersal: 0.59, -0.0003

## plot the predictions
pred_v1 <- data.frame(AnnualDispPotKmY = lat_v1$AnnualDispPotKmY)
pred_v1 <- data.frame(predict(mod_bp_lat_v1, newdata = pred_v1, interval = "confidence"))
ci_v1 <- confint(mod_bp_lat_v1)

bp_y <- predict(mod_bp_lat_v1, newdata = data.frame(AnnualDispPotKmY = mod_bp_lat_v1$psi[,2]))

df_v1 <- data.frame(AnnualDispPotKmY = lat_v1$AnnualDispPotKmY,
                     ShiftKmY = lat_v1$ShiftKmY,
                     pred_shift = pred_v1$fit,
                     min = min(lat_v1$ClimVeloTKmY), 
                     max = max(lat_v1$ClimVeloTKmY),
                     mean = mean(lat_v1$ClimVeloTKmY), 
                     breakpoint_x = mod_bp_lat_v1$psi[,2],
                     breakpoint_y = bp_y,
                     ci_lower = pred_v1$lwr,
                     ci_upper = pred_v1$upr,
                     ci_bp_lower = ci_v1[,2],
                     ci_bp_upper = ci_v1[,3])

## plot the fitted model
df_v1 %>%
  ggplot(aes(x = AnnualDispPotKmY, y = pred_shift)) +
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey70", alpha = 0.7) +
  geom_point(data = lat_v1, aes(colour = ClimVeloTKmY, y = ShiftKmY), alpha = 0.5) + 
  geom_line() +
  theme_bw() +
  stat_function(colour = "black", # add 1:1 line
                fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_pointrange(aes(xmin = ci_bp_lower, xmax = ci_bp_upper, x = breakpoint_x, y = breakpoint_y),
                  size = 0.3, colour = "black") +
  geom_vline(aes(xintercept = mean))  # plot theoretical point 
# geom_rect(aes(xmin = min, xmax = max, ymin = 0, ymax = 41),
#           alpha = 0.002) 


#############################################################################
####    fitting separate breakpoint regressions per climate velocity     ####
#############################################################################
## bin by climate velocity and fit separate regression to each 
## reset data 
lat = filter(mod_data_geo, Gradient == "Latitudinal")

## remove outlier
# lat <- filter(lat, ShiftKmY < 20)

## choose 4 bins so that enough species with different dispersal abilities are sampled across different climate velocities
hist(lat$ClimVeloTKmY)
q = quantile(lat$ClimVeloTKmY, probs = c(0,0.25, 0.5, 0.75,1))
lat$ClimVeloTKmY_cont <- lat$ClimVeloTKmY # save original climate velocity as new variable
lat$ClimVeloTKmY = cut(lat$ClimVeloTKmY,
                       breaks = q,
                       include.lowest = T)
plot(lat$ClimVeloTKmY)

lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\[", "(") 
lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\]", ")") 

lat <- lat %>%
  mutate(quant_max = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,2], "\\)", " "))) %>%
  mutate(quant_min = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,1], "\\(", " "))) %>%
  mutate(quant_mean = (quant_min + quant_max)/2) 

quants <- lat %>%
  ungroup() %>%
  dplyr::select(quant_mean, quant_min, quant_max, ClimVeloTKmY) %>%
  distinct()

mycol <- rev(colorRampPalette(RColorBrewer::brewer.pal(6, "RdBu"))(4))

## plot by binned velocity 
lat %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() +
  facet_grid(~ClimVeloTKmY) +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  labs(x = "Potential dispersal rate (km/y)", y = "Observed range shift rate (km/y)",
       colour = "") +
  scale_colour_manual(values = mycol) +
  stat_function(colour = "black", linetype = "dashed", fun = function(x){x}) 

## now model:
split <- split(lat, 
               f = lat$quant_mean)

mod_lms <- lapply(split, FUN = function(x) {
  mod_lat <- lm(ShiftKmY ~ AnnualDispPotKmY, data = x)
  return(mod_lat)
}
)

mod_fits <- lapply(split, FUN = function(x) {
  mod_lat <- lm(ShiftKmY ~ AnnualDispPotKmY, data = x)
  
  fit = segmented(mod_lat,
                  seg.Z = ~ AnnualDispPotKmY)
  return(fit)
}
)

## compare fit 
summary(mod_lms[[1]])
# Residual standard error: 0.8617 on 402 degrees of freedom
# Multiple R-squared:  0.06083,	Adjusted R-squared:  0.05849 
summary(mod_fits[[1]])
# Residual standard error: 0.7935 on 400 degrees of freedom
# Multiple R-Squared: 0.2077,  Adjusted R-squared: 0.2017 

summary(mod_lms[[2]])
# Residual standard error: 4.834 on 273 degrees of freedom
# Multiple R-squared:  0.01641,	Adjusted R-squared:  0.0128 
summary(mod_fits[[2]])
# Residual standard error: 4.564 on 271 degrees of freedom
# Multiple R-Squared: 0.1295,  Adjusted R-squared: 0.1198 

summary(mod_lms[[3]])
# Residual standard error: 3.316 on 330 degrees of freedom
# Multiple R-squared:  0.002757,	Adjusted R-squared:  -0.0002647 
summary(mod_fits[[3]])
# Residual standard error: 3.258 on 328 degrees of freedom
# Multiple R-Squared: 0.04282,  Adjusted R-squared: 0.03406 

summary(mod_lms[[4]])
# Residual standard error: 2.611 on 332 degrees of freedom
# Multiple R-squared:  0.0005762,	Adjusted R-squared:  -0.002434 
summary(mod_fits[[4]])
# Residual standard error: 2.474 on 330 degrees of freedom
# Multiple R-Squared: 0.1084,  Adjusted R-squared: 0.1003 

## extract break points  
results <- lapply(mod_fits, FUN = function(x) {
  x$psi
})

results <- as.data.frame(do.call(rbind, results))
results$quant_mean <- unique(lat$quant_mean)[order(unique(lat$quant_mean))]
results <- left_join(results, quants)

## get confidence intervals 
cis <- lapply(mod_fits, FUN = function(x) {
  confint(x)
})
cis <- as.data.frame(do.call(rbind, cis))

results <- cis %>%
  select(-Est.) %>% cbind(results, .) %>%
  rename("ci_bp_lower" = `CI(95%).low`,
         "ci_bp_upper" = `CI(95%).up`)

## plot theoretical versus real break points 
results %>% 
  ggplot(aes(y = quant_mean, x = ClimVeloTKmY)) +
  theme_bw() +
  labs(y = "Breakpoint", x = "Climate velocity") +
  geom_pointrange(aes(y = Est., x = ClimVeloTKmY, 
                      ymin = ci_bp_lower,
                      ymax = ci_bp_upper),
                  colour = "blue") +
  scale_y_log10() +
  geom_point() 


## extract slopes
slopes <- lapply(mod_fits, FUN = function(x) {
  d <- do.call(rbind, slope(x))
  return(d)
})

slopes <- as.data.frame(do.call(rbind, slopes))
slopes$quant_mean <- rep(unique(lat$quant_mean)[order(unique(lat$quant_mean))], each = 2)
slopes$section <- rep(c(1,2), 4)
slopes$expected_slope <- rep(c(1,0), 4)
slopes <- left_join(slopes, quants) %>%
  rename("ci_slope_lower" = `CI(95%).l`,
         "ci_slope_upper" = `CI(95%).u`)

## plot slopes 
slopes %>% 
  ggplot(aes(y = expected_slope, x = ClimVeloTKmY)) +
  theme_bw() +
  labs(y = "Slope", x = "Climate velocity", 
       colour = "Expected slope") +
  geom_pointrange(aes(y = Est., x = ClimVeloTKmY,
                      ymin = ci_slope_lower,
                      ymax = ci_slope_upper, 
                      colour = as.factor(expected_slope))) +
  scale_colour_discrete(c("red", "blue")) +
  geom_point() 


## get the model predictions
df_all <- c()
i=1
while(i <= length(mod_fits)) {
  df = split[[i]]
  mod <- mod_fits[[i]]
  
  pred <- data.frame(AnnualDispPotKmY = df$AnnualDispPotKmY)
  pred <- data.frame(predict(mod, newdata = pred, interval = "confidence"))
  
  ci <- confint(mod)
  
  bp_y <- predict(mod, newdata = data.frame(AnnualDispPotKmY = mod$psi[,2]))
  
  df <- data.frame(AnnualDispPotKmY = df$AnnualDispPotKmY,
                   ShiftKmY = df$ShiftKmY,
                   pred_shift = pred$fit,
                   ClimVeloTKmY = df$ClimVeloTKmY,
                   quant_min = min(df$ClimVeloTKmY_cont), 
                   quant_max = max(df$ClimVeloTKmY_cont),
                   quant_mean = mean(df$ClimVeloTKmY_cont), 
                   breakpoint_x = mod$psi[,2],
                   breakpoint_y = bp_y,
                   ci_lower = pred$lwr,
                   ci_upper = pred$upr,
                   ci_bp_lower = ci[,2],
                   ci_bp_upper = ci[,3])
  
  df_all <- rbind(df_all, df)
  i=i+1
}

df <- left_join(df_all, results) %>%
  mutate(section = ifelse(AnnualDispPotKmY <= Est., 1, 2)) %>%
  rename("Est_breakpoint" = Est.) %>%
  left_join(., slopes) 

## plot the fitted model
df %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY)) + 
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey70", alpha = 0.7) +
  geom_point(data = lat, aes(colour = ClimVeloTKmY_cont), alpha = 0.5) + 
  geom_ribbon(aes(ymax = ShiftKmY + St.Err.,ymin = ShiftKmY - St.Err.),
              fill = "grey70", alpha = 0.03) + # add confidence intervals
  geom_line(aes(group = ClimVeloTKmY, y = pred_shift)) +
  theme_bw() +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +  facet_grid(~ClimVeloTKmY) +
  labs(x = "Potential dispersal rate (km/y)", y = "Observed range shift rate (km/y)",
       colour = "Mean climate\nvelocity\nacross study\narea (km/y)") +
  stat_function(colour = "black", # add 1:1 line
                fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  geom_vline(aes(xintercept = quant_mean)) + # plot theoretical point 
  geom_pointrange(aes(xmin = ci_bp_lower, xmax = ci_bp_upper, x = breakpoint_x, y = breakpoint_y),
                  size = 0.3, colour = "black") +
  # geom_rect(aes(xmin = quant_min, xmax = quant_max, ymin = 0, ymax =41), 
  #           alpha = 0.002) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 



#############################################################################
####    fitting separate breakpoint regressions per climate velocity     ####
####                        (only leading edge obs.)                     ####
#############################################################################
## reset data
lat = filter(mod_data, Gradient == "Latitudinal")

## remove outlier
# lat <- filter(lat, ShiftKmY < 20)

## remove centroid obs 
lat <- filter(lat, Position != "Centroid")

## choose 3 bins so that enough species with different dispersal abilities are sampled across different climate velocities
hist(lat$ClimVeloTKmY)
q = quantile(lat$ClimVeloTKmY, probs = c(0,0.33, 0.66,1))
lat$ClimVeloTKmY_cont <- lat$ClimVeloTKmY # save original climate velocity as new variable
lat$ClimVeloTKmY = cut(lat$ClimVeloTKmY,
                       breaks = q,
                       include.lowest = T)
plot(lat$ClimVeloTKmY)

lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\[", "(") 
lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\]", ")") 

lat <- lat %>%
  mutate(quant_max = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,2], "\\)", " "))) %>%
  mutate(quant_min = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,1], "\\(", " "))) %>%
  mutate(quant_mean = (quant_min + quant_max)/2) 

quants <- lat %>%
  ungroup() %>%
  dplyr::select(quant_mean, quant_min, quant_max, ClimVeloTKmY) %>%
  distinct()

mycol <- rev(colorRampPalette(RColorBrewer::brewer.pal(6, "RdBu"))(3))

## plot by binned velocity 
lat %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() +
  facet_grid(~ClimVeloTKmY) +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  labs(x = "Potential dispersal rate (km/y)", y = "Observed range shift rate (km/y)",
       colour = "") +
  scale_colour_manual(values = mycol) +
  stat_function(colour = "black", linetype = "dashed", fun = function(x){x}) 

## now model
split <- split(lat, 
               f = lat$quant_mean)

mod_lms <- lapply(split, FUN = function(x) {
  mod_lat <- lm(ShiftKmY ~ AnnualDispPotKmY, data = x)
  
  return(mod_lat)
}
)
mod_fits <- lapply(split, FUN = function(x) {
  mod_lat <- lm(ShiftKmY ~ AnnualDispPotKmY, data = x)
  
  fit = segmented(mod_lat,
                  seg.Z = ~ AnnualDispPotKmY)
  return(fit)
}
)

## compare fit 
summary(mod_lms[[1]])
# Residual standard error: 1.262 on 83 degrees of freedom
# Multiple R-squared:  0.07806,	Adjusted R-squared:  0.06695 
summary(mod_fits[[1]])
# Residual standard error: 0.9168 on 81 degrees of freedom
# Multiple R-Squared: 0.525,  Adjusted R-squared: 0.5074 

summary(mod_lms[[2]])
# Residual standard error: 4.637 on 81 degrees of freedom
# Multiple R-squared:  0.01382,	Adjusted R-squared:  0.001643 
summary(mod_fits[[2]])
# Residual standard error: 4.383 on 79 degrees of freedom
# Multiple R-Squared: 0.1407,  Adjusted R-squared: 0.108 

summary(mod_lms[[3]])
# Residual standard error: 2.272 on 73 degrees of freedom
# Multiple R-squared:  0.009122,	Adjusted R-squared:  -0.004452 
summary(mod_fits[[3]])
# Residual standard error: 1.839 on 71 degrees of freedom
# Multiple R-Squared: 0.3681,  Adjusted R-squared: 0.3414 


## extract break points  
results <- lapply(mod_fits, FUN = function(x) {
  x$psi
})

results <- as.data.frame(do.call(rbind, results))
results$quant_mean <- unique(lat$quant_mean)[order(unique(lat$quant_mean))]
results <- left_join(results, quants)

## get confidence intervals 
cis <- lapply(mod_fits, FUN = function(x) {
  confint(x)
})
cis <- as.data.frame(do.call(rbind, cis))

results <- cis %>%
  select(-Est.) %>% cbind(results, .) %>%
  rename("ci_bp_lower" = `CI(95%).low`,
         "ci_bp_upper" = `CI(95%).up`)

## plot theoretical versus real break points 
results %>% 
  ggplot(aes(y = quant_mean, x = ClimVeloTKmY)) +
  theme_bw() +
  labs(y = "Breakpoint", x = "Climate velocity") +
  geom_pointrange(aes(y = Est., x = ClimVeloTKmY, 
                      ymin = ci_bp_lower,
                      ymax = ci_bp_upper),
                  colour = "blue") +
  scale_y_log10() +
  geom_point() 


## extract slopes
slopes <- lapply(mod_fits, FUN = function(x) {
  d <- do.call(rbind, slope(x))
  return(d)
})

slopes <- as.data.frame(do.call(rbind, slopes))
slopes$quant_mean <- rep(unique(lat$quant_mean)[order(unique(lat$quant_mean))], each = 2)
slopes$section <- rep(c(1,2), 3)
slopes$expected_slope <- rep(c(1,0), 3)
slopes <- left_join(slopes, quants) %>%
  rename("ci_slope_lower" = `CI(95%).l`,
         "ci_slope_upper" = `CI(95%).u`)

## plot slopes 
slopes %>% 
  ggplot(aes(y = expected_slope, x = ClimVeloTKmY)) +
  theme_bw() +
  labs(y = "Slope", x = "Climate velocity", 
       colour = "Expected slope") +
  geom_pointrange(aes(y = Est., x = ClimVeloTKmY,
                      ymin = ci_slope_lower,
                      ymax = ci_slope_upper, 
                      colour = as.factor(expected_slope))) +
  scale_colour_discrete(c("red", "blue")) +
  geom_point() 


## get the model predictions
df_all <- c()
i=1
while(i <= length(mod_fits)) {
  df = split[[i]]
  mod <- mod_fits[[i]]
  
  pred <- data.frame(AnnualDispPotKmY = df$AnnualDispPotKmY)
  pred <- data.frame(predict(mod, newdata = pred, interval = "confidence"))
  
  ci <- confint(mod)
  
  bp_y <- predict(mod, newdata = data.frame(AnnualDispPotKmY = mod$psi[,2]))
  
  df <- data.frame(AnnualDispPotKmY = df$AnnualDispPotKmY,
                   ShiftKmY = df$ShiftKmY,
                   pred_shift = pred$fit,
                   ClimVeloTKmY = df$ClimVeloTKmY,
                   quant_min = min(df$ClimVeloTKmY_cont), 
                   quant_max = max(df$ClimVeloTKmY_cont),
                   quant_mean = mean(df$ClimVeloTKmY_cont), 
                   breakpoint_x = mod$psi[,2],
                   breakpoint_y = bp_y,
                   ci_lower = pred$lwr,
                   ci_upper = pred$upr,
                   ci_bp_lower = ci[,2],
                   ci_bp_upper = ci[,3])
  
  df_all <- rbind(df_all, df)
  i=i+1
}

df <- left_join(df_all, results) %>%
  mutate(section = ifelse(AnnualDispPotKmY <= Est., 1, 2)) %>%
  rename("Est_breakpoint" = Est.) %>%
  left_join(., slopes) 

## plot the fitted model
df %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY)) + 
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey70", alpha = 0.7) +
  geom_point(data = lat, aes(colour = ClimVeloTKmY_cont), alpha = 0.5) + 
  geom_ribbon(aes(ymax = ShiftKmY + St.Err.,ymin = ShiftKmY - St.Err.),
              fill = "grey70", alpha = 0.03) + # add confidence intervals
  geom_line(aes(group = ClimVeloTKmY, y = pred_shift)) +
  theme_bw() +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  facet_grid(~ClimVeloTKmY) +
  labs(x = "Potential dispersal rate (km/y)", y = "Observed range shift rate (km/y)",
       colour = "Mean climate\nvelocity\nacross study\narea (km/y)") +
  stat_function(colour = "black", # add 1:1 line
                fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  geom_vline(aes(xintercept = quant_mean)) + # plot theoretical point 
  geom_pointrange(aes(xmin = ci_bp_lower, xmax = ci_bp_upper, x = breakpoint_x, y = breakpoint_y),
                  size = 0.3, colour = "black") +
  # geom_rect(aes(xmin = quant_min, xmax = quant_max, ymin = 0, ymax =41), 
  #           alpha = 0.002) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 


###############################################################################
####      analyzing influence of outliers / points with high leverage      ####
###############################################################################
## https://online.stat.psu.edu/stat462/node/171/
## reset data
lat = filter(mod_data, Gradient == "Latitudinal")

## flag outliers and points with high leverage
sres <- rstudent(mod_bp_lat) # extract Studentized residuals (resid / stand dev)
leverage <- hatvalues(mod_bp_lat) # calculate leverage of each point

## plot
mat <- matrix(cbind(leverage, sres), ncol = 2) 
plot(x = mat[,1], y = mat[,2], ylab = "Studentized residuals", xlab = "Leverage", bty = "l", pch = 16, cex = 0.7, ylim = 
       c(-max(abs(sres)),max(abs(sres)) ))
abline(h = c(-3,0,3), lwd = 2, col = "grey", lty = 2)

## mark points with stud resid > 3
out1 <- as.numeric(as.character(rownames(lat)))[abs(sres) > 3] # outlying obs
xout1 <- leverage[abs(sres) > 3]
yout1 <- sres[abs(sres) > 3]
length(xout1) ## 34 points identified
text(y = yout1, xout1 - 0.005, cex = 0.7, labels = as.character(out1),  font = 3)

## mark points with high leverage 
cut <- 3*mean(leverage)   # heuristic rule to define 3 times the mean as  extreme
lev1 <- as.numeric(as.character(rownames(lat)))[leverage > cut] # outlying observations
ylev1 <- sres[leverage > cut]
xlev1 <- leverage[leverage > cut]
length(xlev1) ## 64 points identified 
text(y = ylev1, xlev1-0.005, cex = 0.7, labels = as.character(lev1),  font = 3, col = "red")

## add vars to data flagging outliers & high leverage points
lat$high_lev <- ifelse(leverage > cut, "Y", "N")
lat$outlier <- ifelse(abs(sres) > 3, "Y", "N")

lat_copy <- lat

## first, remove high leverage points 
lat <- filter(lat, high_lev == "N")

## fit normal regression 
mod_lat <- lm(ShiftKmY ~ AnnualDispPotKmY, data = lat)
summary(mod_lat)
# Residual standard error: 3.178 on 1279 degrees of freedom
# Multiple R-squared:  0.03289,	Adjusted R-squared:  0.03214 

## fit breakpoint regression
mod_bp_lev <- segmented(mod_lat, 
                        seg.Z = ~ AnnualDispPotKmY) ## do not set any starting value for break point

summary(mod_bp_lev)
# Residual standard error: 2.934 on 1277 degrees of freedom
# Multiple R-Squared: 0.1773,  Adjusted R-squared: 0.1753 

## get estimated breakpoint and its standard error
mod_bp_lev$psi ## 5.00

## get the slopes
slope(mod_bp_lev) ## 0.481, 0.00065

## plot the predictions
pred <- data.frame(AnnualDispPotKmY = lat$AnnualDispPotKmY)
pred <- data.frame(predict(mod_bp_lev, newdata = pred, interval = "confidence"))
ci <- confint(mod_bp_lev)

bp_y <- predict(mod_bp_lev, newdata = data.frame(AnnualDispPotKmY = mod_bp_lev$psi[,2]))

df <- data.frame(AnnualDispPotKmY = lat$AnnualDispPotKmY,
                 ShiftKmY = lat$ShiftKmY,
                 pred_shift = pred$fit,
                 min = min(lat$ClimVeloTKmY), 
                 max = max(lat$ClimVeloTKmY),
                 mean = mean(lat$ClimVeloTKmY), 
                 breakpoint_x = mod_bp_lev$psi[,2],
                 breakpoint_y = bp_y,
                 ci_lower = pred$lwr,
                 ci_upper = pred$upr,
                 ci_bp_lower = ci[,2],
                 ci_bp_upper = ci[,3])

## plot the fitted model
df %>%
  ggplot(aes(x = AnnualDispPotKmY, y = pred_shift)) +
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey70", alpha = 0.7) +
  geom_point(data = lat, aes(colour = ClimVeloTKmY, y = ShiftKmY), alpha = 0.5) + 
  geom_line() +
  theme_bw() +
  stat_function(colour = "black", # add 1:1 line
                fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_pointrange(aes(xmin = ci_bp_lower, xmax = ci_bp_upper, x = breakpoint_x, y = breakpoint_y),
                  size = 0.3, colour = "black") +
  geom_vline(aes(xintercept = mean))  # plot theoretical point 
# geom_rect(aes(xmin = min, xmax = max, ymin = 0, ymax = 41),
#           alpha = 0.002) 


## now, remove outliers
lat = lat_copy
lat <- filter(lat, outlier == "N")

## fit normal regression 
mod_lat <- lm(ShiftKmY ~ AnnualDispPotKmY, data = lat)
summary(mod_lat)
# Residual standard error: 2.085 on 1309 degrees of freedom
# Multiple R-squared:  0.02122,	Adjusted R-squared:  0.02047 

## fit breakpoint regression
mod_bp_out <- segmented(mod_lat, 
                        seg.Z = ~ AnnualDispPotKmY) ## do not set any starting value for break point

summary(mod_bp_out)
# Residual standard error: 1.816 on 1307 degrees of freedom
# Multiple R-Squared: 0.2584,  Adjusted R-squared: 0.2567 

## get estimated breakpoint and its standard error
mod_bp_out$psi ## 5.004

## get the slopes
slope(mod_bp_out) ## 0.461, -0.00013

## plot the predictions
pred <- data.frame(AnnualDispPotKmY = lat$AnnualDispPotKmY)
pred <- data.frame(predict(mod_bp_out, newdata = pred, interval = "confidence"))
ci <- confint(mod_bp_out)

bp_y <- predict(mod_bp_out, newdata = data.frame(AnnualDispPotKmY = mod_bp_out$psi[,2]))

df <- data.frame(AnnualDispPotKmY = lat$AnnualDispPotKmY,
                 ShiftKmY = lat$ShiftKmY,
                 pred_shift = pred$fit,
                 min = min(lat$ClimVeloTKmY), 
                 max = max(lat$ClimVeloTKmY),
                 mean = mean(lat$ClimVeloTKmY), 
                 breakpoint_x = mod_bp_out$psi[,2],
                 breakpoint_y = bp_y,
                 ci_lower = pred$lwr,
                 ci_upper = pred$upr,
                 ci_bp_lower = ci[,2],
                 ci_bp_upper = ci[,3])

## plot the fitted model
df %>%
  ggplot(aes(x = AnnualDispPotKmY, y = pred_shift)) +
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey70", alpha = 0.7) +
  geom_point(data = lat, aes(colour = ClimVeloTKmY, y = ShiftKmY), alpha = 0.5) + 
  geom_line() +
  theme_bw() +
  stat_function(colour = "black", # add 1:1 line
                fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_pointrange(aes(xmin = ci_bp_lower, xmax = ci_bp_upper, x = breakpoint_x, y = breakpoint_y),
                  size = 0.3, colour = "black") +
  geom_vline(aes(xintercept = mean))  # plot theoretical point 
# geom_rect(aes(xmin = min, xmax = max, ymin = 0, ymax = 41),
#           alpha = 0.002) 





## next steps for this project:
## ____________________________
## - see what people report for break point regressions, write up these results
## - use plant + bird trait models to predict potential dispersal rate and test for dispersal rate/range shift rate relationships 

## notes:
## least squares method not applicable when you want to find the no-effect range
## (full range at which there so relationship between x and y)
# the method to find the no-effect range is progressive partial regression over the range, extending the range with small steps until the regression coefficient gets significantly different from zero
# https://www.waterlog.info/pdf/Range%20of%20no%20effect.pdf
# "range of no effect (plateau)
# may found by a type of segmented regression called partial regression, that is only a regression
# over that part of the X values where no effect occurs"

## then test for independence and normality of residuals and ensure variance is evenly distributed 
## bootstrap to get confidence intervals 

## chose breaks to ensure sufficient sampling of points before and after breakpoint



## try mixed effect model 

## bin by climate velocity and fit separate regression to each 
## reset data 
lat = filter(mod_data, Gradient == "Latitudinal")

## remove outlier
# lat <- filter(lat, ShiftKmY < 20)

## choose 4 bins so that enough species with different dispersal abilities are sampled across different climate velocities
hist(lat$ClimVeloTKmY)
q = quantile(lat$ClimVeloTKmY, probs = c(0,0.25, 0.5, 0.75,1))
lat$ClimVeloTKmY_cont <- lat$ClimVeloTKmY # save original climate velocity as new variable
lat$ClimVeloTKmY = cut(lat$ClimVeloTKmY,
                       breaks = q,
                       include.lowest = T)
plot(lat$ClimVeloTKmY)

lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\[", "(") 
lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\]", ")") 

lat <- lat %>%
  mutate(quant_max = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,2], "\\)", " "))) %>%
  mutate(quant_min = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,1], "\\(", " "))) %>%
  mutate(quant_mean = (quant_min + quant_max)/2) 

quants <- lat %>%
  ungroup() %>%
  dplyr::select(quant_mean, quant_min, quant_max, ClimVeloTKmY) %>%
  distinct()

mycol <- rev(colorRampPalette(RColorBrewer::brewer.pal(6, "RdBu"))(4))

## plot by binned velocity 
lat %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() +
  facet_grid(~ClimVeloTKmY) +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  labs(x = "Potential dispersal rate (km/y)", y = "Observed range shift rate (km/y)",
       colour = "") +
  scale_colour_manual(values = mycol) +
  stat_function(colour = "black", linetype = "dashed", fun = function(x){x}) 

## now model:
mod_lat <- lme(ShiftKmY ~ AnnualDispPotKmY, 
               random = ~1|ClimVeloTKmY,
               data = lat)

summary(mod_lat)

fit = segmented(mod_lat, seg.Z = ~ AnnualDispPotKmY,
                random = list(ClimVeloTKmY = pdDiag(~1 + AnnualDispPotKmY + G0)),
                psi = mean(lat$ClimVeloTKmY_cont))
fit2 = segmented(mod_lat, seg.Z = ~ AnnualDispPotKmY,
                random = list(ClimVeloTKmY = pdDiag(~1 + AnnualDispPotKmY + U + G0)),
                psi = mean(lat$ClimVeloTKmY_cont))

summary(fit)
## low variance among slope U
## high variance among breakpoints G0
fit$psi.i
slope(fit)

AIC(mod_lat, fit, fit2)
## random bp model is best fit 

## bootstrap sample
bootsegMix(fit, B=10)
