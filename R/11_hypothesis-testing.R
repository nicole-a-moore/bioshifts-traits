## Hypotheses:
## 1. Species that are expected to be able keep up with climate change will have smaller range shift lags at the leading edge than species that are not expected to be able to keep up with climate change 
## 2. Range shift lags should increase with the difference between the velocity of climate change and the species' annual dispersal capacity 
## 3. Relationships will be clearer in taxonomic groups that actively disperse because their realized dispersal patterns rely less on external dispersal conditions (i.e., wind and ocean currents) and are less hindered by topography (more likely to actually realize their annual dispersal potential)
## 4. Dispersal potential will be less limiting across elevation than latitude (more species should be able to keep up with climate change across elevation)

library(tidyverse)
library(PNWColors)
library(gridExtra)
library(grid)
library(cowplot)
library(parallel)
library(pbapply)
library(traitdataform)
library(data.table)
library(lme4)
source("R/taxonomic-harmonization/clean_taxa_functions.R")

## read function to harmonize taxonomy 
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
v1 = read.csv("data-processed/v1-corrected-range-shifts.csv")

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

#transform study area
v1$Area <- log(v1$Area)

## subset to species with dispersal scale 
v1 <- filter(v1, scientificName %in% dscale$scientificName)
length(unique(v1$Species)) #623 species 


# ## read in bioshifts v2
# v2 = read.csv("data-raw/bioshiftsv2/Bioshifts.v2.final.csv")
# ## clean names to make them match reported names in species list
# v2$reported_name = Clean_Names(gsub("_", " ", v2$Scientific.Name), return_gen_sps = F)
# 
# spv2 <- filter(sp, v2 == 1) %>%
#   select(reported_name, scientificName)
# 
# ## yay! all are there
# ## except Coprosma 3sp
# which(!v2$reported_name %in% spv2$reported_name)
# 
# ## join to add scientific name column
# v2 = left_join(v2, spv2)
# 
# ## subset to species with dispersal scale 
# v2 <- filter(v2, scientificName %in% dscale$scientificName)
# length(unique(v2$scientificName)) #587 species 
# 
# ## make sure all the species are there:
# length(unique(c(v2$scientificName, v1$scientificName))) # 651 unique species 
# length(unique(dscale$scientificName)) # 651 unique species 

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
length(unique(v1$Species)) #still 623 species


#----------------------
## calculate one dispersal distance value per species per metric
## ex. one mean, one max
# dscale %>%
#   ggplot(aes(x = Code)) + geom_bar() + coord_flip()

v1 = v1 %>% 
  group_by(Species, Code) %>%
  mutate(Code = ifelse(Code %in% c("99thPercentileDispersalDistance",
                                   "90thPercentileDispersalDistance"),
                       "MaxDispersalDistance", 
                       Code)) %>% ## recode 99th and 90th percentile to be "max"
  mutate(DispersalDistanceKm_unique = ifelse(Code == "MaxDispersalDistance",
                                             max(DispersalDistanceKm),
                                             ifelse(Code == "MedianDispersalDistance",
                                                    median(DispersalDistanceKm), 
                                                    ifelse(Code %in% c("MeanDispersalDistance", "DispersalDistance"),
                                                           mean(DispersalDistanceKm),
                                                           NA)))) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(DispersalDistanceKm_max = max(DispersalDistanceKm)) %>%
  select(Code, Species, DispersalDistanceKm, DispersalDistanceKm_unique, DispersalDistanceKm_max, 
         everything()) %>%
  ungroup()

## look at the metrics that were larger than max 
max = v1 %>%
  filter(Code == "MaxDispersalDistance") %>%
  mutate(matches = ifelse(DispersalDistanceKm_unique == DispersalDistanceKm_max, "Y", "N")) %>%
  select(matches, Species) %>%
  unique() %>%
  filter(matches == "N")

length(unique(max$Species))
## 13 species have other dispersal distance metrics that are larger than their so-called maximum 

v1 %>%
  filter(Species %in% max$Species) %>%
  filter(DispersalDistanceKm == DispersalDistanceKm_max) %>% 
  select(Species, Code) %>% 
  unique() 
## most are mean values, some median
## mostly fish, 1 plant, 1 tree, some small mammals

## correct max dispersal distance for these species
v1 = v1 %>%
  group_by(Species, Code) %>% 
  mutate(DispersalDistanceKm_unique = ifelse(DispersalDistanceKm_unique != DispersalDistanceKm_max &
                                               Code == "MaxDispersalDistance", 
                                             DispersalDistanceKm_max,
                                             DispersalDistanceKm_unique)) %>%
  ungroup() %>%
  select(-DispersalDistanceKm_max)

v1 %>%
  select(Code, Species) %>%
  unique() %>%
  group_by(Code) %>%
  tally()

## make sure we don't drop any data
maxsp = unique(v1$Species[which(v1$Code == "MaxDispersalDistance")])
nomaxsp <- unique(v1$Species[which(!v1$Species %in% maxsp)])

## 594 species have max dispersal distance
## start with max!
v1 <- v1 %>%
  filter(Code == "MaxDispersalDistance") %>%
  select(-DispersalDistanceKm, -DispersalDistance, -DispersalUnit, -Field, 
         -DispersalSource, -ObservationTypeSpecific, -Sex) %>%
  unique()

which(nomaxsp %in% v1$Species)
which(maxsp %in% v1$Species)

v1 <- rename(v1, "MaxDispersalDistanceKm" = DispersalDistanceKm_unique)

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

## get rid of the single tunicate (probably passively disperses)
v1 <- filter(v1, scientificName != "Salix alba")
## now have 593

## get rid of duplicates
v1 <- v1 %>%
  select(-Database) %>%
  distinct()

## save dataset 
#write.csv(v1, "data-processed/bioshiftsv1_max-dispersal-distance.csv", row.names = FALSE)

###################################################
####   calculating annual dispersal potential  ####
###################################################
v1 = read.csv("data-processed/bioshiftsv1_max-dispersal-distance.csv")

### join age at maturity and longevity data with dispersal data 
am <- read.csv("data-processed/age-at-maturity.csv")
long <- read.csv("data-processed/longevity.csv")

long %>% 
  ggplot(aes(x = log(LifeSpan), fill = class)) + geom_histogram() +
  theme_bw()

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

long_join = long %>%
  rename("scientificName" = SpeciesChecked) %>%
  group_by(scientificName) %>%
  mutate(LifeSpanYears = as.numeric(as.character(LifeSpan))) %>%
  mutate(LifeSpanYears = max(LifeSpanYears)) %>% # select maximum per species 
  ungroup() %>%
  select(scientificName, LifeSpanYears) %>%
  unique() 

## join to dispersal data:
v1 <- left_join(v1, long_join, by = "scientificName") %>%
  left_join(., am_join, by = "scientificName")

## for animal species that are sessile, use age at maturity 
mob <- read.csv("data-processed/Mobility_sp_bioshiftsv1.csv")
mob$SpeciesChecked <- str_replace_all(mob$SpeciesChecked, "\\_", " ")
mob$SpeciesChecked <- Clean_Names(mob$SpeciesChecked, return_gen_sps = F)

## some manual fixes:
mob$SpeciesChecked[which(mob$SpeciesChecked == "Poecile montanus")] <- "Parus montanus"
mob$SpeciesChecked[which(mob$SpeciesChecked == "Dendrocopos major")] <- "Dendrocopus major"

## harmonize taxonomy
mob_harm <- harmonize(mob$SpeciesChecked)

notfound <- filter(mob_harm, is.na(db_code))

## rename columns 
mob <- mob %>%
  rename("reported_name" = SpeciesChecked) %>%
  mutate(reported_name_fixed = reported_name)

mob <- left_join(mob, mob_harm, by = c("reported_name_fixed" = "species")) %>%
  unique()

## make sure all the animal are there 
animals <- filter(v1, kingdom == "Animalia")
View(animals[which(!animals$scientificName %in% mob$scientificName),])

## 3 are not 

## left join mobility data with v1 
v1_test <- mob %>% 
  select(scientificName, Mobility, Locomotion_mode) %>%
  unique() %>%
  left_join(v1, .)

## see how many animals are sessile 
v1_test %>%
  filter(kingdom == "Animalia") %>%
  select(scientificName, Mobility, Locomotion_mode) %>%
  unique() %>%
  group_by(Mobility, Locomotion_mode) %>%
  tally()
## none - all have active dispersal locomotion types 

## check ones that are NA:
v1_test %>%
  filter(kingdom == "Animalia") %>%
  unique() %>%
  filter(is.na(Mobility)) 
## fish and rodent, active dispersers so it's ok

## get rid of duplicated data (some are same species but have different scientific names)
v1 <- select(v1, -scientificName) %>%
  distinct()

length(unique(v1$Species)) ## still have all the species
length(unique(v1$Species[which(is.na(v1$AgeAtMaturityDays) & is.na(v1$LifeSpanYears))])) 
## 142 / 593 species do not have longevity/age at maturity data 

unique(v1$Species[which(is.na(v1$AgeAtMaturityDays) & is.na(v1$LifeSpanYears))])

## calculate dispersal potential for species with age at maturity/longevity 
v1 = v1 %>%
  mutate(MaxDispersalPotentialKmY = ifelse(!is.na(LifeSpanYears), 
                                           MaxDispersalDistanceKm/LifeSpanYears,
                                           ifelse(!is.na(YearOfMaturity),
                                                  MaxDispersalDistanceKm/YearOfMaturity,
                                                  NA))) %>%
  mutate(MaxDispersalPotentialmY = ifelse(!is.na(LifeSpanYears), 
                                          (MaxDispersalDistanceKm*1000)/LifeSpanYears,
                                          ifelse(!is.na(YearOfMaturity),
                                                 (MaxDispersalDistanceKm*1000)/YearOfMaturity,
                                                 NA))) %>%
  mutate(MaxDispersalDistancem = MaxDispersalDistanceKm*1000)

ggplot(v1, aes(x = log(MaxDispersalPotentialKmY), fill = class)) + geom_histogram()
ggplot(v1, aes(x = log(MaxDispersalPotentialmY), fill = class)) + geom_histogram()
ggplot(v1, aes(x = log(MaxDispersalDistanceKm), fill = class)) + geom_histogram()



############################################
####   calculate lags and expectations  ####
############################################
lags <- v1 %>%
  filter(!is.na(corrected_shift)) %>%
  ## calculate lag
  mutate(lag = ifelse(Gradient == "Elevation",
                      corrected_shift - EleVeloT,
                      ifelse(Gradient == "Latitudinal",
                             corrected_shift - LatVeloT, 
                             NA))) %>%
  mutate(raw_lag = ifelse(Gradient == "Elevation",
                      ShiftR - EleVeloT,
                      ifelse(Gradient == "Latitudinal",
                             ShiftR - LatVeloT, 
                             NA))) %>%
  ## filter out observations where climate velocity is negative at leading edge/optimum (expect contraction)
  filter(LatVeloT >= 0 | EleVeloT >= 0) %>%
  ## get rid of trailing edge 
  filter(., Position != "Trailing edge") %>%
  ## make sure none have empty dispersal potential 
  filter(!is.na(MaxDispersalPotentialKmY)) %>%
  ## create variable to describe whether or not species can keep up with climate change
  mutate(expect_tracking = ifelse(Gradient == "Elevation" & MaxDispersalPotentialmY >= EleVeloT,
                                  "Yes", 
                                  ifelse(Gradient == "Elevation" & MaxDispersalPotentialmY < EleVeloT,
                                         "No",
                                         ifelse(Gradient == "Latitudinal" & MaxDispersalPotentialKmY >= LatVeloT,
                                                "Yes",
                                                ifelse(Gradient == "Latitudinal" & MaxDispersalPotentialKmY < LatVeloT,
                                                       "No",
                                                       NA))))) %>%
  ## reorder factors 
  mutate(expect_tracking = factor(expect_tracking, ordered = TRUE, levels = c("Yes", "No"))) %>%
  mutate(Group = factor(Group, ordered = TRUE, levels = c("Birds", "Plants", "Mammals",
                                                          "Fish", "Amphibians", "Squamates"))) %>%
  ## calculate expected lag rate 
  mutate(expected_lag = ifelse(Gradient == "Elevation",
                                  MaxDispersalPotentialmY - EleVeloT,
                                  ifelse(Gradient == "Latitudinal",
                                         MaxDispersalPotentialKmY - LatVeloT,
                                         NA))) %>%
  ## make one column for annual dispersal potential, climate velo for easier plotting of lat x elev data together
  mutate(annual_dispersal_pot = ifelse(Gradient == "Elevation",
                                       MaxDispersalPotentialmY,
                                       ifelse(Gradient == "Latitudinal",
                                              MaxDispersalPotentialKmY,
                                              NA))) %>%
  mutate(climate_velocity = ifelse(Gradient == "Elevation",
                                       EleVeloT,
                                       ifelse(Gradient == "Latitudinal",
                                              LatVeloT,
                                              NA)))
 
## view how many observations we have 
lags %>%
  ggplot(aes(x = Position)) + geom_bar() + theme_bw() +
  labs(y = "Count", x = "Range shift parameter") + 
  facet_wrap(~Gradient)

lags %>%
  select(MaxDispersalPotentialKmY, Group, Species) %>%
  unique() %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, fill = Group)) + geom_histogram() + theme_bw() +
  scale_x_log10() + 
  labs(y = "Count", x = "Annual dispersal potential (km/y)") 


#################################
####   visualizing the data  ####
#################################
pal = pnw_palette("Bay",7)

## write out dataset for meeting with Brunno and Lise 
write.csv(lags, row.names = F, "data-processed/dispersal-data_modelling-meeting.csv")

## 1. How many species should vs. shouldn't be able to keep up with temperature change?
## 2. Is dispersal potential less limiting across elevation than latitude?
hist_gradient <- lags %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, fill = expect_tracking)) + 
  geom_histogram() +
  theme_classic() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  facet_grid(~Gradient) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), 
                labels = c("0.0001","0.001", "0.01", "0.1", "1", "10", "100", "1000")) +
  labs(x = "Annual dispersal potential (km/y)",
       y = "Number of range shift observations", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?")

ggsave(hist_gradient, height = 3.5, width = 8, device = "png", path = "figures/dispersal", 
       filename = "histogram_dispersal-potential-across-gradients.png")

hist_taxa <- lags %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, fill = expect_tracking)) + 
  geom_histogram() +
  theme_bw() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  facet_grid(Gradient~Group) +
  labs(x = "Annual dispersal potential (km/y)",
       y = "Number of range shift observations", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), 
                labels = c("1e-4    ","1e-3", " 0.01", "0.1", "1", "10", "100", "1000"))  +
  theme(panel.grid = element_blank(), strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none")

ggsave(hist_taxa, height = 3.5, width = 11.5, device = "png", path = "figures/dispersal", 
       filename = "histogram_dispersal-potential-across-groups.png")

bars = lags %>%
  ggplot(aes(x = Group, fill = expect_tracking)) + 
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  coord_flip() +
  facet_wrap(~Gradient) +
  labs(x = "Taxonomic group", 
       y = "Number of range shift observations", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?")

ggsave(bars, height = 3.5, width = 8.5, device = "png", path = "figures/dispersal", 
       filename = "barplot_dispersal-potential-across-groups.png")


## 2. Do species that are expected to be able keep up with climate change have less range shift lags than species that are not expected to be able to keep up with climate change?
ele_box_all <- lags %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (m/y)", 
       y = "Inferred range shift lag (m/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Position~Gradient) +
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(-15, 9))+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

lat_box_all <- lags %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Position~Gradient) +
  scale_y_continuous(limits = c(-15, 9)) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

box_all = grid.arrange(ele_box_all, lat_box_all, nrow = 1, widths = c(0.38, 0.62))

ggsave(box_all, height = 3.5, width = 7, device = "png", path = "figures/dispersal", 
       filename = "dp-vs-inferred-shift_boxplot.png")

ele_box_all <- lags %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (m/y)", 
       y = "Inferred range shift lag (m/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Position~Gradient) +
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(-15, 9))+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

lat_box_all <- lags %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Position~Gradient) +
  scale_y_continuous(limits = c(-15, 9)) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

box_all = grid.arrange(ele_box_all, lat_box_all, nrow = 1, widths = c(0.38, 0.62))

ggsave(box_all, height = 3.5, width = 7, device = "png", path = "figures/dispersal", 
       filename = "dp-vs-inferred-shift_points.png")

## now plot by group
## LEADING EDGE
points_lat = lags %>%
  filter(Position == "Leading edge") %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

points_ele <- lags %>%
  filter(Position == "Leading edge") %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (m/y)", 
       y = "Inferred range shift lag (m/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       title = "Leading edge") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none")

pts_all = grid.arrange(points_ele, points_lat, nrow = 2)

ggsave(pts_all, height = 5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "dp-vs-inferred-shift_points_by-group_le.png")

box_lat = lags %>%
  filter(Position == "Leading edge") %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

box_ele <- lags %>%
  filter(Position == "Leading edge") %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (m/y)", 
       y = "Inferred range shift lag (m/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       title = "Leading edge") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none")

pts_all = grid.arrange(box_ele, box_lat, nrow = 2)

ggsave(pts_all, height = 5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "dp-vs-inferred-shift_boxes_by-group_le.png")

## raw lags
box_lat = lags %>%
  filter(Position == "Leading edge") %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = raw_lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Raw range shift lag (km/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

box_ele <- lags %>%
  filter(Position == "Leading edge") %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = annual_dispersal_pot, y = raw_lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (m/y)", 
       y = "Raw range shift lag (m/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       title = "Leading edge") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none")

pts_all = grid.arrange(box_ele, box_lat, nrow = 2)

ggsave(pts_all, height = 5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "dp-vs-inferred-shift_boxes_by-group_le_raw.png")

## CENTROID
points_lat = lags %>%
  filter(Position == "Centroid") %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

points_ele <- lags %>%
  filter(Position == "Centroid") %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (m/y)", 
       y = "Inferred range shift lag (m/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       title = "Centroid") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none") 

pts_all = grid.arrange(points_ele, points_lat, nrow = 2)

ggsave(pts_all, height = 5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "dp-vs-inferred-shift_points_by-group_cent.png")

box_lat = lags %>%
  filter(Position == "Centroid") %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

box_ele <- lags %>%
  filter(Position == "Centroid") %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (m/y)", 
       y = "Inferred range shift lag (m/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       title = "Centroid") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none")

pts_all = grid.arrange(box_ele, box_lat, nrow = 2)

ggsave(pts_all, height = 5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "dp-vs-inferred-shift_boxes_by-group_cent.png")

## raw lag
box_lat = lags %>%
  filter(Position == "Centroid") %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = raw_lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Raw range shift lag (km/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

box_ele <- lags %>%
  filter(Position == "Centroid") %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = annual_dispersal_pot, y = raw_lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (m/y)", 
       y = "Raw range shift lag (m/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       title = "Centroid") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none")

pts_all = grid.arrange(box_ele, box_lat, nrow = 2)

ggsave(pts_all, height = 5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "dp-vs-inferred-shift_boxes_by-group_cent_raw.png")

## 3. Does observed range shift lag should increase with expected shift lag? (i.e., the difference between the velocity of climate change and the species' annual dispersal capacity)

lat_lags = lags %>%
  ## get rid of species who are expected to track climate change
  filter(expected_lag <= 0) %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = expected_lag, y = lag, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  scale_colour_manual(values = c(pal[6])) + 
  labs(x = "Expected annual range shift lag (km/y)", 
       y = "Inferred annual range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  facet_grid(Position~Gradient) +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = 0, colour = "darkgrey") +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme(panel.grid = element_blank(), 
        strip.background = element_rect(colour="black", fill="white")) +
  labs(y = "") +
  scale_x_continuous(limits = c(-14, 1)) +
  scale_y_continuous(limits = c(-14, 4)) 

ele_lags = lags %>%
  ## get rid of species who are expected to track climate change
  filter(expected_lag <= 0) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = expected_lag, y = lag, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))+
  scale_colour_manual(values = c(pal[6])) + 
  labs(x = "Expected annual range shift lag (m/y)", 
       y = "Inferred annual range shift lag (m/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  facet_grid(Position~Gradient) +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = 0, colour = "darkgrey") +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  scale_x_continuous(limits = c(-14, 1)) +
  scale_y_continuous(limits = c(-14, 4)) 

no = grid.arrange(ele_lags, lat_lags, nrow = 1, widths = c(0.4, 0.6))

ggsave(no, height = 5.5, width = 9, device = "png", path = "figures/dispersal", 
       filename = "expected-vs-inferred-shift_no.png")

## is it because there is structure to the data?
lat_lags_group = lags %>%
  ## get rid of species who are expected to track climate change
  filter(expected_lag <= 0) %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = expected_lag, y = lag, colour = Group)) + 
  geom_point() +
  theme_bw() +
  labs(x = "Expected annual range shift lag (km/y)", 
       y = "Inferred annual range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  facet_grid(Position~Gradient) +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = 0, colour = "darkgrey") +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  labs(y = "") +
  scale_x_continuous(limits = c(-14, 1))+
  scale_y_continuous(limits = c(-14, 4)) +
  scale_colour_manual(values = c(pal))

ele_lags_group = lags %>%
  ## get rid of species who are expected to track climate change
  filter(expected_lag <= 0) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = expected_lag, y = lag, colour = Group)) + 
  geom_point() +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  labs(x = "Expected annual range shift lag (m/y)", 
       y = "Inferred annual range shift lag (m/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  facet_grid(Position~Gradient) +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = 0, colour = "darkgrey") +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  scale_x_continuous(limits = c(-14, 1)) +
  scale_y_continuous(limits = c(-14, 4)) + 
  scale_colour_manual(values = c(pal[2]))

grid.arrange(ele_lags_group, lat_lags_group, nrow = 1, widths = c(0.4, 0.6))


## expect no relationship between excess annual dispersal potential and range shift lag
lat_excess = lags %>%
  ## get rid of species who aren't expected to track climate change
  filter(expected_lag > 0) %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = expected_lag, y = lag, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  scale_colour_manual(values = c(pal[4])) + 
  labs(x = "Excess annual dispersal potential (km/y)", 
       y = "Inferred annual range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  facet_grid(Position~Gradient) +
  geom_vline(xintercept = 0, colour = "darkgrey") +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  labs(y = "") +
  scale_y_continuous(limits = c(-15, 10))

ele_excess = lags %>%
  ## get rid of species who aren't expected to track climate change
  filter(expected_lag > 0) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = expected_lag, y = lag, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  scale_colour_manual(values = c(pal[4])) + 
  labs(x = "Excess annual dispersal potential (m/y)", 
       y = "Inferred annual range shift lag (m/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  facet_grid(Position~Gradient) +
  geom_vline(xintercept = 0, colour = "darkgrey") +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  scale_y_continuous(limits = c(-15, 10))

yes = grid.arrange(ele_excess, lat_excess, nrow = 1, widths = c(0.4, 0.6))

ggsave(yes, height = 5.5, width = 9, device = "png", path = "figures/dispersal", 
       filename = "expected-vs-inferred-shift_yes.png")

## how similar are the distributions of inferred lags and expected dispersal lags?
lags %>%
  ## get rid of species who are expected to track climate change
  filter(expected_lag <= 0) %>%
  gather(key = "lag_type", value = "lag", c(expected_lag, lag)) %>%
  ggplot(aes(x = lag, fill = lag_type)) + 
  geom_histogram() +
  theme_bw() +
  scale_fill_manual(values = pal) + 
  labs(x = "Expected annual range shift lag (km/y)", 
       y = "Inferred annual range shift lag (km/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  facet_grid(lag_type~Position~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

  
######################
####   modelling  ####
######################
## make more automatic 
## get rid of groups that have only species who can or only species who can't keep up with climate change 
lat <- filter(lags, 
              Group %in% c("Birds", "Plants", "Fish") & Gradient == "Latitudinal") %>%
  mutate(Group = factor(.$Group, ordered = F))
ele <- filter(lags, 
              Group %in% c("Plants") & Gradient == "Elevation") 

## now model
lat_mod <- lm(lag ~ expect_tracking*Group*Position, data = lat)
summary(lat_mod)
plot(lat_mod)


ncvTest(lat_mod)
## residual error does not have constant variance
## zero bound data - Poisson?
## ask model to allow variance to get smaller depending on value of y 


ele_mod <- lm(lag ~ expect_tracking*Position, data = ele)
summary(ele_mod)
plot(ele_mod)

## now plot predictions 
lat_pred <- expand.grid(Group = unique(lat$Group), expect_tracking = unique(lat$expect_tracking),
                        Position = unique(ele$Position)) 
pred <- predict(lat_mod, lat_pred, se.fit = TRUE)
lat_pred$predicted_val = pred$fit
lat_pred$predicted_se = pred$se.fit
lat_pred$expect_tracking = factor(lat_pred$expect_tracking, levels = c("No", "Yes"),
                                  ordered = T)

ele_pred <- expand.grid(Group = unique(ele$Group), expect_tracking = unique(ele$expect_tracking),
                        Position = unique(ele$Position)) 
pred <- predict(ele_mod, ele_pred, se.fit = TRUE)
ele_pred$predicted_val = pred$fit
ele_pred$predicted_se = pred$se.fit
ele_pred$expect_tracking = factor(ele_pred$expect_tracking, levels = c("No", "Yes"),
                                  ordered = T)

predplot_l <- ggplot(data = lat_pred, aes(x = expect_tracking, colour = expect_tracking, y = predicted_val,
                                          shape = Position)) + 
  geom_point() +
  geom_errorbar(aes(ymax = predicted_val + predicted_se,
                ymin = predicted_val - predicted_se)) +
  facet_wrap(~Group) +
  theme_classic() +
  labs(x = "Expect dispersal potential to allow species to keep up with temp change?",
       y = "Model-predicted inferred range shift lag (km/y)",
       title = "Latitude") +
  scale_colour_manual(values = c(pal[6], pal[4])) +
  theme(legend.position = "none")

predplot_e <- ggplot(data = ele_pred, aes(x = expect_tracking, colour = expect_tracking, y = predicted_val,
                                          shape = Position)) + 
  geom_point() +
  geom_errorbar(aes(ymax = predicted_val + predicted_se,
                    ymin = predicted_val - predicted_se)) +
  facet_wrap(~Group) +
  theme_classic() +
  labs(x = "", title = "Elevation",
       y = "Model-predicted inferred range shift lag (m/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) +
  guides(colour = "none")

gr <- grid.arrange(predplot_l, predplot_e, nrow = 1, widths = c(0.6, 0.4))

ggsave(gr, height = 4, width = 9, device = "png", path = "figures/dispersal", 
       filename = "pred_inferred.png")




#######################
####    garbage   #####
#######################




lags %>%
  ggplot(aes(x = annual_dispersal_pot, y = corrected_shift, colour = climate_velocity)) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  facet_grid(~Gradient)

lags %>%
  filter(expect_tracking == "No") %>%
  ggplot(aes(x = annual_dispersal_pot, y = corrected_shift, colour = Group)) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  facet_grid(Group~Gradient)


## look at trends within studies
lags %>%
  filter(Position == "Leading edge") %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Birds") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_boxplot(aes(fill = expect_tracking)) +
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Reference) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

lags %>%
  filter(Position == "Leading edge") %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Plants") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  geom_boxplot(aes(fill = expect_tracking)) +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Reference) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 


lags %>%
  filter(Position == "Leading edge") %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Fish") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  geom_boxplot(aes(fill = expect_tracking)) +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Reference) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 


lags %>%
  filter(Position == "Centroid") %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Birds") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_boxplot(aes(fill = expect_tracking)) +
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Reference) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

lags %>%
  filter(Position == "Centroid") %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Plants") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  geom_boxplot(aes(fill = expect_tracking)) +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Reference) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 


lags %>%
  filter(Position == "Centroid") %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Fish") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  geom_boxplot(aes(fill = expect_tracking)) +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Reference) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

lags %>%
  filter(Position == "Leading edge") %>%
  filter(Gradient == "Elevation") %>%
  filter(Group == "Plants") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  geom_boxplot(aes(fill = expect_tracking)) +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Reference) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

lags %>%
  filter(Position == "Centroid") %>%
  filter(Gradient == "Elevation") %>%
  filter(Group == "Plants") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  geom_boxplot(aes(fill = expect_tracking)) +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Reference) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 


## potential problem: velocity of temperature is the same within study areas
## within studies, pattern doesn't hold but across studies it does 
## I think only a problem if there is sampling bias (if species with different dispersal abilities are spread unevenly across sites with different climate velocity)
## it's bad if this relationship is positive (ie. species with farther dispersal are sampled in places where climate is moving faster)
## it's also bad if this relationship is negative (ie. species with shorter dispersal are sampled in places where climate is moving faster - this would cause us to see larger lags in species )

lags %>%
  ggplot(aes(x = climate_velocity, y = annual_dispersal_pot, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  facet_wrap(~Gradient) +
  scale_y_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) 

lags %>%
  ggplot(aes(x = climate_velocity, y = annual_dispersal_pot, colour = Reference)) + 
  geom_point() +
  theme_bw() +
  facet_grid(Position~Gradient) +
  scale_y_log10() +
  geom_smooth(method = "lm") 

## is climate velocity lower in sp. we expect to be able to track?
lags %>%
  ggplot(aes(x = expect_tracking, y = climate_velocity)) + 
  geom_violin() +
  theme_bw() +
  facet_wrap(Position~Gradient) + 
  geom_jitter(height = 0, width = 0.1)


## try binning by climate velocity and seeing whether relationship holds within places with similar climate velocity?
lags %>%
  mutate(cv_binned = ifelse(climate_velocity > 5, "high", "low")) %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  facet_wrap(Position~Gradient~cv_binned) + 
  scale_colour_manual(values = c(pal[4], pal[6])) +
  scale_x_log10() 
 

test = filter(lags, Gradient == "Elevation", Position == "Leading edge")

q1=quantile(test$climate_velocity,probs=c(0,0.25,0.5, 0.75,1))
cut(test$climate_velocity,breaks=q1,include.lowest=T)

min = min(test$climate_velocity)
max = max(test$climate_velocity)
cut(test$climate_velocity,breaks=c(min, (min+max)/2, max),include.lowest=T)


test %>%
  mutate(expect_tracking = factor(expect_tracking, levels = c("No", "Yes"), ordered = TRUE)) %>%
  mutate(cv_binned = cut(test$climate_velocity, breaks=c(min, (min+max)/3, 2*(min+max)/3, max), include.lowest=T)) %>%
  ggplot(aes(x = cv_binned, y = lag, colour = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  facet_wrap(Position~Gradient~Group) + 
  scale_colour_manual(values = c(pal[6], pal[4])) 

  
  
  








lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Plants") %>%
  ggplot(., aes(x = lag, fill = Source)) + geom_histogram() 

lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Plants") %>%
  ggplot(., aes(x = corrected_shift, fill = Source)) + geom_histogram() 

lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Birds") %>%
  ggplot(., aes(x = lag, fill = Source)) + geom_histogram() 

lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Birds") %>%
  ggplot(., aes(x = corrected_shift, fill = Source)) + geom_histogram() 

lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Fish") %>%
  ggplot(., aes(x = lag, fill = Source)) + geom_histogram() 

lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Fish") %>%
  ggplot(., aes(x = corrected_shift, fill = Source)) + geom_histogram() 

lags %>%
  filter(Gradient == "Elevation") %>%
  filter(Group == "Plants") %>%
  ggplot(aes(x = lag, fill = Source)) + geom_histogram()

lags %>%
  filter(Gradient == "Elevation") %>%
  filter(Group == "Plants") %>%
  ggplot(aes(x = corrected_shift, fill = Source)) + geom_histogram()
