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
library(MuMIn)
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

#transform study area
v1$Area <- log(v1$Area)

## subset to species with dispersal scale 
v1 <- filter(v1, scientificName %in% dscale$scientificName)
length(unique(v1$Species)) #696 species 


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
#write.csv(v1, "data-processed/corrected-bioshiftsv1_max-dispersal-distance.csv", row.names = FALSE)

###################################################
####   calculating annual dispersal potential  ####
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

## for now, calculate the mean and max dispersal potential within species 
v1 = v1 %>%
  group_by(Species) %>%
  mutate(MeanDispersalPotentialKmY = mean(DispersalPotentialKmY)) %>%
  mutate(MeanDispersalPotentialmY = mean(DispersalPotentialmY)) %>%
  mutate(MaxDispersalPotentialKmY = max(DispersalPotentialKmY)) %>%
  mutate(MaxDispersalPotentialmY = max(DispersalPotentialmY)) 

v1 %>%
  gather(key = "type", value = "disp_pot", c(MaxDispersalPotentialKmY, MeanDispersalPotentialKmY)) %>%
  ggplot(aes(x = disp_pot, fill = type)) + geom_histogram() + scale_x_log10() +
  facet_wrap(~class)

## difference between max and mean for birds is sometimes large 
v1 %>%
  mutate(diff = MaxDispersalPotentialKmY - MeanDispersalPotentialKmY) %>%
  ggplot(aes(x = diff, fill = class)) + geom_histogram() + scale_x_log10() 

## make dataframe that has one row per shift 
v1 <- v1 %>%
  select(-ObservationTypeSpecific, -ObservationTypeGeneral, -DispersalDistanceKm,
         -DispersalDistancem, -Code, -Field, -Sex, -DispersalUnit, -Database, -DispersalSource,
         -DispersalPotentialKmY, -DispersalPotentialmY, -DispersalDistance) %>%
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
####   calculate lags and expectations  ####
############################################
## CHOOSE WHETHER MEAN OR MAX HERE
v1$dispersal_potential_kmY = v1$MeanDispersalPotentialKmY
v1$dispersal_potential_mY = v1$MeanDispersalPotentialmY

lags <- v1 %>%
  filter(!is.na(CorrShift)) %>%
  ## calculate lag
  mutate(lag = ifelse(Gradient == "Elevation",
                      CorrShift - EleVeloT,
                      ifelse(Gradient == "Latitudinal",
                             CorrShift - LatVeloT, 
                             NA))) %>%
  mutate(raw_lag = ifelse(Gradient == "Elevation",
                      ShiftR - EleVeloT,
                      ifelse(Gradient == "Latitudinal",
                             ShiftR - LatVeloT, 
                             NA))) %>%
  ## calculate percentage tracking 
  mutate(percent_lag = ifelse(Gradient == "Elevation",
                          (ShiftR/EleVeloT),
                          ifelse(Gradient == "Latitudinal",
                                 (ShiftR/LatVeloT), 
                                 NA))) %>%
  ## filter out observations where climate velocity is negative at leading edge/optimum (expect contraction)
  filter(LatVeloT >= 0 | EleVeloT >= 0) %>%
  ## get rid of negative shifts 
  filter(ShiftR > 0) %>%
  ## get rid of trailing edge 
  filter(., Position != "Trailing edge") %>%
  ## make sure none have empty dispersal potential 
  filter(!is.na(dispersal_potential_kmY)) %>%
  ## create variable to describe whether or not species can keep up with climate change
  mutate(expect_tracking = ifelse(Gradient == "Elevation" & dispersal_potential_mY >= EleVeloT,
                                  "Yes", 
                                  ifelse(Gradient == "Elevation" & dispersal_potential_mY < EleVeloT,
                                         "No",
                                         ifelse(Gradient == "Latitudinal" & dispersal_potential_kmY >= LatVeloT,
                                                "Yes",
                                                ifelse(Gradient == "Latitudinal" & dispersal_potential_kmY < LatVeloT,
                                                       "No",
                                                       NA))))) %>%
  ## reorder factors 
  mutate(expect_tracking = factor(expect_tracking, ordered = TRUE, levels = c("Yes", "No"))) %>%
  mutate(Group = factor(Group, ordered = TRUE, levels = c("Birds", "Plants", "Mammals",
                                                          "Fish", "Amphibians", "Squamates"))) %>%
  ## make one column for annual dispersal potential, climate velo for easier plotting of lat x elev data together
  mutate(annual_dispersal_pot = ifelse(Gradient == "Elevation",
                                       dispersal_potential_mY,
                                       ifelse(Gradient == "Latitudinal",
                                              dispersal_potential_kmY,
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
  select(dispersal_potential_kmY, Group, Species) %>%
  unique() %>%
  ggplot(aes(x = dispersal_potential_kmY, fill = Group)) + geom_histogram() + theme_bw() +
  scale_x_log10() + 
  labs(y = "Count", x = "Annual dispersal potential (km/y)") 


#################################
####   visualizing the data  ####
#################################
pal = pnw_palette("Bay",7)

## 1. How many species should vs. shouldn't be able to keep up with temperature change?
## 2. Is dispersal potential less limiting across elevation than latitude?
hist_gradient <- lags %>%
  ggplot(aes(x = dispersal_potential_kmY, fill = expect_tracking)) + 
  geom_histogram() +
  theme_classic() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  facet_grid(Group~Gradient) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), 
                labels = c("0.0001","0.001", "0.01", "0.1", "1", "10", "100", "1000")) +
  labs(x = "Annual dispersal potential (km/y)",
       y = "Number of range shift observations", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?")

ggsave(hist_gradient, height = 3.5, width = 8, device = "png", path = "figures/corrected-shifts/", 
       filename = "histogram_dispersal-potential-across-gradients-and-taxa.png")

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

ggsave(bars, height = 3.5, width = 8.5, device = "png", path = "figures/corrected-shifts/", 
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
       y = "Corrected range shift lag (m/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Position~Gradient) +
  scale_y_continuous(limits = c(-14, 45)) +
  theme(legend.position = "none")+
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
       y = "Corrected range shift lag (km/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Position~Gradient) +
  scale_y_continuous(limits = c(-14, 45)) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

box_all = grid.arrange(ele_box_all, lat_box_all, nrow = 1, widths = c(0.38, 0.62))

ggsave(box_all, height = 3.5, width = 7, device = "png", path = "figures/corrected-shifts/", 
       filename = "dp-vs-raw-shift_boxplot.png")

ele_box_all <- lags %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (m/y)", 
       y = "Corrected range shift lag (m/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Position~Gradient) +
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(-14, 45))+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

lat_box_all <- lags %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = CorrShift, colour = expect_tracking)) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  scale_colour_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Corrected range shift lag (km/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Position~Gradient) +
  scale_y_continuous(limits = c(-14, 45)) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"))

box_all = grid.arrange(ele_box_all, lat_box_all, nrow = 1, widths = c(0.38, 0.62))

ggsave(box_all, height = 3.5, width = 7, device = "png", path = "figures/corrected-shifts/", 
       filename = "dp-vs-raw-shift_points.png")

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
       y = "Corrected range shift lag (km/y)", 
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
       y = "Corrected range shift lag (m/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       title = "Leading edge") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none")

pts_all = grid.arrange(points_ele, points_lat, nrow = 2)

ggsave(pts_all, height = 5, width = 11, device = "png", path = "figures/corrected-shifts", 
       filename = "dp-vs-raw-shift_points_by-group_le.png")

box_lat = lags %>%
  filter(Position == "Leading edge") %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Corrected range shift lag (km/y)", 
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
       y = "Corrected range shift lag (m/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       title = "Leading edge") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none")

pts_all = grid.arrange(box_ele, box_lat, nrow = 2)

ggsave(pts_all, height = 5, width = 11, device = "png", path = "figures/corrected-shifts", 
       filename = "dp-vs-raw-shift_boxes_by-group_le.png")


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
       y = "Corrected range shift lag (km/y)", 
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
       y = "Corrected range shift lag (m/y)", 
       colour = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       title = "Centroid") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none") 

pts_all = grid.arrange(points_ele, points_lat, nrow = 2)

ggsave(pts_all, height = 5, width = 11, device = "png", path = "figures/corrected-shifts", 
       filename = "dp-vs-raw-shift_points_by-group_cent.png")

box_lat = lags %>%
  filter(Position == "Centroid") %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Corrected range shift lag (km/y)", 
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
       y = "Corrected range shift lag (m/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?",
       title = "Centroid") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Gradient~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none")

pts_all = grid.arrange(box_ele, box_lat, nrow = 2)

ggsave(pts_all, height = 5, width = 11, device = "png", path = "figures/corrected-shifts", 
       filename = "dp-vs-raw-shift_boxes_by-group_cent.png")


## make plot of dispersal potential vs. climate velocity coloured by lag within studies
# lags %>%
#   filter(Gradient == "Latitudinal") %>%
#   filter(Source %in% unique(.$Source)[1]) %>%
#   ggplot(aes(x = annual_dispersal_pot, y = lag, colour = climate_velocity)) + 
#   geom_point() +
#   theme_bw() +
#   scale_x_log10() +
#   geom_smooth(method = "lm") +
#   #geom_abline(slope = 1, intercept = 0) + 
#   facet_wrap(Source~Group)

lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Plants") %>%
  ggplot(aes(x = annual_dispersal_pot, y = CorrShift, colour = climate_velocity,
             group = climate_velocity
             )) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  geom_smooth(method = "lm", se = F) + 
  facet_wrap(~Position)+
  labs(colour = "Climate velocity")

lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Birds") %>%
  ggplot(aes(x = annual_dispersal_pot, y = CorrShift, colour = climate_velocity,
             group = climate_velocity
             )) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  geom_smooth(method = "lm", se = F)  + 
  facet_wrap(~Position)

lags %>%
  filter(Gradient == "Elevation") %>%
  filter(Group == "Plants") %>%
  ggplot(aes(x = annual_dispersal_pot, y = CorrShift, colour = climate_velocity,
             group = climate_velocity
             )) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  geom_smooth(method = "lm", se = F) + 
  facet_wrap(~Position)

lags %>%
  filter(Gradient == "Elevation") %>%
  filter(Group == "Birds") %>%
  ggplot(aes(x = annual_dispersal_pot, y = CorrShift, colour = climate_velocity,
             group = climate_velocity
             )) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  geom_smooth(method = "lm", se = F) + 
  facet_wrap(~Position)


## modelling approach:
## climate velocity as interaction with dispersal distance 
##    high climate velocity + high dispersal potential = high range shift 
##    low climate velocity + high dispersal potential = low range shift
##    high climate velocity + low dispersal potential = low range shift 
##    low climate velocity + low dispersal potential = low range shift
## birds and plants separate 







######################
####   modelling  ####
######################
## get rid of groups that have only species who can or only species who can't keep up with climate change 
# groups_keep <- lags %>% filter(Gradient == "Latitudinal") %>% group_by(expect_tracking, Group, Gradient) %>% tally() 
# lat_keep <- groups_keep$Group[which(duplicated(groups_keep$Group))]
# groups_keep <- lags %>% filter(Gradient == "Elevation")  %>% group_by(expect_tracking, Group, Gradient) %>% tally() 
# ele_keep <- groups_keep$Group[which(duplicated(groups_keep$Group))]
# 
# lat <- filter(lags, Group %in% lat_keep & Gradient == "Latitudinal")%>%
#   mutate(Group = factor(Group, ordered = F))
# ele <-  filter(lags, Group %in% ele_keep & Gradient == "Elevation") %>%
#   mutate(Group = factor(Group, ordered = F))

lat <- filter(lags,Gradient == "Latitudinal")%>%
  mutate(Group = factor(Group, ordered = F))
ele <-  filter(lags,Gradient == "Elevation") %>%
  mutate(Group = factor(Group, ordered = F))

## try modelling shift ~ dispersal potential*Group, while controlling for study ID
library(nlme)

## start with only lat x birds 
lat_birds <- filter(lat, Group == "Birds")
mod_lat = lme(data = lat_birds, ShiftR ~ log(annual_dispersal_pot)*Position, random = ~1|Source)
summary(mod_lat)

## plot predictions
lat_pred <- expand.grid(Position = unique(lat_birds$Position), 
                        annual_dispersal_pot = seq(min(lat_birds$annual_dispersal_pot), max(lat_birds$annual_dispersal_pot),
                                                   by = 1), 
                        Source = unique(lat_birds$Source))
pred <- predict(mod_lat, lat_pred, se.fit = T, re.form = NA)
lat_pred$predicted_val = pred

mod_data <- lat_pred %>%
  group_by(annual_dispersal_pot, Position) %>%
  mutate(predicted_val = mean(predicted_val, na.rm = T)) %>%
  ungroup() %>% 
  select(-Source) %>%
  distinct() 

mod_data %>%
  ggplot(data = ., aes(x = annual_dispersal_pot), y = predicted_val, 
                       colour = Position) +
  geom_line() +
  # geom_errorbar(aes(ymax = predicted_val + predicted_se,
  #               ymin = predicted_val - predicted_se), alpha = 0.03) +
  # facet_wrap(Position~Group) +
  theme_classic() +
  labs(x = "Maximum dispersal potential (km/y)",
       y = "Predicted range shift (km/y)",
       title = "Latitude") +
  geom_point(data = filter(lat, Group == "Birds"), aes(y = ShiftR, x = annual_dispersal_pot)) +
  scale_x_log10()


## now elev x birds
ele_birds <- filter(ele, Group == "Birds")
mod_ele = lme(data = ele_birds, ShiftR ~ log(annual_dispersal_pot)*Position, random = ~1|Source)
summary(mod_ele)

## plot predictions
ele_pred <- expand.grid(Position = unique(ele_birds$Position), 
                        annual_dispersal_pot = seq(min(ele_birds$annual_dispersal_pot), max(ele_birds$annual_dispersal_pot),
                                                   by = 5), 
                        Source = unique(ele_birds$Source))
pred <- predict(mod_ele, ele_pred, se.fit = T, re.form = NA)
ele_pred$predicted_val = pred

mod_data <- ele_pred %>%
  group_by(annual_dispersal_pot, Position) %>%
  mutate(predicted_val = mean(predicted_val, na.rm = T)) %>%
  ungroup() %>% 
  select(-Source) %>%
  distinct() 

mod_data %>%
  ggplot(data = ., aes(x = annual_dispersal_pot, y = predicted_val, 
                       colour = Position)) +
  geom_line() +
  # geom_errorbar(aes(ymax = predicted_val + predicted_se,
  #               ymin = predicted_val - predicted_se), alpha = 0.03) +
  # facet_wrap(Position~Group) +
  theme_classic() +
  labs(x = "Maximum dispersal potential (m/y)",
       y = "Predicted range shift (m/y)",
       title = "Elevation") +
  geom_point(data = filter(ele, Group == "Birds"), aes(y = ShiftR, x = annual_dispersal_pot)) +
  scale_x_log10()


##now with lat x plants 
lat_plants <- filter(lat, Group == "Plants")
mod_lat = lme(data = lat_plants, ShiftR ~ log(annual_dispersal_pot)*Position, random = ~1|Source)
summary(mod_lat)

## plot predictions
lat_pred <- expand.grid(Position = unique(lat_plants$Position), 
                        annual_dispersal_pot = seq(min(lat_plants$annual_dispersal_pot), max(lat_plants$annual_dispersal_pot),
                                                   by = 1), 
                        Source = unique(lat_plants$Source))
pred <- predict(mod_lat, lat_pred, se.fit = T, re.form = NA)
lat_pred$predicted_val = pred

mod_data <- lat_pred %>%
  group_by(annual_dispersal_pot, Position) %>%
  mutate(predicted_val = mean(predicted_val, na.rm = T)) %>%
  ungroup() %>% 
  select(-Source) %>%
  distinct() 

mod_data %>%
  ggplot(data = ., aes(x = annual_dispersal_pot, y = predicted_val, 
         colour = Position)) +
  geom_line() +
  # geom_errorbar(aes(ymax = predicted_val + predicted_se,
  #               ymin = predicted_val - predicted_se), alpha = 0.03) +
  # facet_wrap(Position~Group) +
  theme_classic() +
  labs(x = "Maximum dispersal potential (km/y)",
       y = "Predicted range shift (km/y)",
       title = "Latitude") +
  geom_point(data = filter(lat, Group == "Plants"), aes(y = ShiftR, x = annual_dispersal_pot)) +
  scale_x_log10()


## now elev x plants
ele_plants <- filter(ele, Group == "Plants")
mod_ele = lme(data = ele_plants, ShiftR ~ log(annual_dispersal_pot)*Position, random = ~1|Source)
summary(mod_ele)

## plot predictions
ele_pred <- expand.grid(Position = unique(ele_plants$Position), 
                        annual_dispersal_pot = seq(min(ele_plants$annual_dispersal_pot), max(ele_plants$annual_dispersal_pot),
                                                   by = 5), 
                        Source = unique(ele_plants$Source))
pred <- predict(mod_ele, ele_pred, se.fit = T, re.form = NA)
ele_pred$predicted_val = pred

mod_data <- ele_pred %>%
  group_by(annual_dispersal_pot, Position) %>%
  mutate(predicted_val = mean(predicted_val, na.rm = T)) %>%
  ungroup() %>% 
  select(-Source) %>%
  distinct() 

mod_data %>%
  ggplot(data = ., aes(x = annual_dispersal_pot, y = predicted_val, 
                       colour = Position)) +
  geom_line() +
  # geom_errorbar(aes(ymax = predicted_val + predicted_se,
  #               ymin = predicted_val - predicted_se), alpha = 0.03) +
  # facet_wrap(Position~Group) +
  theme_classic() +
  labs(x = "Maximum dispersal potential (m/y)",
       y = "Predicted range shift (m/y)",
       title = "Elevation") +
  geom_point(data = filter(ele, Group == "Plants"), aes(y = ShiftR, x = annual_dispersal_pot)) +
  scale_x_log10()


### try one model for everything
#######################
###### Latitude ######
#######################
mod_lat = lmer(data = lat, ShiftR ~ log(annual_dispersal_pot)*Position*Group + (1|Source))
mod_lat_tax = lmer(data = lat, ShiftR ~ log(annual_dispersal_pot)*Position*Group + (1|Source) + (1|family/Genus))

summary(mod_lat)
summary(mod_lat_tax)

R2m = r.squaredGLMM(mod_lat)[[1]] # marginal R2, ie part of variation explain by fixed effect 
R2c = r.squaredGLMM(mod_lat)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect 
R2m_tax = r.squaredGLMM(mod_lat_tax)[[1]] 
R2c_tax = r.squaredGLMM(mod_lat_tax)[[2]]

## plot predictions
lat_pred <- expand.grid(Group = unique(lat$Group),
                        Position = unique(lat$Position), 
                        annual_dispersal_pot = seq(min(lat$annual_dispersal_pot), max(lat$annual_dispersal_pot),
                                                   by = 1))
pred <- predict(mod_lat_tax, lat_pred, se.fit = T, re.form = NA)
lat_pred$predicted_val = pred

mod_data <- lat_pred %>%
  group_by(annual_dispersal_pot, Position, Group) %>%
  mutate(predicted_val = mean(predicted_val, na.rm = T)) %>%
  ungroup() %>% 
  distinct() 

## restrict to proper data range
minmax <- lat %>%
  group_by(Group, Position) %>%
  summarize(max = max(annual_dispersal_pot), 
            min = min(annual_dispersal_pot))

lat_plots <- mod_data %>%
  left_join(., minmax) %>%
  filter(annual_dispersal_pot >= min, annual_dispersal_pot <= max) %>%
  ggplot(data = ., aes(x = annual_dispersal_pot, y = predicted_val, 
                       colour = Group))  +
  # geom_errorbar(aes(ymax = predicted_val + predicted_se,
  #               ymin = predicted_val - predicted_se), alpha = 0.03) +
  # facet_wrap(Position~Group) +
  theme_classic() +
  labs(x = "Maximum dispersal potential (km/y)",
       y = "Predicted range shift (km/y)",
       title = "Latitude") +
  geom_point(data = lat, aes(y = ShiftR, x = annual_dispersal_pot, colour = Group),
             alpha = 0.1) +
  scale_x_log10() + 
  facet_wrap(~Position) +
  geom_line()

#######################
###### Elevation ######
#######################
mod_ele = lmer(data = ele, ShiftR ~ log(annual_dispersal_pot)*Position*Group + (1|Source))
mod_ele_tax = lmer(data = ele, ShiftR ~ log(annual_dispersal_pot)*Position*Group + (1|Source) + (1|family/Genus))

summary(mod_ele)
summary(mod_ele_tax)

R2m = r.squaredGLMM(mod_ele)[[1]] # marginal R2, ie part of variation explain by fixed effect 
R2c = r.squaredGLMM(mod_ele)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect 
R2m_tax = r.squaredGLMM(mod_ele_tax)[[1]] 
R2c_tax = r.squaredGLMM(mod_ele_tax)[[2]]


## plot predictions
ele_pred <- expand.grid(Group = unique(ele$Group),
                        Position = unique(ele$Position), 
                        annual_dispersal_pot = seq(min(ele$annual_dispersal_pot), max(ele$annual_dispersal_pot),
                                                   by = 1))
pred <- predict(mod_ele_tax, ele_pred, se.fit = T, re.form = NA)
ele_pred$predicted_val = pred

mod_data <- ele_pred %>%
  group_by(annual_dispersal_pot, Position, Group) %>%
  mutate(predicted_val = mean(predicted_val, na.rm = T)) %>%
  ungroup() %>% 
  distinct() 

## restrict to proper data range
minmax <- ele %>%
  group_by(Group, Position) %>%
  summarize(max = max(annual_dispersal_pot), 
            min = min(annual_dispersal_pot))

ele_plots <- mod_data %>%
  left_join(., minmax) %>%
  filter(annual_dispersal_pot >= min, annual_dispersal_pot <= max) %>%
  ggplot(data = ., aes(x = annual_dispersal_pot, y = predicted_val, 
                       colour = Group))  +
  # geom_errorbar(aes(ymax = predicted_val + predicted_se,
  #               ymin = predicted_val - predicted_se), alpha = 0.03) +
  # facet_wrap(Position~Group) +
  theme_classic() +
  labs(x = "Maximum dispersal potential (m/y)",
       y = "Predicted range shift (m/y)",
       title = "Elevation") +
  geom_point(data = ele, aes(y = ShiftR, x = annual_dispersal_pot, colour = Group),
             alpha = 0.1) +
  scale_x_log10() + 
  facet_wrap(~Position) +
  geom_line()




## try other models that control for climate velocity 
#######################
### Latitude (cv) #####
#######################
mod_lat_cv = lmer(data = lat, ShiftR ~ log(annual_dispersal_pot)*Position*Group + climate_velocity +
                     (1|Source))
mod_lat_tax_cv = lmer(data = lat, ShiftR ~ log(annual_dispersal_pot)*Position*Group + climate_velocity +
                     (1|Source) + (1|family/Genus))

summary(mod_lat_cv)
summary(mod_lat_tax_cv)

R2m_cv = r.squaredGLMM(mod_lat_cv)[[1]] # marginal R2, ie part of variation explain by fixed effect 
R2c_cv = r.squaredGLMM(mod_lat_cv)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect 
R2m_tax_cv = r.squaredGLMM(mod_lat_cv)[[1]] 
R2c_tax_cv = r.squaredGLMM(mod_lat_tax_cv)[[2]]
## this model has better R2s 

## plot predictions
lat_pred <- expand.grid(Group = unique(lat$Group),
                        Position = unique(lat$Position), 
                        annual_dispersal_pot = seq(min(lat$annual_dispersal_pot), max(lat$annual_dispersal_pot),
                                                   by = 1),
                        climate_velocity = mean(lat$climate_velocity))
pred <- predict(mod_lat_tax_cv, lat_pred, se.fit = T, re.form = NA)
lat_pred$predicted_val = pred

mod_data <- lat_pred %>%
  group_by(annual_dispersal_pot, Position, Group) %>%
  mutate(predicted_val = mean(predicted_val, na.rm = T)) %>%
  ungroup() %>% 
  distinct() 

## restrict to proper data range
minmax <- lat %>%
  group_by(Group, Position) %>%
  summarize(max = max(annual_dispersal_pot), 
            min = min(annual_dispersal_pot))

mod_data %>%
  left_join(., minmax) %>%
  filter(annual_dispersal_pot >= min, annual_dispersal_pot <= max) %>%
  ggplot(data = ., aes(x = annual_dispersal_pot, y = predicted_val, 
                       colour = Group))  +
  # geom_errorbar(aes(ymax = predicted_val + predicted_se,
  #               ymin = predicted_val - predicted_se), alpha = 0.03) +
  # facet_wrap(Position~Group) +
  theme_classic() +
  labs(x = "Maximum dispersal potential (km/y)",
       y = "Predicted range shift (km/y)",
       title = "Latitude") +
  geom_point(data = lat, aes(y = ShiftR, x = annual_dispersal_pot, colour = Group),
             alpha = 0.1) +
  scale_x_log10() + 
  facet_wrap(~Position) +
  geom_line()


## plot predictions across climate velocity
lat_pred <- expand.grid(Group = unique(lat$Group),
                        Position = unique(lat$Position), 
                        annual_dispersal_pot = mean(lat$annual_dispersal_pot), 
                        climate_velocity = seq(min(lat$climate_velocity), max(lat$climate_velocity),
                                               by = 0.1))
pred <- predict(mod_lat_tax_cv, lat_pred, se.fit = T, re.form = NA)
lat_pred$predicted_val = pred

mod_data <- lat_pred %>%
  group_by(climate_velocity, Position, Group) %>%
  mutate(predicted_val = mean(predicted_val, na.rm = T)) %>%
  ungroup() %>% 
  distinct() 

## restrict to proper data range
minmax <- lat %>%
  group_by(Group, Position) %>%
  summarize(max = max(climate_velocity), 
            min = min(climate_velocity))

mod_data %>%
  left_join(., minmax) %>%
  filter(climate_velocity >= min, climate_velocity <= max) %>%
  ggplot(data = ., aes(x = climate_velocity, y = predicted_val, 
                       colour = Group)) +
  # geom_errorbar(aes(ymax = predicted_val + predicted_se,
  #               ymin = predicted_val - predicted_se), alpha = 0.03) +
  # facet_wrap(Position~Group) +
  theme_classic() +
  labs(x = "Climate velocity (km/y)",
       y = "Predicted range shift (km/y)",
       title = "Latitude") +
  geom_point(data = lat, aes(y = ShiftR, x = climate_velocity, colour = Group),
             alpha = 0.1) +
  #scale_x_log10() + 
  facet_wrap(~Position) +
  geom_line()


#######################
### Elevation (cv) ####
#######################
mod_ele_cv = lmer(data = ele, ShiftR ~ log(annual_dispersal_pot)*Position*Group + climate_velocity +
                    (1|Source))
mod_ele_tax_cv = lmer(data = ele, ShiftR ~ log(annual_dispersal_pot)*Position*Group + climate_velocity +
                        (1|Source) + (1|family/Genus))

summary(mod_ele_cv)
summary(mod_ele_tax_cv)

R2m_cv = r.squaredGLMM(mod_ele_cv)[[1]] # marginal R2, ie part of variation explain by fixed effect 
R2c_cv = r.squaredGLMM(mod_ele_cv)[[2]] # conditional R2, ie part of variation explain by fixed effect and random effect 
R2m_tax_cv = r.squaredGLMM(mod_ele_cv)[[1]] 
R2c_tax_cv = r.squaredGLMM(mod_ele_tax_cv)[[2]]
## this model has better R2s but by a bit

## plot predictions
ele_pred <- expand.grid(Group = unique(ele$Group),
                        Position = unique(ele$Position), 
                        annual_dispersal_pot = seq(min(ele$annual_dispersal_pot), max(ele$annual_dispersal_pot),
                                                   by = 5), 
                        climate_velocity = mean(ele$climate_velocity))
pred <- predict(mod_ele_cv, ele_pred, se.fit = T, re.form = NA)
ele_pred$predicted_val = pred

mod_data <- ele_pred %>%
  group_by(annual_dispersal_pot, Position, Group) %>%
  mutate(predicted_val = mean(predicted_val, na.rm = T)) %>%
  ungroup() %>% 
  distinct() 

## restrict to proper data range
minmax <- ele %>%
  group_by(Group, Position) %>%
  summarize(max = max(annual_dispersal_pot), 
            min = min(annual_dispersal_pot))

mod_data %>%
  left_join(., minmax) %>%
  filter(annual_dispersal_pot >= min, annual_dispersal_pot <= max) %>%
  ggplot(data = ., aes(x = annual_dispersal_pot, y = predicted_val, 
                       colour = Group))  +
  # geom_errorbar(aes(ymax = predicted_val + predicted_se,
  #               ymin = predicted_val - predicted_se), alpha = 0.03) +
  # facet_wrap(Position~Group) +
  theme_classic() +
  labs(x = "Maximum dispersal potential (km/y)",
       y = "Predicted range shift (km/y)",
       title = "Elevation") +
  geom_point(data = ele, aes(y = ShiftR, x = annual_dispersal_pot, colour = Group),
             alpha = 0.1) +
  scale_x_log10() + 
  facet_wrap(~Position) +
  geom_line()


## plot predictions across climate velocity 
ele_pred <- expand.grid(Group = unique(ele$Group),
                        Position = unique(ele$Position), 
                        annual_dispersal_pot = mean(ele$annual_dispersal_pot), 
                        climate_velocity = seq(min(ele$climate_velocity), max(ele$climate_velocity),
                                               by = 0.1))
pred <- predict(mod_ele_cv, ele_pred, se.fit = T, re.form = NA)
ele_pred$predicted_val = pred

mod_data <- ele_pred %>%
  group_by(climate_velocity, Position, Group) %>%
  mutate(predicted_val = mean(predicted_val, na.rm = T)) %>%
  ungroup() %>% 
  distinct() 

## restrict to proper data range
minmax <- ele %>%
  group_by(Group, Position) %>%
  summarize(max = max(climate_velocity), 
            min = min(climate_velocity))

mod_data %>%
  left_join(., minmax) %>%
  filter(climate_velocity >= min, climate_velocity <= max) %>%
  ggplot(data = ., aes(x = climate_velocity, y = predicted_val, 
                       colour = Group)) +
  # geom_errorbar(aes(ymax = predicted_val + predicted_se,
  #               ymin = predicted_val - predicted_se), alpha = 0.03) +
  # facet_wrap(Position~Group) +
  theme_classic() +
  labs(x = "Climate velocity (km/y)",
       y = "Predicted range shift (km/y)",
       title = "Elevation") +
  geom_point(data = ele, aes(y = ShiftR, x = climate_velocity, colour = Group),
             alpha = 0.1) +
  #scale_x_log10() + 
  facet_wrap(~Position) +
  geom_line()


## newest plan: compete models
## one set for latitude, one set for elevation
## source  and nested taxonomy as random effect 
## with/without climate velocity 
## with/without different methods 

## problem:
## lmer doesn't report AIC 








#######################
####    garbage   #####
#######################
## now model
# lat_mod <- lm(lag ~ expect_tracking*Group*Position, data = lat)
# summary(lat_mod)
# plot(lat_mod)


# ncvTest(lat_mod)
## residual error does not have constant variance
## zero bound data - Poisson?
## ask model to allow variance to get smaller depending on value of y 

# ele_mod <- lm(lag ~ expect_tracking*Position, data = ele)
# summary(ele_mod)
# plot(ele_mod)

## try new model
lat$raw_lag <- as.numeric(as.character(lat$raw_lag))
lat_mod <- lm(raw_lag ~ dispersal_potential_kmY*Group, data = lat)
summary(lat_mod)
#plot(lat_mod)

ele_mod <- lm(raw_lag ~ dispersal_potential_mY*Position, data = ele)
summary(ele_mod)
#plot(ele_mod)

## now plot predictions 
# lat_pred <- expand.grid(Group = unique(lat$Group), expect_tracking = unique(lat$expect_tracking),
#                         Position = unique(ele$Position)) 
# pred <- predict(lat_mod, lat_pred, se.fit = TRUE)
# lat_pred$predicted_val = pred$fit
# lat_pred$predicted_se = pred$se.fit
# lat_pred$expect_tracking = factor(lat_pred$expect_tracking, levels = c("No", "Yes"),
#                                   ordered = T)
# 
# ele_pred <- expand.grid(Group = unique(ele$Group), expect_tracking = unique(ele$expect_tracking),
#                         Position = unique(ele$Position)) 
# pred <- predict(ele_mod, ele_pred, se.fit = TRUE)
# ele_pred$predicted_val = pred$fit
# ele_pred$predicted_se = pred$se.fit
# ele_pred$expect_tracking = factor(ele_pred$expect_tracking, levels = c("No", "Yes"),
#                                   ordered = T)
# 
# predplot_l <- ggplot(data = lat_pred, aes(x = expect_tracking, colour = expect_tracking, y = predicted_val,
#                                           shape = Position)) + 
#   geom_point() +
#   geom_errorbar(aes(ymax = predicted_val + predicted_se,
#                 ymin = predicted_val - predicted_se)) +
#   facet_wrap(~Group) +
#   theme_classic() +
#   labs(x = "Expect dispersal potential to allow species to keep up with temp change?",
#        y = "Model-predicted inferred range shift lag (km/y)",
#        title = "Latitude") +
#   scale_colour_manual(values = c(pal[6], pal[4])) +
#   theme(legend.position = "none")
# 
# predplot_e <- ggplot(data = ele_pred, aes(x = expect_tracking, colour = expect_tracking, y = predicted_val,
#                                           shape = Position)) + 
#   geom_point() +
#   geom_errorbar(aes(ymax = predicted_val + predicted_se,
#                     ymin = predicted_val - predicted_se)) +
#   facet_wrap(~Group) +
#   theme_classic() +
#   labs(x = "", title = "Elevation",
#        y = "Model-predicted inferred range shift lag (m/y)") +
#   scale_colour_manual(values = c(pal[6], pal[4])) +
#   guides(colour = "none")
# 
# gr <- grid.arrange(predplot_l, predplot_e, nrow = 1, widths = c(0.6, 0.4))
# 
# ggsave(gr, height = 4, width = 9, device = "png", path = "figures/dispersal", 
#        filename = "pred_inferred.png")



lat_pred <- expand.grid(#Position = unique(lat$Position), 
  Group = unique(lat$Group), 
  dispersal_potential_kmY = seq(min(lat$dispersal_potential_kmY), max(lat$dispersal_potential_kmY),
                                by = 1))
pred <- predict(lat_mod, lat_pred, se.fit = TRUE)
lat_pred$predicted_val = pred$fit
lat_pred$predicted_se = pred$se.fit

ele_pred <- expand.grid(Position = unique(ele$Position), 
                        dispersal_potential_mY = seq(min(ele$dispersal_potential_mY), max(ele$dispersal_potential_mY),
                                                     by = 1))
pred <- predict(ele_mod, ele_pred, se.fit = TRUE)
ele_pred$predicted_val = pred$fit
ele_pred$predicted_se = pred$se.fit

ggplot(data = lat_pred, aes(x = dispersal_potential_kmY, y = predicted_val,
                            colour = Group)) +
  geom_point() +
  # geom_errorbar(aes(ymax = predicted_val + predicted_se,
  #               ymin = predicted_val - predicted_se), alpha = 0.03) +
  #facet_wrap(Position~Group) +
  theme_classic() +
  labs(x = "Maximum dispersal potential (km/y)",
       y = "Predicted range shift lag (km/y)",
       title = "Latitude") 

ggplot(data = ele_pred, aes(x = expect_tracking, colour = expect_tracking, y = predicted_val,
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


lags %>%
  ggplot(aes(x = annual_dispersal_pot, y = CorrShift, colour = climate_velocity)) + 
  geom_point() +
  theme_bw() +
  scale_x_log10() +
  facet_grid(~Gradient)

lags %>%
  filter(expect_tracking == "No") %>%
  ggplot(aes(x = annual_dispersal_pot, y = CorrShift, colour = Group)) + 
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
  ggplot(., aes(x = CorrShift, fill = Source)) + geom_histogram() 

lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Birds") %>%
  ggplot(., aes(x = lag, fill = Source)) + geom_histogram() 

lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Birds") %>%
  ggplot(., aes(x = CorrShift, fill = Source)) + geom_histogram() 

lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Fish") %>%
  ggplot(., aes(x = lag, fill = Source)) + geom_histogram() 

lags %>%
  filter(Gradient == "Latitudinal") %>%
  filter(Group == "Fish") %>%
  ggplot(., aes(x = CorrShift, fill = Source)) + geom_histogram() 

lags %>%
  filter(Gradient == "Elevation") %>%
  filter(Group == "Plants") %>%
  ggplot(aes(x = lag, fill = Source)) + geom_histogram()

lags %>%
  filter(Gradient == "Elevation") %>%
  filter(Group == "Plants") %>%
  ggplot(aes(x = CorrShift, fill = Source)) + geom_histogram()
