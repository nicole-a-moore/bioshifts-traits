## seeing whether empirical dispersal scale explains variation in range shifts
library(tidyverse)
library(PNWColors)
library(gridExtra)
library(grid)
library(cowplot)
library(lme4)
source("R/taxonomic-harmonization/clean_taxa_functions.R")


#############################
####   data preparation  ####
#############################
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


## read in bioshifts v2
v2 = read.csv("data-raw/bioshiftsv2/Bioshifts.v2.final.csv")
## clean names to make them match reported names in species list
v2$reported_name = Clean_Names(gsub("_", " ", v2$Scientific.Name), return_gen_sps = F)

spv2 <- filter(sp, v2 == 1) %>%
  select(reported_name, scientificName)

## yay! all are there
## except Coprosma 3sp
which(!v2$reported_name %in% spv2$reported_name)

## join to add scientific name column
v2 = left_join(v2, spv2)

## subset to species with dispersal scale 
v2 <- filter(v2, scientificName %in% dscale$scientificName)
length(unique(v2$scientificName)) #587 species 

## make sure all the species are there:
length(unique(c(v2$scientificName, v1$scientificName))) # 651 unique species 
length(unique(dscale$scientificName)) # 651 unique species 

#----------------------
### relate velocity of shift to dispersal scale ###
colnames(v1)
colnames(v2)


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

## 593 species have max dispersal distance
## start with max!
v1 <- v1 %>%
  filter(Code == "MaxDispersalDistance") %>%
  select(-DispersalDistanceKm, -DispersalDistance, -DispersalUnit, -Field, 
         -DispersalSource, -Database, -ObservationType, -Sex) %>%
  unique()

which(nomaxsp %in% v1$Species)
which(maxsp %in% v1$Species)

v1 <- rename(v1, "MaxDispersalDistanceKm" = DispersalDistanceKm_unique)

## save dataset 
#write.csv(v1, "data-processed/bioshiftsv1_max-dispersal-distance.csv", row.names = FALSE)
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

#----------------------
#take log of dispersal scale and dispersal potential 
v1$MaxDispersalDistanceKm = log(v1$MaxDispersalDistanceKm)
hist(v1$MaxDispersalDistanceKm)
v1$MaxDispersalPotentialKmY = log(v1$MaxDispersalPotentialKmY)
hist(v1$MaxDispersalPotentialKmY)
v1$MaxDispersalPotentialmY = log(v1$MaxDispersalPotentialmY)
v1$MaxDispersalDistancem = log(v1$MaxDispersalDistancem)

v1 %>%
  ggplot(aes(x = MaxDispersalDistanceKm, y = MaxDispersalPotentialKmY, colour = class)) + geom_point() +
  theme_bw()

#----------------------
#prepare a dataset for each gradient

#elevation
ele <- v1[v1$Gradient == "Elevation",]
ele <- ele[which(is.na(ele$EleVeloT)==F),]
ele$lags <- ele$ShiftR - ele$EleVeloT
hist(ele$lags)

#latitude
lat <- v1[v1$Gradient == "Latitudinal",]
lat$lags <- lat$ShiftR - lat$LatVeloT
hist(lat$lags)

## how many shifts for species with maximum dispersal distance in v1? 
nrow(ele) + nrow(lat) #5284

## how many leading edge shifts for species with maximum dispersal distance in v1? 
length(which(ele$Position == "Leading edge")) + length(which(lat$Position == "Leading edge")) #1393

lat = lat %>%
  filter(LatVeloT >= 0) %>%
  ## get rid of cases where you expect contraction at leading edge (negative velocity)
  mutate(expect_tracking_dist = ifelse(MaxDispersalDistanceKm > log(LatVeloT),"Yes", "No")) %>%
  mutate(expect_tracking_pot = ifelse(MaxDispersalPotentialKmY > log(LatVeloT),"Yes", "No")) %>%
  ## get rid of trailing edge 
  filter(., Position != "Trailing edge") %>%
  filter(!is.na(MaxDispersalPotentialKmY))

ele = ele %>%
  filter(EleVeloT >= 0) %>%
  ## get rid of cases where you expect contraction at leading edge (negative velocity)
  mutate(expect_tracking_dist = ifelse(MaxDispersalDistancem > log(EleVeloT),"Yes", "No")) %>%
  mutate(expect_tracking_pot = ifelse(MaxDispersalPotentialmY > log(EleVeloT),"Yes", "No"))  %>%
  ## get rid of trailing edge 
  filter(., Position != "Trailing edge") %>%
  filter(!is.na(MaxDispersalPotentialKmY))


lat %>% rbind(., ele) %>%
  ggplot(aes(x = Position)) + geom_bar() + theme_bw() +
  labs(y = "Count", x = "Range shift parameter") + 
  facet_wrap(~Gradient)

lat %>% rbind(., ele) %>%
  group_by(Gradient, Position) %>%
  tally()

lat %>% rbind(., ele) %>%
  select(Species, class) %>%
  unique() %>%
  tally()

###############################
####   data visualization  ####
###############################
#----------------------
#prepare a dataset for each range shift parameter
ele_te <- ele[ele$Position=="Trailing edge",]
ele_le <- ele[ele$Position=="Leading edge",]
ele_o <- ele[ele$Position=="Centroid",]
lat_te <- lat[lat$Position=="Trailing edge",]
lat_le <- lat[lat$Position=="Leading edge",]
lat_o <- lat[lat$Position=="Centroid",]

pal = pnw_palette("Bay",7)


## how many species should be able to keep up with climate change across elev and lat?
pol1 <- data.frame(x = c(-10, 10, 10), y = c(-10, -10, 10))
pol2 <- data.frame(x = c(-10, -10, 10), y = c(-10, 10, 10))
pol3 <- data.frame(x = c(-4, 15, 15), y = c(-4, -4, 15))
pol4 <- data.frame(x = c(-4, -4, 15), y = c(-4, 15, 15))

plot1 = lat %>% 
  ggplot(aes(y = MaxDispersalDistanceKm, x = log(LatVeloT))) + 
  geom_polygon(data = pol1, aes(x = x, y = y), fill = pal[6], alpha = 0.5) +
  geom_polygon(data = pol2, aes(x = x, y = y), fill = pal[4], alpha = 0.5) +
  geom_point(aes(colour = expect_tracking_dist)) +
  scale_x_continuous(limits = c(-10, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(-10, 10), expand = c(0,0)) +
  theme_light() + theme(panel.grid = element_blank()) + 
  labs(y = "Log maximum dispersal distance (km)", 
       x = "Log mean velocity of temperature (km/y)", 
       title = "Latitudinal shifts")+
  theme_classic() +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme(legend.position = "none")

plot2 = ele %>% 
  ggplot(aes(y = MaxDispersalDistancem, x = log(EleVeloT))) +
  geom_polygon(data = pol3, aes(x = x, y = y), fill = pal[6], alpha = 0.5) +
  geom_polygon(data = pol4, aes(x = x, y = y), fill = pal[4], alpha = 0.5) +
  geom_point(aes(colour = expect_tracking_dist)) +
  scale_x_continuous(limits = c(-4, 15), expand = c(0,0)) +
  scale_y_continuous(limits = c(-4, 15), expand = c(0,0)) +
  theme(panel.grid = element_blank())  +
  theme_classic() + 
  labs(colour = "Does dispersal\npotential allow\nspecies to keep up\nwith temp change?",
       y = "Log maximum dispersal distance (m)", 
       x = "Log mean velocity of temperature (m/y)", 
       title = "Elevational shifts")  +
  scale_colour_manual(values = c(pal[6], pal[4]))

p = plot_grid(plot1, plot2, nrow = 1, rel_widths = c(4/10, 6/10))
## more species should be able to keep up climate change across elevation 

ggsave(p, height = 3.5, width = 9, device = "png", path = "figures/dispersal", 
       filename = "maxdispersal-versus-inferred-lag.png")

## try with dispersal potential 
plot3 = lat %>% 
  ggplot(aes(y = MaxDispersalPotentialKmY, x = log(LatVeloT))) + 
  geom_polygon(data = pol1, aes(x = x, y = y), fill = pal[6], alpha = 0.5) +
  geom_polygon(data = pol2, aes(x = x, y = y), fill = pal[4], alpha = 0.5) +
  geom_point(aes(colour = expect_tracking_pot)) +
  scale_x_continuous(limits = c(-10, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(-10, 10), expand = c(0,0)) +
  theme_light() + theme(panel.grid = element_blank()) + 
  labs(y = "Log maximum dispersal potential (km/y)", 
       x = "Log mean velocity of temperature (km/y)", 
       title = "Latitudinal shifts")+
  theme_classic() +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme(legend.position = "none")

plot4 = ele %>% 
  ggplot(aes(y = MaxDispersalPotentialmY, x = log(EleVeloT))) +
  geom_polygon(data = pol3, aes(x = x, y = y), fill = pal[6], alpha = 0.5) +
  geom_polygon(data = pol4, aes(x = x, y = y), fill = pal[4], alpha = 0.5) +
  geom_point(aes(colour = expect_tracking_pot)) +  
  scale_x_continuous(limits = c(-4, 15), expand = c(0,0)) +
  scale_y_continuous(limits = c(-4, 15), expand = c(0,0)) +
  theme(panel.grid = element_blank())  +
  theme_classic() + 
  labs(colour = "Does dispersal\npotential allow\nspecies to keep up\nwith temp change?",
       y = "Log maximum dispersal potential (m/y)", 
       x = "Log mean velocity of temperature (m/y)", 
       title = "Elevational shifts")  +
  scale_colour_manual(values = c(pal[6], pal[4]))

p2 = plot_grid(plot3, plot4, nrow = 1, rel_widths = c(4/10, 6/10))

ggsave(p2, height = 3.5, width = 9, device = "png", path = "figures/dispersal", 
       filename = "maxpotential-versus-inferred-lag.png")

ele %>%
  group_by(expect_tracking_pot) %>%
  tally()

lat %>%
  group_by(expect_tracking_pot) %>%
  tally()

lat %>%
  rbind(., ele) %>%
  ggplot(aes(x = MaxDispersalDistanceKm, fill = class)) + geom_histogram() +
  theme_bw() + facet_wrap(~Gradient)

## take a closer look 
ele_le %>%
  filter(abs(lags) < 30) %>% ## get rid of outliers
  ggplot(aes(x = MaxDispersalDistancem, y = lags, colour = expect_tracking_dist, shape = Sampling)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  facet_grid(~class) +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Elevational shifts",
       x = "Log maximum dispersal distance (m)", 
       y = "Range shift lag (m/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-3, 14)) +
  scale_y_continuous(limits = c(-24, 30))

ggsave(height = 3.5, width = 9, device = "png", path = "figures/dispersal", 
       filename = "ele-le_points_distance_inferred.png")

ele_le %>%
  filter(abs(lags) < 30) %>% ## get rid of outliers
  ggplot(aes(x = MaxDispersalDistancem, y = lags, colour = expect_tracking_dist)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Elevational shifts",
       x = "Log maximum dispersal distance (m)", 
       y = "Range shift lag (m/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  facet_wrap(~class, nrow = 1) +
  scale_x_continuous(limits = c(-3, 14)) +
  scale_y_continuous(limits = c(-24, 30))

ggsave(height = 3.5, width = 9, device = "png", path = "figures/dispersal", 
       filename = "ele-le_boxplot_distance_inferred.png")

lat_le %>%
  ggplot(aes(x = MaxDispersalDistanceKm, y = lags, colour = expect_tracking_dist, shape = Sampling)) + 
  geom_hline(yintercept = 0) +
  geom_point() +
  facet_grid(~class) +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Latitudinal shifts",
       x = "Log maximum dispersal distance (km)", 
       y = "Range shift lag (km/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-10, 9.5)) +
  scale_y_continuous(limits = c(-17, 34))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "lat-le_points_distance_inferred.png")

lat_le %>%
 ggplot(aes(x = MaxDispersalDistanceKm, y = lags, colour = expect_tracking_dist)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Latitudinal shifts",
       x = "Log maximum dispersal distance (km)", 
       y = "Range shift lag (km/y)") +
  facet_grid(~class) +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-10, 9.5)) +
  scale_y_continuous(limits = c(-17, 34))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "lat-le_boxplot_distance_inferred.png")

lat_o %>%
  ggplot(aes(x = MaxDispersalDistanceKm, y = lags, colour = expect_tracking_dist, shape = Sampling)) + 
  geom_hline(yintercept = 0) +
  geom_point() +
  facet_grid(~class) +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Latitudinal shifts",
       x = "Log maximum dispersal distance (km)", 
       y = "Range shift lag (km/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_continuous(limits = c(-20, 40))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "lat-o_points_distance_inferred.png")

lat_o %>%
  ggplot(aes(x = MaxDispersalDistanceKm, y = lags, colour = expect_tracking_dist)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Latitudinal shifts",
       x = "Log maximum dispersal distance (km)", 
       y = "Range shift lag (km/y)") +
  facet_grid(~class) +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_continuous(limits = c(-20, 40))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "lat-o_boxplot_distance_inferred.png")


ele_o %>%
  ggplot(aes(x = MaxDispersalDistancem, y = lags, colour = expect_tracking_dist, shape = Sampling)) + 
  geom_hline(yintercept = 0) +
  geom_point() +
  facet_grid(~class) +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Elevational shifts",
       x = "Log maximum dispersal distance (m)", 
       y = "Range shift lag (m/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-4, 14.5)) +
  scale_y_continuous(limits = c(-30, 20))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "ele-o_points_distance_inferred.png")

ele_o %>%
  ggplot(aes(x = MaxDispersalDistancem, y = lags, colour = expect_tracking_dist)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Latitudinal shifts",
       x = "Log maximum dispersal distance (m)", 
       y = "Range shift lag (m/y)") +
  facet_grid(~class) +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-4, 14.5)) +
  scale_y_continuous(limits = c(-30, 20))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "ele_o_boxplot_distance_inferred.png")

## now with dispersal potential 
ele_le %>%
  filter(abs(lags) < 30) %>% ## get rid of outliers
  ggplot(aes(x = MaxDispersalPotentialmY, y = lags, colour = expect_tracking_pot, shape = Sampling)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  facet_grid(~class) +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Elevational shifts",
       x = "Log maximum dispersal potential (m/y)", 
       y = "Range shift lag (m/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank())  +
  scale_x_continuous(limits = c(-3, 14)) +
  scale_y_continuous(limits = c(-24, 30))

ggsave(height = 3.5, width = 9, device = "png", path = "figures/dispersal", 
       filename = "ele-le_points_potential_inferred.png")

ele_le %>%
  filter(abs(lags) < 30) %>% ## get rid of outliers
  ggplot(aes(x = MaxDispersalPotentialmY, y = lags, colour = expect_tracking_pot)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Elevational shifts",
       x = "Log maximum dispersal potential (m/y)", 
       y = "Range shift lag (m/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  facet_wrap(~class, nrow = 1)  +
  scale_x_continuous(limits = c(-3, 14)) +
  scale_y_continuous(limits = c(-24, 30))

ggsave(height = 3.5, width = 9, device = "png", path = "figures/dispersal", 
       filename = "ele-le_boxplot_potential_inferred.png")

lat_le %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, y = lags, colour = expect_tracking_pot, shape = Sampling)) + 
  geom_hline(yintercept = 0) +
  geom_point() +
  facet_grid(~class) +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Latitudinal shifts",
       x = "Log maximum dispersal potential (km/y)", 
       y = "Range shift lag (km/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-10, 9.5)) +
  scale_y_continuous(limits = c(-17, 34))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "lat-le_points_potential_inferred.png")

lat_le %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, y = lags, colour = expect_tracking_pot)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Latitudinal shifts",
       x = "Log maximum dispersal potential (km/y)", 
       y = "Range shift lag (km/y)") +
  facet_grid(~class) +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank())+
  scale_x_continuous(limits = c(-10, 9.5)) +
  scale_y_continuous(limits = c(-17, 34))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "lat-le_boxplot_potential_inferred.png")

lat_o %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, y = lags, colour = expect_tracking_pot, shape = Sampling)) + 
  geom_hline(yintercept = 0) +
  geom_point() +
  facet_grid(~class) +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Latitudinal shifts",
       x = "Log maximum dispersal potential (km/y)", 
       y = "Range shift lag (km/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_continuous(limits = c(-20, 40))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "lat-o_points_potential_inferred.png")

lat_o %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, y = lags, colour = expect_tracking_pot)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Latitudinal shifts",
       x = "Log maximum dispersal potential (km/y)", 
       y = "Range shift lag (km/y)") +
  facet_grid(~class) +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_continuous(limits = c(-20, 40))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "lat-o_boxplot_potential_inferred.png")

ele_o %>%
  ggplot(aes(x = MaxDispersalPotentialmY, y = lags, colour = expect_tracking_pot, shape = Sampling)) + 
  geom_hline(yintercept = 0) +
  geom_point() +
  facet_grid(~class) +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Elevational shifts",
       x = "Log maximum dispersal potential (m/y)", 
       y = "Range shift lag (m/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank())+
  scale_x_continuous(limits = c(-4, 14.5)) +
  scale_y_continuous(limits = c(-30, 20))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "ele-o_points_potential_inferred.png")

ele_o %>%
  ggplot(aes(x = MaxDispersalPotentialmY, y = lags, colour = expect_tracking_pot)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Elevational shifts",
       x = "Log maximum dispersal potential (m/y)", 
       y = "Range shift lag (m/y)") +
  facet_grid(~class) +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-4, 14.5)) +
  scale_y_continuous(limits = c(-30, 20))

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "ele_o_boxplot_potential_inferred.png")

## so interesting!


## now plot inferred lags 

## get rid of species we could not estimate a corrected shift for 
lat <- filter(lat, !is.na(corrected_shift))
ele <- filter(ele, !is.na(corrected_shift))

## calculate inferred lag
lat$inferred_lag =  lat$corrected_shift - lat$LatVeloT
ele$inferred_lag =  ele$corrected_shift - ele$EleVeloT
 
#prepare a dataset for each range shift parameter
ele_te <- ele[ele$Position=="Trailing edge",]
ele_le <- ele[ele$Position=="Leading edge",]
ele_o <- ele[ele$Position=="Centroid",]
lat_te <- lat[lat$Position=="Trailing edge",]
lat_le <- lat[lat$Position=="Leading edge",]
lat_o <- lat[lat$Position=="Centroid",]

ele_le %>%
  ggplot(aes(x = MaxDispersalPotentialmY, y = inferred_lag, colour = expect_tracking_dist)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Elevational shifts",
       x = "Log maximum dispersal potential (m/y)", 
       y = "Inferred range shift lag (m/y)") +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) +
  facet_wrap(~class, nrow = 1)

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "ele_le_boxplot_potential_inferred_lag.png")

lat_le %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, y = inferred_lag, colour = expect_tracking_pot)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Latitudinal shifts",
       x = "Log maximum dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)") +
  facet_grid(~class) +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) 

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "lat_le_boxplot_potential_inferred_lag.png")

lat_o %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, y = inferred_lag, colour = expect_tracking_dist)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Latitudinal shifts",
       x = "Log maximum dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)") +
  facet_grid(~class) +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank())

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "lat_o_boxplot_potential_inferred_lag.png")

ele_o %>%
  ggplot(aes(x = MaxDispersalPotentialmY, y = inferred_lag, colour = expect_tracking_pot)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       title = "Elevational shifts",
       x = "Log maximum dispersal potential (m/y)", 
       y = "Inferred range shift lag (m/y)") +
  facet_grid(~class) +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme_bw() + 
  guides(colour = "legend", shape = "none") +
  theme(legend.position = "bottom", panel.grid = element_blank()) 

ggsave(height = 3.5, width = 11, device = "png", path = "figures/dispersal", 
       filename = "ele_o_boxplot_potential_inferred_lag.png")

###############################
####     model fitting     ####
###############################
## first: are range shift lags are larger in species whose dispersal potential is less than the speed of climate change?

## filter out classes with fewer than 5 sp that have dispersal potential
## elevation:
dat.MaxDispPot.ele <- data.frame(ele[which(is.na(ele$MaxDispersalPotentialmY)==F),])

dsp <- unique(dat.MaxDispPot.ele[,c(which(colnames(dat.MaxDispPot.ele) == "class"),
                                     which(colnames(dat.MaxDispPot.ele) == "Species"))])
cl <- table(dsp$class)
nn <- names(which(cl >= 5))
dat.MaxDispPot.ele <- dat.MaxDispPot.ele[which(dat.MaxDispPot.ele$class %in% nn),]
dat.MaxDispPot.ele <-droplevels(dat.MaxDispPot.ele)
length(unique(dat.MaxDispPot.ele$class)) # 6

## latitude:
dat.MaxDispPot.lat <- data.frame(lat[which(is.na(lat$MaxDispersalPotentialKmY)==F),])

dsp <- unique(dat.MaxDispPot.lat[,c(which(colnames(dat.MaxDispPot.lat) == "class"),
                                    which(colnames(dat.MaxDispPot.lat) == "Species"))])
cl <- table(dsp$class)
nn <- names(which(cl >= 5))
dat.MaxDispPot.lat <- dat.MaxDispPot.lat[which(dat.MaxDispPot.lat$class %in% nn),]
dat.MaxDispPot.lat <-droplevels(dat.MaxDispPot.lat)
length(unique(dat.MaxDispPot.lat$class)) # 7

#scale variables
# elevation:
dat.MaxDispPot.ele$max_disp_pot <- scale(dat.MaxDispPot.ele$MaxDispersalPotentialmY)
dat.MaxDispPot.ele$shift <- scale(dat.MaxDispPot.ele$ShiftR)
dat.MaxDispPot.ele$lags_scaled <- scale(dat.MaxDispPot.ele$lags)
dat.MaxDispPot.ele$inferred_lag_scaled <- scale(dat.MaxDispPot.ele$inferred_lag)
dat.MaxDispPot.ele$velocity_scaled <- scale(dat.MaxDispPot.ele$EleVeloT)

# latitude:
dat.MaxDispPot.lat$max_disp_pot <- scale(dat.MaxDispPot.lat$MaxDispersalPotentialKmY)
dat.MaxDispPot.lat$shift <- scale(dat.MaxDispPot.lat$ShiftR)
dat.MaxDispPot.lat$lags_scaled <- scale(dat.MaxDispPot.lat$lags)
dat.MaxDispPot.lat$inferred_lag_scaled <- scale(dat.MaxDispPot.lat$inferred_lag)
dat.MaxDispPot.lat$velocity_scaled <- scale(dat.MaxDispPot.lat$LatVeloT)


## split by type of shift 
ele_te <- dat.MaxDispPot.ele[dat.MaxDispPot.ele$Position=="Trailing edge",]
ele_le <- dat.MaxDispPot.ele[dat.MaxDispPot.ele$Position=="Leading edge",]
ele_o <- dat.MaxDispPot.ele[dat.MaxDispPot.ele$Position=="Centroid",]
lat_te <- dat.MaxDispPot.lat[dat.MaxDispPot.lat$Position=="Trailing edge",]
lat_le <- dat.MaxDispPot.lat[dat.MaxDispPot.lat$Position=="Leading edge",]
lat_o <- dat.MaxDispPot.lat[dat.MaxDispPot.lat$Position=="Centroid",]

## model inferred range shift as a function of whether or not we expect the species to be able to track the climate 
mod_lat_le = lm(inferred_lag ~ expect_tracking_pot, 
                  data = lat_le)

summary(mod_lat_le)

lat_le %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, y = inferred_lag, fill = expect_tracking_pot)) + geom_boxplot()

ele_le %>%
  ggplot(aes(x = expect_tracking_pot, y = inferred_lag)) + geom_boxplot()

lat_o %>%
  ggplot(aes(x = expect_tracking_pot, y = inferred_lag)) + geom_boxplot()
ele_o %>%
  ggplot(aes(x = expect_tracking_pot, y = inferred_lag)) + geom_boxplot()


## elevation leading edge
mod_ele_le = lm(inferred_lag ~ expect_tracking_pot, 
                  data = ele_le)

summary(mod_ele_le)

## latitude optimum
mod_lat_o = lm(inferred_lag ~ expect_tracking_pot,
                 data = lat_o)

summary(mod_lat_o)

## elevation optimum
mod_ele_o = lm(inferred_lag ~ expect_tracking_pot, 
                 data = ele_o)

summary(mod_ele_o)




## try fitting one model 
all_data <- rbind(lat_le, ele_le) %>%
  rbind(., lat_o) %>%
  rbind(., ele_o)

model_all <- lm(inferred_lag ~ expect_tracking_pot*Gradient + expect_tracking_pot*Position,
  data = all_data)

summary(model_all)

all_leading = filter(all_data, Position == "Leading edge")

model_all_leading <- lm(inferred_lag ~ expect_tracking_pot*Gradient,
                data = all_leading)

summary(model_all_leading)


ggplot(aes(x = Gradient, y = inferred_lag, colour = expect_tracking_pot), data = all_data) + 
  geom_boxplot() +
  facet_grid(Position ~ class) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  scale_colour_manual(values = c(pal[6], pal[4])) +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       x = "Log maximum dispersal potential (km/y or m/y)", 
       y = "Inferred range shift lag (km/y or m/y)") +
  theme(legend.position = "bottom", panel.grid = element_blank())

all_data %>%
  filter(abs(lags)< 60) %>%
  ggplot(aes(x = Gradient, y = lags, colour = expect_tracking_pot), data = .) + 
  geom_boxplot() +
  facet_grid(Position ~ class) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  scale_colour_manual(values = c(pal[6], pal[4])) +
  labs(colour = "Does dispersal potential allow species to keep up with temp change?", 
       x = "Log maximum dispersal potential (km/y or m/y)", 
       y = "Range shift lag (km/y or m/y)") +
  theme(legend.position = "bottom", panel.grid = element_blank())




## feedback:
## impute lifespan and age at maturity (try in a small subset (birds) and see what makes best imputation)
## Jenn: no matter what you will have to show raw data, so beginning with modelling that isn't a waste of time
## see if v2 has velocity?
## see whether relationships do not exist at trailing edge
## add trailing edge + negative velocity shifts
## ask Sam for dispersal data
## see if we can add mean/mode other dispersal metrics
## try using taxonomic group instead of class (group plants together)
## co linearity between class and max dispersal distance?
## compare species with both elev and lat shifts to see if dispersal potential keeps up with climate velocity better across elev


## modelling options:
## create model that predicts range shift values for each observation you need after accounting for methodological variables
## use these as your dependent variable
## easier to control for methodology when you can use ALL of the points 
## avoids signularities 
## OR
## include the methodological factors in your model 

## look within articles and see if the expected relationship is there 

## jonathon wants to use all dispersal data to look at the evolution of dispersal distance/potential

## Brown methodological effects

## high climate velocity and low range shift = high lag
## see 

## remove log scale 


## non-log versions
pol1 <- data.frame(x = c(0.05, 7, 7, 0.05), y = c(0.00005, 0.00005, 6, 0.05))
pol2 <- data.frame(x = c(0.05, 0.05, 7, 7), y = c(0.05, 150, 150, 6))

min(exp(1)^lat$MaxDispersalPotentialmY)
max(exp(1)^lat$MaxDispersalPotentialmY)
min(lat$LatVeloT)
max(lat$LatVeloT)

lat %>% 
  ggplot(aes(y = exp(1)^MaxDispersalPotentialKmY, x = LatVeloT)) + 
  geom_polygon(data = pol1, aes(x = x, y = y), fill = pal[6], alpha = 0.5) +
  geom_polygon(data = pol2, aes(x = x, y = y), fill = pal[4], alpha = 0.5) +
  geom_point(aes(colour = expect_tracking_pot)) +
  theme_light() + theme(panel.grid = element_blank()) + 
  labs(y = "Maximum dispersal potential (km/y)", 
       x = "Mean velocity of temperature (km/y)", 
       title = "Latitudinal shifts")+
  theme_classic() +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme(legend.position = "none") +
  scale_x_log10(limits = c(0.05, 7), expand = c(0,0)) +
  scale_y_log10(limits = c(0.00005, 150), expand = c(0,0))

ele_exp = filter(ele, expect_tracking_pot == "No")

ele %>% 
  ggplot(aes(y = exp(1)^MaxDispersalPotentialmY, x = EleVeloT)) + 
  geom_point(aes(colour = expect_tracking_pot), size = 0.75) +
  theme_light() + theme(panel.grid = element_blank()) + 
  labs(y = "Maximum dispersal potential (m/y)", 
       x = "Mean velocity of temperature (m/y)", 
       title = "Latitudinal shifts")+
  theme_classic() +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme(legend.position = "none") +
  geom_point(data = ele_exp, aes(colour = expect_tracking_pot), size = 0.75) 

lat_exp = filter(lat, expect_tracking_pot == "No")

lat %>% 
  ggplot(aes(y = exp(1)^MaxDispersalPotentialKmY, x = LatVeloT)) + 
  geom_point(aes(colour = expect_tracking_pot), size = 0.75, position = "jitter") +
  theme_light() + theme(panel.grid = element_blank()) + 
  labs(y = "Maximum dispersal potential (km/y)", 
       x = "Mean velocity of temperature (km/y)", 
       title = "Latitudinal shifts")+
  theme_classic() +
  scale_colour_manual(values = c(pal[6], pal[4])) + 
  theme(legend.position = "none") +
  geom_point(data = lat_exp, aes(colour = expect_tracking_pot), size = 0.75) 
