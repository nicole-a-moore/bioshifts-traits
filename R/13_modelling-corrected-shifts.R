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
library(ggallin)
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
v1$dispersal_potential_kmY = v1$MaxDispersalPotentialKmY
v1$dispersal_potential_mY = v1$MaxDispersalPotentialmY

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


##############################################################
####   test evidence for dispersal x climate interaction  ####
##############################################################
## convert units of latitudinal and elevation shifts & dispersal & climate velocity to km/y:
lags$ShiftKmY <- ifelse(lags$Gradient == "Elevation", lags$ShiftR / 1000, lags$ShiftR)
lags$ClimVeloTKmY <- ifelse(lags$Gradient == "Elevation", lags$EleVeloT / 1000, lags$LatVeloT)
lags$AnnualDispPotKmY <- ifelse(lags$Gradient == "Elevation", lags$annual_dispersal_pot / 1000, lags$annual_dispersal_pot)

## plot corrected shifts 
lags %>%
  ggplot(., aes(x = annual_dispersal_pot, y = CorrShift, colour = Gradient)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10()
## use uncorrected shifts for now since we want to compare elevation and latitude
## might need to correct shifts using a single model (converting units) if we want to plot latitude and elevation together 

## now plot uncorrected shifts 
## colour by gradient 
lags %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = Gradient)) + 
  theme_bw() + 
  geom_point() + 
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "")

## most species have very high dispersal potential compared to the climate velocity and range shift
lags %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ClimVeloTKmY, colour = Gradient)) + 
  theme_bw() + 
  geom_point() + 
  labs(x = "Annual dispersal potential (km/y)", y = "Annual climate velocity (km/y)",
       colour = "") +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits = c(0, 1500))+
  scale_y_continuous(limits = c(0, 1500))

## zoom in: 
lags %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ClimVeloTKmY, colour = Gradient)) + 
  theme_bw() + 
  geom_point() + 
  labs(x = "Annual dispersal potential (km/y)", y = "Annual climate velocity (km/y)",
       colour = "") +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(limits = c(0, 8))+
  scale_y_continuous(limits = c(0, 8))

## species with higher dispersal potential have higher observed shifts 
## line = 1:1 line
lags %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = Gradient)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "") +
  stat_function(colour = "black", fun = function(x){x}) +
  scale_y_continuous(limits = c(0, 40))

## across elevation, no effect of dispersal on shift 
## across latitude, shift increases with dispersal potential 
lags %>%
  filter(ShiftKmY > 0.0001) %>%
  filter(AnnualDispPotKmY < 8) %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = Gradient)) + 
  theme_bw() + 
  geom_point() + 
  # scale_x_log10() +
  # scale_y_log10() +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "") +
  stat_function(colour = "black", fun = function(x){x}) +
  geom_smooth(method = "lm", colour = "black", aes(group = Gradient))
## but lots of species have observed shift > dispersal potential 

## colour by taxonomic group 
lags %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = Group)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "")

## colour by climate velocity 
lags %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() + 
  facet_grid(~Gradient) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "Climate velocity (km/y)")

## just elevation shifts 
lags %>%
  filter(Gradient == "Elevation") %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() + 
  facet_grid(~Gradient) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "Climate velocity (km/y)") +
  geom_smooth(method = "lm", se = FALSE, colour = "black") 

## just latitudinal shifts 
lags %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() + 
  facet_grid(~Gradient) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "Climate velocity (km/y)") +
  geom_smooth(method = "lm", se = FALSE, colour = "black")


######################
####   modelling  ####
######################
## bin climate velocity 
ele <- filter(lags, Gradient == "Elevation") %>%
  mutate(Group = factor(Group, ordered = F))

hist(ele$ClimVeloTKmY)
q = quantile(ele$ClimVeloTKmY, probs = c(0,0.2,0.4,0.6,0.8,1))
ele$ClimVeloTKmY = cut(ele$ClimVeloTKmY,
                       breaks = q,
                       include.lowest = T)
plot(ele$ClimVeloTKmY)

ele$ClimVeloTKmY <- str_replace_all(ele$ClimVeloTKmY, "\\[", "(") 
ele$ClimVeloTKmY <- str_replace_all(ele$ClimVeloTKmY, "\\]", ")") 

ele <- ele %>%
  mutate(quant_max = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,2], "\\)", " "))) %>%
  mutate(quant_min = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,1], "\\(", " "))) %>%
  mutate(quant_mean = (quant_min + quant_max)/2) 

## model everything beyond where dispersal = climate velocity 
## model log-log relationship 
## log shift ~ log annual dispersal potential + climate velocity 
## expect flat relationship, and expect higher CV to have higher intercept 
lat <- filter(lags, Gradient == "Latitudinal") %>%
  mutate(Group = factor(Group, ordered = F))

hist(lat$ClimVeloTKmY)
q = quantile(lat$ClimVeloTKmY, probs = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
lat$ClimVeloTKmY = cut(lat$ClimVeloTKmY,
                       breaks = q,
                       include.lowest = T)
plot(lat$ClimVeloTKmY)

lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\[", "(") 
lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\]", ")") 

mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)

lat <- lat %>%
  mutate(quant_max = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,2], "\\)", " "))) %>%
  mutate(quant_min = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,1], "\\(", " "))) %>%
  mutate(quant_mean = (quant_min + quant_max)/2) %>%
  mutate(ClimVeloTKmY = paste(quant_min, " km/y - ", quant_max," km/y" , sep = ""))

lat %>%
  filter(ShiftR > 0.0001) %>%
  filter(annual_dispersal_pot < 8) %>%
  ## add lines to indicate where inflection point is expected (where mean cv = annual dp = observed shift)
  ggplot(., aes(x = annual_dispersal_pot, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  #scale_x_log10() +
 # scale_y_log10() +
  # scale_y_continuous(trans = pseudolog10_trans) +
  facet_wrap(~ClimVeloTKmY, nrow = 2) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "") +
  scale_colour_manual(values = mycolours) +
  geom_vline(aes(xintercept = quant_min)) +
  geom_hline(aes(yintercept = quant_min)) +
  #geom_abline(intercept = 0, slope = 1, colour = "black") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

ele %>%
  filter(ShiftR > 0.0001) %>%
  filter(annual_dispersal_pot < 5) %>%
  ## add lines to indicate where inflection point is expected (where mean cv = annual dp = observed shift)
  ggplot(., aes(x = annual_dispersal_pot, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  #scale_x_log10() +
  # scale_y_log10() +
  # scale_y_continuous(trans = pseudolog10_trans) +
  facet_wrap(~ClimVeloTKmY, nrow = 2) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "") +
  scale_colour_manual(values = mycolours) +
  geom_vline(aes(xintercept = quant_min)) +
  geom_hline(aes(yintercept = quant_min)) +
  #geom_abline(intercept = 0, slope = 1, colour = "black") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(path = "figures/sotm", filename = "all-data.png", 
       device = "png", height = 4, width = 8.5)


model_data <- lat %>%
  filter(ShiftR > 0.0001) %>%
  filter(annual_dispersal_pot >= quant_min) %>%
  arrange(quant_min) %>%
  ungroup() %>%
  mutate(ClimVeloTKmY = factor(ClimVeloTKmY, 
                               levels = unique(ClimVeloTKmY),
                               ordered = TRUE))


## make log version of dispersal potential
model_data$log_annual_dispersal_pot = log(model_data$annual_dispersal_pot)
model_data$log_shift = log(model_data$ShiftKmY)

mod <- lm(log_shift ~ log_annual_dispersal_pot + ClimVeloTKmY,
          data = model_data)
summary(mod)
#plot(mod)

pred_data <- expand.grid(log_annual_dispersal_pot = seq(from = min(model_data$log_annual_dispersal_pot),
                                                        to = max(model_data$log_annual_dispersal_pot), 
                                                        by = 0.01),
                         ClimVeloTKmY = unique(model_data$ClimVeloTKmY))

pred <- predict(mod, pred_data, se.fit = TRUE)
pred_data$predicted_val = pred$fit
pred_data$predicted_se = pred$se.fit


## restrict X values 
bounds <- select(model_data, ClimVeloTKmY, log_annual_dispersal_pot) %>%
  group_by(ClimVeloTKmY) %>%
  mutate(min_log_annual_dispersal_pot = min(log_annual_dispersal_pot),
         max_log_annual_dispersal_pot = max(log_annual_dispersal_pot)) %>%
  select(-log_annual_dispersal_pot) %>%
  distinct()

pred_data = left_join(pred_data, bounds) %>%
  filter(log_annual_dispersal_pot >= min_log_annual_dispersal_pot & 
           log_annual_dispersal_pot <= max_log_annual_dispersal_pot)
  
model_data %>%
  ## add lines to indicate where inflection point is expected (where mean cv = annual dp = observed shift)
  ggplot(., aes(x = annual_dispersal_pot, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() +
  scale_y_log10() +
  # scale_y_continuous(trans = pseudolog10_trans) +
  facet_wrap(~ClimVeloTKmY, nrow = 2) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "") +
  scale_colour_manual(values = mycolours) +
  geom_vline(aes(xintercept = quant_min)) +
  geom_hline(aes(yintercept = quant_min)) +
  #geom_abline(intercept = 0, slope = 1, colour = "black") + 
  geom_smooth(data = pred_data, 
         aes(x = exp(1)^log_annual_dispersal_pot, 
             y = exp(1)^predicted_val,
             col = ClimVeloTKmY),
         colour = "black", linewidth = 1.25) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
# +
#   geom_ribbon(data = pred_data, 
#               aes(x = exp(1)^log_annual_dispersal_pot, 
#                   y = exp(1)^predicted_val,
#                   ymax = exp(1)^predicted_val+exp(1)^predicted_se,
#                   ymin = exp(1)^predicted_val-exp(1)^predicted_se),
#               colour = "#ffffff00",
#               alpha = 0.3)
  
ggsave(path = "figures/sotm", filename = "model-predictions_not-dispersal-limited.png", 
       device = "png", height = 4, width = 8.5)

## plot model predictions 
ggplot(pred_data, aes(x = ClimVeloTKmY, y = exp(1)^predicted_val, colour = ClimVeloTKmY)) +
  geom_point() +
  scale_y_log10() +
  scale_colour_manual(values = mycolours) +
  theme_bw() +
  labs(x = "Annual climate velocity (km/y)", y = "Predicted range shift (km/y)")


## plot predicted value at inflection point against predicted inflection point 
pred_data <- model_data %>%
  mutate(log_annual_dispersal_pot = log(quant_min)) %>%
  select(ClimVeloTKmY, log_annual_dispersal_pot) %>%
  distinct()

pred <- predict(mod, pred_data, se.fit = TRUE)
pred_data$predicted_val = pred$fit
pred_data$predicted_se = pred$se.fit

ggplot(pred_data, aes(x = ClimVeloTKmY, y = exp(1)^predicted_val, colour = ClimVeloTKmY))  +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = exp(1)^predicted_val - exp(1)^predicted_se,
                    ymax = exp(1)^predicted_val + exp(1)^predicted_se)) +
  geom_point(data = model_data, aes(y = quant_min), colour = "black", size = 3) +
  scale_colour_manual(values = mycolours) +
  theme_bw() +
  scale_y_continuous(trans = pseudolog10_trans, limits = c(-1.1, 3.1)) +
  labs(x = "Annual climate velocity (km/y)", y = "Predicted range shift (km/y)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 1.5), "cm")) 

ggsave(path = "figures/sotm", filename = "model-predictions_not-dispersal-limited_intercept.png", 
       device = "png", height = 4, width = 7)






## now model relationship before inflection point 
## expect 1:1 relationship, expect intercept to be the same across climate velocities
model_data <- lat %>%
  filter(annual_dispersal_pot < quant_min) %>%
  arrange(quant_min) %>%
  ungroup() %>%
  mutate(ClimVeloTKmY = factor(ClimVeloTKmY, 
                               levels = unique(ClimVeloTKmY),
                               ordered = TRUE)) 
  


## make log version of dispersal potential
model_data$log_annual_dispersal_pot = log(model_data$annual_dispersal_pot)
model_data$log_shift = log(model_data$ShiftKmY)

mod <- lm(log_shift ~ log_annual_dispersal_pot + ClimVeloTKmY,
          data = model_data)
summary(mod)
#plot(mod)

pred_data <- expand.grid(log_annual_dispersal_pot = seq(from = min(model_data$log_annual_dispersal_pot),
                                                        to = max(model_data$log_annual_dispersal_pot), 
                                                        by = 0.01),
                         ClimVeloTKmY = unique(model_data$ClimVeloTKmY))

pred <- predict(mod, pred_data, se.fit = TRUE)
pred_data$predicted_val = pred$fit
pred_data$predicted_se = pred$se.fit

## restrict X values 
bounds <- select(model_data, ClimVeloTKmY, log_annual_dispersal_pot) %>%
  group_by(ClimVeloTKmY) %>%
  mutate(min_log_annual_dispersal_pot = min(log_annual_dispersal_pot),
         max_log_annual_dispersal_pot = max(log_annual_dispersal_pot)) %>%
  select(-log_annual_dispersal_pot) %>%
  distinct()

pred_data = left_join(pred_data, bounds) %>%
  filter(log_annual_dispersal_pot >= min_log_annual_dispersal_pot & 
           log_annual_dispersal_pot <= max_log_annual_dispersal_pot)

pred_data %>%
  ggplot(aes(x = exp(1)^log_annual_dispersal_pot, 
             y = exp(1)^predicted_val,
             col = ClimVeloTKmY)) +
  geom_point() +
  scale_colour_manual(values = mycolours) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10()

lat %>%
  filter(ShiftR > 0.0001) %>%
  #filter(annual_dispersal_pot >= ShiftR) %>%
  ## add lines to indicate where inflection point is expected (where mean cv = annual dp = observed shift)
  filter(annual_dispersal_pot < quant_min) %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  facet_wrap(~ClimVeloTKmY, nrow = 2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "") +
  scale_colour_manual(values = mycolours) +
  geom_vline(aes(xintercept = quant_min)) +
  geom_hline(aes(yintercept = quant_min)) +
  geom_abline(intercept = 0, slope = 1, colour = "black") +
  geom_smooth(data = pred_data, 
              aes(x = exp(1)^log_annual_dispersal_pot, 
                  y = exp(1)^predicted_val, 
                  group = ClimVeloTKmY),
              colour = "black") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(path = "figures/sotm", filename = "model-predictions_dispersal-limited.png", 
       device = "png", height = 4, width = 8.5)
## in most cases where dispersal potential was less than climate velocity, observed shift was greater than dispersal potential 

## grey out very erroneous ones 
alpha = filter(lat, ShiftKmY >= AnnualDispPotKmY + 0.05) 
lat %>%
  filter(ShiftR > 0.0001) %>%
  #filter(annual_dispersal_pot >= ShiftR) %>%
  ## add lines to indicate where inflection point is expected (where mean cv = annual dp = observed shift)
  filter(annual_dispersal_pot < quant_min) %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  geom_point(data = alpha, colour = "grey") + 
  #facet_wrap(~ClimVeloTKmY, nrow = 2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "") +
  scale_colour_manual(values = mycolours) +
  geom_vline(aes(xintercept = quant_min)) +
  geom_hline(aes(yintercept = quant_min)) +
  geom_abline(intercept = 0, slope = 1, colour = "black") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))


## now model all together to see if shifts are greater in spp with larger dispersal potential 
lat %>%
  filter(ShiftR > 0.0001) %>%
  mutate(expect_limited = ifelse(expect_tracking == "Yes", "No", "Yes")) %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = expect_limited)) + 
  theme_bw() + 
  geom_point() + 
  facet_wrap(~ClimVeloTKmY, nrow = 2) +
  scale_x_log10() +
  scale_y_log10(limits = c(0.00001, 40)) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "Expect dispersal limitation?") +
  geom_boxplot() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(path = "figures/sotm", filename = "latitude.png", 
       device = "png", height = 4, width = 8.5)

ele %>%
  filter(ShiftKmY > 0.0001) %>%
  mutate(expect_limited = ifelse(expect_tracking == "Yes", "No", "Yes")) %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = expect_limited)) + 
  theme_bw() + 
  geom_point() + 
  facet_wrap(~ClimVeloTKmY, nrow = 2) +
  scale_x_log10() +
  scale_y_log10(limits = c(0.00001, 40)) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "Expect dispersal limitation?") +
  geom_boxplot() +x
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(path = "figures/sotm", filename = "elevation.png", 
       device = "png", height = 4, width = 6.25)


lat %>%
  filter(ShiftR > 0.0001) %>%
  mutate(expect_limited = ifelse(expect_tracking == "Yes", "No", "Yes")) %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = expect_limited)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() +
  scale_y_log10(limits = c(0.00001, 40)) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "Expect dispersal limitation?") +
  #geom_boxplot() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(path = "figures/sotm", filename = "latitude_yes-no.png", 
       device = "png", height = 4, width = 4.5)

ele %>%
  filter(ShiftR > 0.0001) %>%
  mutate(expect_limited = ifelse(expect_tracking == "Yes", "No", "Yes")) %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = expect_limited)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() +
  scale_y_log10(limits = c(0.00001, 40)) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "Expect dispersal limitation?") +
  #geom_boxplot() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(path = "figures/sotm", filename = "elevation_yes-no.png", 
       device = "png", height = 4, width = 4.5)


## for species whose dispersal potential < obsevred shift, how big is the difference?
lags %>%
  filter(ShiftR > 0.0001) %>%
  filter(ShiftKmY > AnnualDispPotKmY) %>%
  mutate(diff = ShiftKmY - AnnualDispPotKmY) %>%
  mutate(mean_diff = mean(.$diff)) %>%
  ggplot(., aes(x = diff)) + 
  theme_bw() + 
  geom_histogram() +
  labs(x = "Difference between observed expansion rate and dispersal potential",
       y = "Frequency")

lags %>%
  filter(ShiftR > 0.0001) %>%
  filter(ShiftKmY > AnnualDispPotKmY) %>%
  mutate(diff = ShiftKmY - AnnualDispPotKmY) %>%
  mutate(mean_diff = mean(.$diff)) %>%
  ggplot(., aes(x = diff)) + 
  theme_bw() + 
  geom_histogram() +
  scale_x_log10() +
  labs(x = "Difference between observed expansion rate and dispersal potential",
       y = "Frequency")

## how many?
lags %>%
  mutate(underestimate = ifelse(ShiftKmY > AnnualDispPotKmY, "Yes", "No")) %>%
  group_by(underestimate) %>%
  tally()

## how many can still disperse faster than temperature change?
lags %>%
  mutate(expect_limited = ifelse(expect_tracking == "Yes", "No", "Yes")) %>%
  group_by(expect_limited) %>%
  tally()

lags %>%
  filter(ShiftR > 0.0001) %>%
  filter(ShiftKmY <= AnnualDispPotKmY + 0.001) %>%
  mutate(expect_limited = ifelse(expect_tracking == "Yes", "No", "Yes")) %>%
  ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() +
  scale_y_log10(limits = c(0.00001, 40)) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "Expect dispersal limitation?") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1)) +
  geom_abline(intercept = 0, slope = 1) 


  



## modelling concerns:
## shift is 0 bounded 
hist(model_data$ShiftKmY)
## modeling error structure might be important here - some studies have methodology that causes more measurement error, climate velocity bins are highly related to study
## also expect more signal:noise where climate velocity is high



model_data <- lat %>%
  filter(ShiftR > 0.0001) %>%
  arrange(quant_min) %>%
  ungroup() %>%
  mutate(ClimVeloTKmY = factor(ClimVeloTKmY, 
                               levels = unique(ClimVeloTKmY),
                               ordered = TRUE))


## make log version of dispersal potential
model_data$log_annual_dispersal_pot = log(model_data$annual_dispersal_pot)
model_data$log_shift = log(model_data$ShiftKmY)

mod <- lm(log_shift ~ log_annual_dispersal_pot + ClimVeloTKmY,
          data = model_data)
summary(mod)
#plot(mod)

pred_data <- expand.grid(log_annual_dispersal_pot = seq(from = min(model_data$log_annual_dispersal_pot),
                                                        to = max(model_data$log_annual_dispersal_pot), 
                                                        by = 0.01),
                         ClimVeloTKmY = unique(model_data$ClimVeloTKmY))

pred <- predict(mod, pred_data, se.fit = TRUE)
pred_data$predicted_val = pred$fit
pred_data$predicted_se = pred$se.fit


## restrict X values 
bounds <- select(model_data, ClimVeloTKmY, log_annual_dispersal_pot) %>%
  group_by(ClimVeloTKmY) %>%
  mutate(min_log_annual_dispersal_pot = min(log_annual_dispersal_pot),
         max_log_annual_dispersal_pot = max(log_annual_dispersal_pot)) %>%
  select(-log_annual_dispersal_pot) %>%
  distinct()

pred_data = left_join(pred_data, bounds) %>%
  filter(log_annual_dispersal_pot >= min_log_annual_dispersal_pot & 
           log_annual_dispersal_pot <= max_log_annual_dispersal_pot)

model_data %>%
  ## add lines to indicate where inflection point is expected (where mean cv = annual dp = observed shift)
  ggplot(., aes(x = annual_dispersal_pot, y = ShiftKmY, colour = ClimVeloTKmY)) + 
  theme_bw() + 
  geom_point() + 
  scale_x_log10() +
  scale_y_log10() +
  # scale_y_continuous(trans = pseudolog10_trans) +
  facet_wrap(~ClimVeloTKmY, nrow = 2) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "") +
  scale_colour_manual(values = mycolours) +
  geom_vline(aes(xintercept = quant_min)) +
  geom_hline(aes(yintercept = quant_min)) +
  #geom_abline(intercept = 0, slope = 1, colour = "black") + 
  geom_smooth(data = pred_data, 
              aes(x = exp(1)^log_annual_dispersal_pot, 
                  y = exp(1)^predicted_val,
                  col = ClimVeloTKmY),
              colour = "black", linewidth = 1.25) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
# +
#   geom_ribbon(data = pred_data, 
#               aes(x = exp(1)^log_annual_dispersal_pot, 
#                   y = exp(1)^predicted_val,
#                   ymax = exp(1)^predicted_val+exp(1)^predicted_se,
#                   ymin = exp(1)^predicted_val-exp(1)^predicted_se),
#               colour = "#ffffff00",
#               alpha = 0.3)

ggsave(path = "figures/sotm", filename = "model-predictions_not-dispersal-limited_full-dataset.png", 
       device = "png", height = 4, width = 8.5)

## plot model predictions 
ggplot(pred_data, aes(x = ClimVeloTKmY, y = exp(1)^predicted_val, colour = ClimVeloTKmY)) +
  geom_point() +
  scale_y_log10() +
  scale_colour_manual(values = mycolours) +
  theme_bw() +
  labs(x = "Annual climate velocity (km/y)", y = "Predicted range shift (km/y)")


## plot predicted value at inflection point against predicted inflection point 
pred_data <- model_data %>%
  mutate(log_annual_dispersal_pot = log(quant_min)) %>%
  select(ClimVeloTKmY, log_annual_dispersal_pot) %>%
  distinct()

pred <- predict(mod, pred_data, se.fit = TRUE)
pred_data$predicted_val = pred$fit
pred_data$predicted_se = pred$se.fit

ggplot(pred_data, aes(x = ClimVeloTKmY, y = exp(1)^predicted_val, colour = ClimVeloTKmY))  +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = exp(1)^predicted_val - exp(1)^predicted_se,
                    ymax = exp(1)^predicted_val + exp(1)^predicted_se)) +
  geom_point(data = model_data, aes(y = quant_min), colour = "black", size = 3) +
  scale_colour_manual(values = mycolours) +
  theme_bw() +
  scale_y_continuous(trans = pseudolog10_trans, 
                     limits = c(-1.1, 3.1)) +
  labs(x = "Annual climate velocity (km/y)", y = "Predicted range shift (km/y)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 1.5), "cm")) 

ggsave(path = "figures/sotm", filename = "model-predictions_not-dispersal-limited_intercept_full-dataset.png", 
       device = "png", height = 4, width = 7)


