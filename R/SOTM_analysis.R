## SOTM figures
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

## convert units of latitudinal and elevation shifts & dispersal & climate velocity to km/y:
lags$ShiftKmY <- ifelse(lags$Gradient == "Elevation", lags$ShiftR / 1000, lags$ShiftR)
lags$ClimVeloTKmY <- ifelse(lags$Gradient == "Elevation", lags$EleVeloT / 1000, lags$LatVeloT)
lags$AnnualDispPotKmY <- ifelse(lags$Gradient == "Elevation", lags$annual_dispersal_pot / 1000, lags$annual_dispersal_pot)

mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)

lags <- select(lags, -c(lag, SLDiff, CorrShift, PredSLShift, study_level_shift)) %>%
  distinct()

## make plots to summarize subset 
nrow(lags) # 1935 range shifts
length(unique(lags$Species)) # 409 species 
  
grad <- lags %>%
  ggplot(aes(x = Gradient, fill = Gradient)) +
  geom_bar() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs (x = "", y = "Number of\nrange expansion\nobservations") +
  scale_x_discrete(labels = c("Elevation studies", "Latitudinal studies")) +
  scale_fill_manual(values = c("#28587B", "#9FB798")) 

ggsave(grad, path = "figures/sotm", filename = "barplot-gradients.png", 
       device = "png", height = 2, width = 4)

lags %>%
  group_by(Gradient) %>% tally()
  
groups <- lags %>%
  ggplot(aes(x = Group)) +
  geom_bar() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs (x = "", y = "Number of\nrange expansion\nobservations") 

ggsave(groups, path = "figures/sotm", filename = "barplot-groups.png", 
       device = "png", height = 2, width = 4)

lags %>%
  group_by(Group) %>% tally()

#################################
##    make expectation plot    ##
#################################
hist(lags$ClimVeloTKmY)

## rug data
rug <- lags %>%
  mutate(x = 0, y = 0)
  
## climate velocity 
rug %>%
  ggplot(aes(x = x, y = ClimVeloTKmY, colour = Gradient)) +
  geom_point(shape = "_", size = 8) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 7), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 7), expand = c(0,0)) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)")

## make histograms
lags %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = ClimVeloTKmY, fill = ..x..)) +
  geom_histogram() + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-1, 15), expand = c(0,0)) +
  labs(x = "Mean climate velocity across study area (km/y)",
       y = "Number of range expansion observations") +
  facet_grid(~Gradient) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"),
        legend.position = "none") +
  scale_fill_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 4.5) 

ggsave(path = "figures/sotm", filename = "hists.png", 
       device = "png", height = 4, width = 8)

## make new dataframe for expectation plot
clim_velos <- c(lags$ClimVeloTKmY)

y_values = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 15, by = 0.01)) 
  return(rep(x, reps))
}))
reps = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 15, by = 0.01)) 
  return(reps)
}))
x_values = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = x, 15, by = 0.01))
}))

y_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.01))
}))
x_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.01))
}))
velos <- unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = 0, to = x, by = 0.01)) 
  return(rep(x, reps))
}))

clim_data <- data.frame(y_values = append(y_values_11, y_values), 
                        x_values = append(x_values_11, x_values), 
                        velos = append(velos, y_values))

gradients <- select(ungroup(lags), ClimVeloTKmY, Gradient) %>%
  distinct()

clim_data <- left_join(clim_data, gradients, by = c("velos" = "ClimVeloTKmY"))
  
one <- clim_data %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 15), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 15), expand = c(0,0),
                     labels = c("   0", "   5", " 10", "   15")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "Mean climate\nvelocity\nacross study\narea (km/y)") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 

one_legend <- get_legend(one)

one <- one +
  theme(legend.position = "none")

ggsave(one, path = "figures/sotm", filename = "expectation-climate-1.png", 
       device = "png", height = 4, width = 5)
ggsave(one_legend, path = "figures/sotm", filename = "expectation-climate-1-legend.png", 
       device = "png", height = 2, width = 1)

## split by latitude versus elevation
two = clim_data %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 15), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 15), expand = c(0,0),
                     labels = c("   0", "   5", " 10", "   15")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "Mean climate\nvelocity\nacross study\narea (km/y)") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  facet_grid(~Gradient) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"))

two_legend = get_legend(two)

two <- two +
  theme(legend.position = "none")

ggsave(two, path = "figures/sotm", filename = "expectation-climate-2.png", 
       device = "png", height = 4, width = 8)
ggsave(two_legend, path = "figures/sotm", filename = "expectation-climate-2-legend.png", 
       device = "png", height = 2.4, width = 1.5)

## zoom in
filtered_lags <- filter(lags, ClimVeloTKmY < 0.5)
clim_velos <- filtered_lags$ClimVeloTKmY

y_values = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 1, by = 0.0001)) 
  return(rep(x, reps))
}))
reps = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 1, by = 0.0001)) 
  return(reps)
}))
x_values = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = x, 1, by = 0.0001))
}))

y_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.0001))
}))
x_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.0001))
}))
velos <- unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = 0, to = x, by = 0.0001)) 
  return(rep(x, reps))
}))

clim_data_zoom <- data.frame(y_values = append(y_values_11, y_values), 
                        x_values = append(x_values_11, x_values), 
                        velos = append(velos, y_values))

gradients <- select(ungroup(filtered_lags), ClimVeloTKmY, Gradient) %>%
  distinct()

clim_data_zoom <- left_join(clim_data_zoom, gradients, by = c("velos" = "ClimVeloTKmY"))

three = clim_data_zoom %>%
  filter(x_values <= 0.05) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 0.5, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  scale_x_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  labs(x = "",
       y = "", 
       colour = "") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 

three_legend = get_legend(three)

three <- three +
  theme(legend.position = "none")

ggsave(three, path = "figures/sotm", filename = "expectation-climate-3.png", 
       device = "png",  height = 2.8, width = 3.1)
ggsave(three, path = "figures/sotm", filename = "expectation-climate-3-large.png", 
       device = "png", height = 4, width = 4.2)
ggsave(three_legend, path = "figures/sotm", filename = "expectation-climate-3-legend.png", 
       device = "png", height = 1, width = 1.5)


## now split by latitude and elevation
ele_split <- clim_data_zoom %>%
  filter(x_values <= 0.05) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none")

ggsave(ele_split, path = "figures/sotm", filename = "expectation-elevation.png", 
       device = "png", height = 4, width = 4.2)

lat_split = clim_data %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 15), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 15), expand = c(0,0),
                     labels = c("    0", "   5", " 10", "   15")) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  facet_grid(~Gradient) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"))
  
ggsave(lat_split, path = "figures/sotm", filename = "expectation-latitude.png", 
       device = "png", height = 4, width = 8)


## add dispersal potential rug plot 
rug_lat <- filter(rug, Gradient == "Latitudinal") %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE))
rug_ele <- filter(rug, Gradient == "Elevation") %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE))

rugplot_lat <- clim_data %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  facet_grid(~Gradient) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 15), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 15), expand = c(0,0),
                     labels = c("    0", "   5", " 10", "   15")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  geom_rug(data = rug_lat,
             aes(x = AnnualDispPotKmY),
             inherit.aes = FALSE) +
  geom_rug(data = rug_ele,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"))

ggsave(rugplot_lat, path = "figures/sotm", filename = "expectation-latitude-rug-zoom.png", 
       device = "png", height = 4, width = 8)

rugplot_ele_zoom_squish <- clim_data%>%
  #filter(x_values <= 0.05) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 1400), expand = expansion(mult = c(0, 0), 
                                                             add = c(20, 0))) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_colour_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_rug(data = rug_ele, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) 

ggsave(rugplot_ele_zoom_squish, path = "figures/sotm", filename = "expectation-elevation-rug.png", 
       device = "png",  height = 2.8, width = 3.1)

rugplot_ele <- clim_data_zoom %>%
  filter(x_values <= 0.05) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 0.5, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  geom_rug(data = rug_ele, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) 

ggsave(rugplot_ele, path = "figures/sotm", filename = "expectation-elevation-rug-zoom.png", 
       device = "png",  height = 2.8, width = 3.1)

rugplot_lat_inset <- clim_data %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 15), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,15), expand = c(0,0))+
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_rug(data = rug_lat, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) 

ggsave(rugplot_lat_inset, path = "figures/sotm", filename = "expectation-latitude-rug-zoom.png", 
       device = "png",  height = 2.8, width = 3)

## zoom out
## make new dataframe for expectation plot
clim_velos <- c(lags$ClimVeloTKmY)

y_values = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 1400, by = 1)) 
  return(rep(x, reps))
}))
reps = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 15, by = 0.01)) 
  return(reps)
}))
x_values = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = x, 1400, by = 1))
}))

y_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.01))
}))
x_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.01))
}))
velos <- unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = 0, to = x, by = 0.01)) 
  return(rep(x, reps))
}))

clim_data <- data.frame(y_values = append(y_values_11, y_values), 
                        x_values = append(x_values_11, x_values), 
                        velos = append(velos, y_values))

gradients <- select(ungroup(lags), ClimVeloTKmY, Gradient) %>%
  distinct()

clim_data <- left_join(clim_data, gradients, by = c("velos" = "ClimVeloTKmY"))


rugplot_lat_unzoom <- clim_data %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  facet_grid(~Gradient) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 1400), expand = expansion(mult = c(0, 0), 
                                                    add = c(20, 0))) +
  scale_y_continuous(limits = c(0, 15), expand = c(0,0),
                     labels = c("    0", "   5", " 10", "   15")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  geom_rug(data = rug_lat,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  geom_rug(data = rug_ele,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"))

ggsave(rugplot_lat_unzoom, path = "figures/sotm", filename = "expectation-latitude-rug-unzoom.png", 
       device = "png", height = 4, width = 8)

rugplot_ele_unzoom <- clim_data %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_path(linewidth = 0.5, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(limits = c(0, 1400), expand = expansion(mult = c(0, 0), 
                                                    add = c(20, 0))) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none") +
  geom_rug(data = rug_ele, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) 

ggsave(rugplot_ele_unzoom, path = "figures/sotm", filename = "expectation-elevation-rug.png", 
       device = "png", height = 4, width = 4.2)

clim_data_test <- clim_data %>%
  arrange(velos, x_values, y_values)

## make log version
rugplot_lat_unzoom_log <- clim_data %>%
  arrange(x_values, y_values) %>% 
  distinct() %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  facet_grid(~Gradient) +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  geom_path(linewidth = 1, alpha = 1, linejoin = "mitre") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(limits = c(0.0001, 1400), 
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_y_continuous(limits = c(0, 15), expand = c(0,0),
                     labels = c("    0", "   5", " 10", "   15")) +
  theme(legend.position = "none") +
  geom_rug(data = rug_lat,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  geom_rug(data = rug_ele,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "") 

ggsave(rugplot_lat_unzoom_log, path = "figures/sotm", filename = "expectation-latitude-rug-unzoom-log.png", 
       device = "png", height = 4, width = 8)

rugplot_lat_unzoom_log_line <- clim_data %>%
  arrange(x_values, y_values) %>% 
  distinct() %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = x_values, y = y_values, colour = y_values, group = velos)) +
  facet_grid(~Gradient) +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  geom_path(linewidth = 1, alpha = 1, linejoin = "mitre") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(limits = c(0.0001, 1400), 
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_y_continuous(limits = c(0, 15), expand = c(0,0),
                     labels = c("    0", "   5", " 10", "   15")) +
  theme(legend.position = "none") +
  geom_rug(data = rug_lat,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  geom_rug(data = rug_ele,
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Expected range expansion rate (km/y)", 
       colour = "") +
  geom_vline(aes(xintercept = max(velos)))

ggsave(rugplot_lat_unzoom_log_line, path = "figures/sotm", filename = "expectation-latitude-rug-unzoom-log-line.png", 
       device = "png", height = 4, width = 8)


clim_velos <- c(lags$ClimVeloTKmY)

y_values = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 1400, by = 1)) 
  return(rep(x, reps))
}))
reps = unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = x, 15, by = 0.0001)) 
  return(reps)
}))
x_values = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = x, 1400, by = 1))
}))

y_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.0001))
}))
x_values_11 = unlist(sapply(clim_velos, FUN = function(x) {
  return(seq(from = 0, to = x, by = 0.0001))
}))
velos <- unlist(sapply(clim_velos, FUN = function(x) {
  reps = length(seq(from = 0, to = x, by = 0.0001)) 
  return(rep(x, reps))
}))

clim_data <- data.frame(y_values = append(y_values_11, y_values), 
                        x_values = append(x_values_11, x_values), 
                        velos = append(velos, y_values))

gradients <- select(ungroup(lags), ClimVeloTKmY, Gradient) %>%
  distinct()

clim_data <- left_join(clim_data, gradients, by = c("velos" = "ClimVeloTKmY"))


rugplot_ele_unzoom_log <- clim_data %>%
  filter(Gradient == "Elevation") %>%
  arrange(x_values, y_values) %>% 
  distinct() %>%
  ggplot(aes(x = x_values, y = y_values, colour = velos, group = velos)) +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(limits = c(0.0001, 1400), expand = c(0,0),
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_rug(data = rug_ele, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) 

ggsave(rugplot_ele_unzoom_log, path = "figures/sotm", filename = "expectation-elevation-rug-log.png", 
       device = "png", height = 2.8, width = 3.1)

rugplot_ele_unzoom_log_line <- clim_data %>%
  filter(Gradient == "Elevation") %>%
  arrange(x_values, y_values) %>% 
  distinct() %>%
  ggplot(aes(x = x_values, y = y_values, colour = velos, group = velos)) +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  geom_path(linewidth = 1, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_log10(limits = c(0.0001, 1400), expand = c(0,0),
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  scale_y_continuous(limits = c(0, 0.02), expand = c(0,0)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_rug(data = rug_ele, 
           aes(x = AnnualDispPotKmY),
           inherit.aes = FALSE) +
  geom_vline(aes(xintercept = max(velos)))

ggsave(rugplot_ele_unzoom_log_line, path = "figures/sotm", filename = "expectation-elevation-rug-log-line.png", 
       device = "png", height = 4, width = 4.2)

## compare distributions directly 
dens_grad <- lags %>%
  gather(key = "Measure", value = "Measurement", c(AnnualDispPotKmY, ClimVeloTKmY)) %>%
  ggplot(aes(x = Measurement, fill = Measure)) +
  geom_density(alpha = 0.1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~Gradient) +
  scale_x_log10() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.8)) +
  scale_fill_discrete(labels = c("Annual dispersal potential",
                               "Climate velocity")) +
  labs(x = "", y = "Density", fill = "")
  

ggsave(dens_grad, path = "figures/sotm", filename = "density_gradient.png", 
       device = "png", height = 3.5, width = 8.5)

dens_taxa <- lags %>%
  gather(key = "Measure", value = "Measurement", c(AnnualDispPotKmY, ClimVeloTKmY)) %>%
  ggplot(aes(x = Measurement, fill = Measure)) +
  geom_density(alpha = 0.1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(Group~Gradient) +
  scale_x_log10() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2)) +
  scale_fill_discrete(labels = c("Annual dispersal potential",
                                 "Climate velocity")) +
  labs(x = "", y = "Density", fill = "")

ggsave(dens_taxa, path = "figures/sotm", filename = "density_taxa.png", 
       device = "png", height = 3.5, width = 8.5)

cv_vs_dp <- lags %>%
  ggplot(aes(y = ClimVeloTKmY, x = AnnualDispPotKmY, colour = Gradient)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_colour_manual(labels = c("Elevation study", 
                                 "Latitudinal study"),
                      values = c("#28587B", "#9FB798")) +
  labs(x = "Dispersal potential (km/y)", y = "Climate velocity (km/y)") +
  theme(legend.position = "none")

ggsave(cv_vs_dp, path = "figures/sotm", filename = "cv_vs_dp.png", 
       device = "png", height = 3.5, width = 4.5)

cv_vs_dp_zoom <- lags %>%
  ggplot(aes(y = ClimVeloTKmY, x = AnnualDispPotKmY, colour = Gradient)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +  
  scale_x_continuous(limits = c(0, 7)) +  
  scale_y_continuous(limits = c(0, 7)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_colour_manual(labels = c("Elevation study", 
                                 "Latitudinal study"),
                      values = c("#28587B", "#9FB798")) +
  labs(x = "Dispersal potential (km/y)", y = "Climate velocity (km/y)") +
  theme(legend.position = "none")

ggsave(cv_vs_dp_zoom, path = "figures/sotm", filename = "cv_vs_dp_zoom.png", 
       device = "png", height = 3.5, width = 4.5)

cv_vs_dp_zoomzoom <- lags %>%
  ggplot(aes(y = ClimVeloTKmY, x = AnnualDispPotKmY, colour = Gradient)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +  
  scale_x_continuous(limits = c(0, 0.015)) +  
  scale_y_continuous(limits = c(0, 0.015)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_colour_manual(labels = c("Elevation study", 
                                 "Latitudinal study"),
                      values = c("#28587B", "#9FB798")) +
  labs(x = "Dispersal potential (km/y)", y = "Climate velocity (km/y)") +
  theme(legend.position = "none")

ggsave(cv_vs_dp_zoomzoom, path = "figures/sotm", filename = "cv_vs_dp_zoomzoom.png", 
       device = "png", height = 3.5, width = 4.5)

## make a bar plot
bar_group <- lags %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  mutate(expect_limitation = ifelse(ClimVeloTKmY <= AnnualDispPotKmY, "No",
                                    "Yes")) %>%
  mutate(expect_limitation = factor(expect_limitation, levels = c("Yes", "No"), 
                                    ordered = TRUE)) %>%
  ggplot(aes(x = expect_limitation, fill = Gradient)) +
  geom_bar() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(Group~Gradient) +
  scale_fill_manual(values = c("#9FB798", "#28587B")) +
  labs(x = "Expect dispersal limitation?", y = "Number of range expansion observations") +
  theme(legend.position = "none")

ggsave(bar_group, path = "figures/sotm", filename = "bar_group.png", 
       device = "png", height = 3.5, width = 4.5)

## get numbers 
lags %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  mutate(expect_limitation = ifelse(ClimVeloTKmY <= AnnualDispPotKmY, "No",
                                    "Yes")) %>%
  mutate(expect_limitation = factor(expect_limitation, levels = c("Yes", "No"), 
                                    ordered = TRUE)) %>%
  group_by(Group, Gradient, expect_limitation) %>%
  tally()

bar_gradient <- lags %>%
  mutate(expect_limitation = ifelse(ClimVeloTKmY <= AnnualDispPotKmY, "No",
                                    "Yes")) %>%
  ggplot(aes(x = expect_limitation, fill = Gradient)) +
  geom_bar() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_grid(~Gradient) +
  scale_fill_manual(values = c("#28587B", "#9FB798")) +
  labs(x = "Expect dispersal limitation?", y = "Number of observed expansions") +
  theme(legend.position = "none")

ggsave(bar_gradient, path = "figures/sotm", filename = "bar_gradient.png", 
       device = "png", height = 3.5, width = 4.5)



## now plot the real points 
lat <- filter(lags, Gradient == "Latitudinal") %>%
  mutate(Group = factor(Group, ordered = F))

lat$ClimVeloTKmY_og <- lat$ClimVeloTKmY
hist(lat$ClimVeloTKmY)
q = quantile(lat$ClimVeloTKmY, probs = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
lat$ClimVeloTKmY = cut(lat$ClimVeloTKmY,
                       breaks = q,
                       include.lowest = T)
plot(lat$ClimVeloTKmY)

lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\[", "(") 
lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\]", ")") 

mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)

# lat <- lat %>%
#   mutate(quant_max = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,2], "\\)", " "))) %>%
#   mutate(quant_min = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,1], "\\(", " "))) %>%
#   mutate(quant_mean = (quant_min + quant_max)/2) %>%
#   mutate(ClimVeloTKmY = paste(quant_min, " km/y - ", quant_max," km/y" , sep = ""))
# 
# obs_lat <- lat %>%
#   filter(ShiftR > 0.0001) %>%
#   ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
#   theme_bw() + 
#   geom_point() + 
#   facet_wrap(~ClimVeloTKmY, nrow = 2) +
#   labs(x = "Potential dispersal rate (km/y)", y = "Observed range expansion rate (km/y)",
#        colour = "") +
#   scale_colour_manual(values = mycolours) +
#   scale_y_continuous(limits = c(0, 40)) +
#   geom_vline(aes(xintercept = quant_min)) +
#   geom_hline(aes(yintercept = quant_min)) +
#   #geom_abline(intercept = 0, slope = 1, colour = "black") + 
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 20, hjust = 1)) +
#   theme(panel.grid = element_blank())
# 
# ggsave(obs_lat, path = "figures/sotm", filename = "observations_lat.png", 
#        device = "png", height = 3.5, width = 8.5)
# 
# 
# obs_lat_zoom <- lat %>%
#   filter(ShiftR > 0.0001) %>%
#   ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
#   theme_bw() + 
#   geom_point() + 
#   facet_wrap(~ClimVeloTKmY, nrow = 2) +
#   labs(x = "Potential dispersal rate (km/y)", y = "Observed range expansion rate (km/y)",
#        colour = "") +
#   scale_colour_manual(values = mycolours) +
#   scale_y_continuous(limits = c(0, 10)) +
#   scale_x_continuous(limits = c(0, 7)) +
#   geom_vline(aes(xintercept = quant_min)) +
#   geom_hline(aes(yintercept = quant_min)) +
#   #geom_abline(intercept = 0, slope = 1, colour = "black") + 
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 20, hjust = 1)) +
#   theme(panel.grid = element_blank())
# 
# 
# 
# ggsave(obs_lat_zoom, path = "figures/sotm", filename = "observations_lat_zoom.png", 
#        device = "png", height = 3.5, width = 8.5)
# 
# 
# disp_lim_lat <- lat %>%
#   filter(Gradient == "Latitudinal") %>%
#   filter(ShiftR > 0.0001) %>%
#   filter(ClimVeloTKmY_og > AnnualDispPotKmY) %>%
#   ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
#   theme_bw() + 
#   geom_point() + 
#   facet_wrap(~ClimVeloTKmY, nrow = 2) +
#   labs(x = "Potential dispersal rate (km/y)", y = "Observed range expansion rate (km/y)",
#        colour = "") +
#   scale_colour_manual(values = mycolours) +
#   #scale_y_continuous(limits = c(0, 40)) +
#   #scale_x_continuous(limits = c(0, 5)) +
#   geom_vline(aes(xintercept = quant_max)) +
#   geom_hline(aes(yintercept = quant_max)) +
#   #geom_abline(intercept = 0, slope = 1, colour = "black") + 
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 20, hjust = 1)) +
#   theme(panel.grid = element_blank()) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") 
# 
# ggsave(disp_lim_lat, path = "figures/sotm", filename = "disp_lim_lat.png", 
#        device = "png", height = 3.5, width = 8.5)
# 
# not_disp_lim_lat <- lat %>%
#   filter(ShiftR > 0.0001) %>%
#   filter(ClimVeloTKmY_og <= AnnualDispPotKmY) %>%
#   ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
#   theme_bw() + 
#   geom_point() + 
#   facet_wrap(~ClimVeloTKmY, nrow = 2) +
#   labs(x = "Potential dispersal rate (km/y)", y = "Observed range expansion rate (km/y)",
#        colour = "") +
#   scale_colour_manual(values = mycolours) +
#   scale_y_continuous(limits = c(0, 40)) +
#   geom_vline(aes(xintercept = quant_max)) +
#   geom_hline(aes(yintercept = quant_max)) +
#   #geom_abline(intercept = 0, slope = 1, colour = "black") + 
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 20, hjust = 1)) +
#   theme(panel.grid = element_blank())
# 
# ggsave(not_disp_lim_lat, path = "figures/sotm", filename = "not_disp_lim_lat.png", 
#        device = "png", height = 3.5, width = 8.5)
# 
# ## elevation
# ele <- filter(lags, Gradient == "Elevation") %>%
#   mutate(Group = factor(Group, ordered = F))
# 
# ele$ClimVeloTKmY_og <- ele$ClimVeloTKmY
# hist(ele$ClimVeloTKmY)
# q = quantile(ele$ClimVeloTKmY, probs = c(0,0.2,0.4,0.6,0.8,1))
# ele$ClimVeloTKmY = cut(ele$ClimVeloTKmY,
#                        breaks = q,
#                        include.lowest = T)
# plot(ele$ClimVeloTKmY)
# 
# ele$ClimVeloTKmY <- str_replace_all(ele$ClimVeloTKmY, "\\[", "(") 
# ele$ClimVeloTKmY <- str_replace_all(ele$ClimVeloTKmY, "\\]", ")") 
# 
# ele <- ele %>%
#   mutate(quant_max = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,2], "\\)", " "))) %>%
#   mutate(quant_min = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,1], "\\(", " "))) %>%
#   mutate(quant_mean = (quant_min + quant_max)/2) 
# 
# ele <- ele %>%
#   mutate(quant_max = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,2], "\\)", " "))) %>%
#   mutate(quant_min = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,1], "\\(", " "))) %>%
#   mutate(quant_mean = (quant_min + quant_max)/2) %>%
#   mutate(ClimVeloTKmY = paste(quant_min, " km/y - ", quant_max," km/y" , sep = "")) 
# 
# 
# obs_ele <- ele %>%
#   filter(ShiftR > 0.0001) %>%
#   ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
#   theme_bw() + 
#   geom_point() + 
#   facet_wrap(~ClimVeloTKmY, nrow = 1) +
#   labs(x = "Potential dispersal rate (km/y)", y = "Observed range expansion rate (km/y)",
#        colour = "") +
#   scale_colour_manual(values = mycolours) +
#   geom_vline(aes(xintercept = quant_min)) +
#   geom_hline(aes(yintercept = quant_min)) +
#   #geom_abline(intercept = 0, slope = 1, colour = "black") + 
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 20, hjust = 1)) +
#   theme(panel.grid = element_blank())
# 
# ggsave(obs_ele, path = "figures/sotm", filename = "observations_ele.png", 
#        device = "png", height = 2.5, width = 10.5)
# 
# 
# disp_lim_ele <- ele %>%
#   filter(ShiftR > 0.0001) %>%
#   filter(ClimVeloTKmY_og > AnnualDispPotKmY) %>%
#   ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
#   theme_bw() + 
#   geom_point() + 
#   facet_wrap(~ClimVeloTKmY, nrow = 1) +
#   labs(x = "Potential dispersal rate (km/y)", y = "Observed range expansion rate (km/y)",
#        colour = "") +
#   scale_colour_manual(values = mycolours) +
#   #scale_y_continuous(limits = c(0, 40)) +
#   #scale_x_continuous(limits = c(0, 5)) +
#   geom_vline(aes(xintercept = quant_max)) +
#   geom_hline(aes(yintercept = quant_max)) +
#   #geom_abline(intercept = 0, slope = 1, colour = "black") + 
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 20, hjust = 1)) +
#   theme(panel.grid = element_blank()) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") 
# 
# ggsave(disp_lim_ele, path = "figures/sotm", filename = "disp_lim_ele.png", 
#        device = "png", height = 2.5, width = 10.5)
# 
# not_disp_lim_ele <- ele %>%
#   filter(ShiftR > 0.0001) %>%
#   filter(ClimVeloTKmY_og <= AnnualDispPotKmY) %>%
#   ggplot(., aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) + 
#   theme_bw() + 
#   geom_point() + 
#   facet_wrap(~ClimVeloTKmY, nrow = 1) +
#   labs(x = "Potential dispersal rate (km/y)", y = "Observed range expansion rate (km/y)",
#        colour = "") +
#   scale_colour_manual(values = mycolours) +
#   geom_vline(aes(xintercept = quant_max)) +
#   geom_hline(aes(yintercept = quant_max)) +
#   #geom_abline(intercept = 0, slope = 1, colour = "black") + 
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 20, hjust = 1)) +
#   theme(panel.grid = element_blank())
# 
# ggsave(not_disp_lim_ele, path = "figures/sotm", filename = "not_disp_lim_ele.png", 
#        device = "png", height = 2.5, width = 10.5)
# 
# 
# lags %>%
#   mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
#   ggplot(aes(x = ClimVeloTKmY, fill = ..x..)) +
#   geom_histogram() +
#   theme_bw() +
#   facet_grid(~Gradient) 
# +
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         panel.spacing = unit(2 , "lines")) +
#   labs(x = "Potential dispersal rate (km/y)",
#        y = "Observed range expansion rate (km/y)", 
#        colour = "") +
#   theme(legend.position = "none")+
#   scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 



#### hypothesis testing 
data_unlog <- lags %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) +
  geom_point() +
  theme_bw() +
  facet_grid(~Gradient) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0), 
                                        add = c(20, 0)),
                     limits = c(0, 1400)) +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0),
                     labels = c("    0", "   10", " 20", "   30", "   40")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none")+
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 

ggsave(data_unlog, path = "figures/sotm", filename = "data-unlog.png", 
       device = "png", height = 4, width = 8)

data_log <- lags %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
    labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
    limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0),
                     labels = c("    0", "   10", " 20", "   30", "   40")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 

ggsave(data_log, path = "figures/sotm", filename = "data-log.png", 
       device = "png", height = 4, width = 8)

## make elevation one with diff scale 
data_ele_log <- lags %>%
  filter(ShiftR > 0.0001) %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)  + 
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) 

ggsave(data_ele_log, path = "figures/sotm", filename = "data-ele-log.png", 
       device = "png", height = 2.8, width = 3.1)


## now, plot line representing max climate velocity across elev and lat:
data_log <- data_log +
  geom_vline(aes(xintercept = max(ClimVeloTKmY)))

ggsave(data_log, path = "figures/sotm", filename = "data-log-line.png", 
       device = "png", height = 4, width = 8)

data_ele_log <- data_ele_log +
  geom_vline(aes(xintercept = max(ClimVeloTKmY)))

ggsave(data_ele_log, path = "figures/sotm", filename = "data-ele-log-line.png", 
       device = "png", height = 4, width = 4.2)

lags %>%
  filter(ShiftR > 0.0001) %>%
  filter(Gradient == "Latitudinal") %>%
  filter(AnnualDispPotKmY >= ClimVeloTKmY) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = Gradient)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  theme(panel.grid = element_blank()) +
  # scale_x_continuous(expand = expansion(mult = c(0, 0), 
  #                                       add = c(20, 0)),
  #                    limits = c(0, 1400)) +
  scale_x_log10(limits = c(0.001, 1400)) +
  scale_y_log10() +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0)) +
  scale_colour_manual(labels = c("Elevation study", 
                                 "Latitudinal study"),
                      values = c("#9FB798")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") + 
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  geom_vline(aes(xintercept = max(ClimVeloTKmY)))

ggsave(path = "figures/sotm", filename = "data-lat-log-gt.png", 
       device = "png", height = 3.5, width = 4.5)


## model data 
notlimited_lat <- lags %>%
  filter(ShiftR > 0.0001) %>%
  filter(Gradient == "Latitudinal") %>%
  filter(AnnualDispPotKmY >= ClimVeloTKmY) %>%
  mutate(ClimVeloTKmY = factor(ClimVeloTKmY, levels = unique(.$ClimVeloTKmY),
                               ordered = TRUE))

limited_lat <- lags %>%
  filter(ShiftR > 0.0001) %>%
  filter(Gradient == "Latitudinal") %>%
  filter(AnnualDispPotKmY < ClimVeloTKmY)  %>%
  mutate(ClimVeloTKmY = factor(ClimVeloTKmY, levels = unique(.$ClimVeloTKmY),
                               ordered = TRUE))


library(nlme)
mod_nl_lat <- lme(ShiftKmY ~ AnnualDispPotKmY, random = ~1|ClimVeloTKmY, 
                 data = notlimited_lat)
summary(mod_nl_lat)

R2m=r.squaredGLMM(mod_nl_lat)[[1]] 
R2c=r.squaredGLMM(mod_nl_lat)[[2]]
R2m
R2c

mod_l_lat <- lme(ShiftKmY ~ AnnualDispPotKmY, random = ~1|ClimVeloTKmY, 
                 data = limited_lat)
summary(mod_l_lat)

R2m=r.squaredGLMM(mod_l_lat)[[1]] 
R2c=r.squaredGLMM(mod_l_lat)[[2]]
R2m
R2c

mod_l_lat_lm <- lm(ShiftKmY ~ AnnualDispPotKmY, 
                data = limited_lat)
summary(mod_l_lat_lm)

R2m=r.squaredGLMM(mod_l_lat_lm)[[1]] 
R2m

AIC(mod_l_lat_lm, 
    mod_l_lat) ## this is better fit 

## elevation
notlimited_ele <- lags %>%
  filter(ShiftR > 0.0001) %>%
  filter(Gradient == "Elevation") %>%
  filter(AnnualDispPotKmY >= ClimVeloTKmY) %>%
  mutate(ClimVeloTKmY = factor(ClimVeloTKmY, levels = unique(.$ClimVeloTKmY),
                               ordered = TRUE))

limited_ele <- lags %>%
  filter(ShiftR > 0.0001) %>%
  filter(Gradient == "Elevation") %>%
  filter(AnnualDispPotKmY < ClimVeloTKmY)  %>%
  mutate(ClimVeloTKmY = factor(ClimVeloTKmY, levels = unique(.$ClimVeloTKmY),
                               ordered = TRUE))


library(nlme) 
mod_nl_ele <- lme(ShiftKmY ~ AnnualDispPotKmY, random = ~1|ClimVeloTKmY,
                  data = notlimited_ele)
summary(mod_nl_ele)

R2m=r.squaredGLMM(mod_nl_ele)[[1]]
R2c=r.squaredGLMM(mod_nl_ele)[[2]]
R2m
R2c

mod_l_ele <- lme(ShiftKmY ~ AnnualDispPotKmY, random = ~1|ClimVeloTKmY,
                 data = limited_ele)
summary(mod_l_ele)

R2m=r.squaredGLMM(mod_l_ele)[[1]]
R2c=r.squaredGLMM(mod_l_ele)[[2]]
R2m
R2c

mod_l_ele_lm <- lm(ShiftKmY ~ AnnualDispPotKmY,
                   data = limited_ele)
summary(mod_l_ele_lm)

R2m=r.squaredGLMM(mod_l_ele_lm)[[1]]
R2m

AIC(mod_l_ele_lm, ## this is better fit
    mod_l_ele)

# function to find mode of each predictor
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## plot predictions  
new_data <- expand.grid(AnnualDispPotKmY = seq(min(notlimited_lat$AnnualDispPotKmY),
                                               max(notlimited_lat$AnnualDispPotKmY), 
                                               by = 0.01))

pred_nl_lat <- predict(mod_nl_lat, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_nl_lat <- new_data %>%
  mutate(pred_expansion = pred_nl_lat$fit) %>%
  mutate(pred_expansion_SE = pred_nl_lat$se.fit)


nl_lat_nopred <- lags %>%
  filter(AnnualDispPotKmY >= ClimVeloTKmY) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0),
                     labels = c("    0", "   10", " 20", "   30", "   40")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 

nl_lat_pred <- nl_lat_nopred +
  geom_ribbon(data = fitted_pred_nl_lat, 
              aes(x = AnnualDispPotKmY, 
                  ymax = pred_expansion + pred_expansion_SE,
                  ymin = pred_expansion - pred_expansion_SE),
              inherit.aes = FALSE,
              colour = "transparent", 
              fill = "grey", 
              alpha = 0.5) + # add SE
  geom_path(data = fitted_pred_nl_lat, 
              aes(x = AnnualDispPotKmY, y = pred_expansion),
              inherit.aes = FALSE,
              size = 1,
              colour = "black") # add prediction

ggsave(nl_lat_nopred, path = "figures/sotm", filename = "data-notlimited-lat-nopredictions.png", 
       device = "png", height = 4, width = 8)
ggsave(nl_lat_pred, path = "figures/sotm", filename = "data-notlimited-lat-predictions.png", 
       device = "png", height = 4, width = 8)


## now ele
new_data <- expand.grid(AnnualDispPotKmY = seq(min(notlimited_ele$AnnualDispPotKmY),
                                               max(notlimited_ele$AnnualDispPotKmY), 
                                               by = 0.01))

pred_nl_ele <- predict(mod_nl_ele, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_nl_ele <- new_data %>%
  mutate(pred_expansion = pred_nl_ele$fit) %>%
  mutate(pred_expansion_SE = pred_nl_ele$se.fit)

nl_ele_nopred <- lags %>%
  filter(ShiftR > 0.0001) %>%
  filter(Gradient == "Elevation") %>%
  filter(AnnualDispPotKmY >= ClimVeloTKmY) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) + 
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none")  +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 

nl_ele_pred <- nl_ele_nopred +
  geom_ribbon(data = fitted_pred_nl_ele, 
              aes(x = AnnualDispPotKmY, 
                  ymax = pred_expansion + pred_expansion_SE,
                  ymin = pred_expansion - pred_expansion_SE),
              inherit.aes = FALSE,
              colour = "transparent", 
              fill = "grey", 
              alpha = 0.5) + # add SE
  geom_path(data = fitted_pred_nl_ele, 
            aes(x = AnnualDispPotKmY, y = pred_expansion),
            inherit.aes = FALSE,
            size = 1,
            colour = "black") # add prediction

ggsave(nl_ele_nopred, path = "figures/sotm", filename = "data-notlimited-ele-nopredictions.png", 
       device = "png", height = 2.8, width = 3.1)
ggsave(nl_ele_pred, path = "figures/sotm", filename = "data-notlimited-ele-predictions.png", 
       device = "png", height = 2.8, width = 3.1)




new_data <- expand.grid(AnnualDispPotKmY = seq(min(limited_lat$AnnualDispPotKmY),
                                               max(limited_lat$AnnualDispPotKmY), 
                                               by = 0.0001))

pred_l_lat <- predict(mod_l_lat, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_l_lat <- new_data %>%
  mutate(pred_expansion = pred_l_lat$fit) %>%
  mutate(pred_expansion_SE = pred_l_lat$se.fit)


l_lat_nopred <- lags %>%
  filter(AnnualDispPotKmY < ClimVeloTKmY) %>%
  mutate(Gradient = factor(Gradient, levels = c("Latitudinal", "Elevation"), ordered = TRUE)) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  facet_grid(~Gradient) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines")) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 40), expand = c(0,0),
                     labels = c("    0", "   10", " 20", "   30", "   40")) +
  labs(x = "Potential dispersal rate (km/y)",
       y = "Observed range expansion rate (km/y)", 
       colour = "") +
  theme(legend.position = "none") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 

l_lat_pred <- l_lat_nopred +
  geom_ribbon(data = fitted_pred_l_lat, 
              aes(x = AnnualDispPotKmY, 
                  ymax = pred_expansion + pred_expansion_SE,
                  ymin = pred_expansion - pred_expansion_SE),
              inherit.aes = FALSE,
              colour = "transparent", 
              fill = "grey", 
              alpha = 0.5) + # add SE
  geom_path(data = fitted_pred_l_lat, 
            aes(x = AnnualDispPotKmY, y = pred_expansion),
            inherit.aes = FALSE,
            size = 1,
            colour = "black") # add prediction

ggsave(l_lat_nopred, path = "figures/sotm", filename = "data-limited-lat-nopredictions.png", 
       device = "png", height = 4, width = 8)
ggsave(l_lat_pred, path = "figures/sotm", filename = "data-limited-lat-predictions.png", 
       device = "png", height = 4, width = 8)


## now ele
new_data <- expand.grid(AnnualDispPotKmY = seq(min(limited_ele$AnnualDispPotKmY),
                                               max(limited_ele$AnnualDispPotKmY), 
                                               by = 0.0001))

pred_l_ele <- predict(mod_l_ele, new_data, level = 0, se.fit = T, re.form = NA)

fitted_pred_l_ele <- new_data %>%
  mutate(pred_expansion = pred_l_ele$fit) %>%
  mutate(pred_expansion_SE = pred_l_ele$se.fit)

l_ele_nopred <- lags %>%
  filter(ShiftR > 0.0001) %>%
  filter(Gradient == "Elevation") %>%
  filter(AnnualDispPotKmY < ClimVeloTKmY) %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  stat_function(colour = "black", fun = function(x){x},
                linetype = "dashed") + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2 , "lines"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                limits = c(0.0001, 1400)) +
  scale_y_continuous(limits = c(0, 0.04), expand = c(0,0)) + 
  labs(x = "",
       y = "", 
       colour = "") +
  theme(legend.position = "none") +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) 

l_ele_pred <- l_ele_nopred +
  geom_ribbon(data = fitted_pred_l_ele, 
              aes(x = AnnualDispPotKmY, 
                  ymax = pred_expansion + pred_expansion_SE,
                  ymin = pred_expansion - pred_expansion_SE),
              inherit.aes = FALSE,
              colour = "transparent", 
              fill = "grey", 
              alpha = 0.5) + # add SE
  geom_path(data = fitted_pred_l_ele, 
            aes(x = AnnualDispPotKmY, y = pred_expansion),
            inherit.aes = FALSE,
            size = 1,
            colour = "black") # add prediction

ggsave(l_ele_nopred, path = "figures/sotm", filename = "data-limited-ele-nopredictions.png", 
       device = "png", height = 2.8, width = 3.1)
ggsave(l_ele_pred, path = "figures/sotm", filename = "data-limited-ele-predictions.png", 
       device = "png", height = 2.8, width = 3.1)

