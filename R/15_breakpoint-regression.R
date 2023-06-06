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
####     subset data, convert units     ####
############################################
## use max dispersal potential per species
v1$dispersal_potential_kmY = v1$MaxDispersalPotentialKmY
v1$dispersal_potential_mY = v1$MaxDispersalPotentialmY

mod_data <- v1 %>%
  ## filter out observations where climate velocity is negative at leading edge/optimum (expect contraction)
  filter(LatVeloT >= 0 | EleVeloT >= 0) %>%
  ## get rid of negative shifts 
  filter(ShiftR > 0) %>%
  ## get rid of trailing edge 
  filter(., Position != "Trailing edge") %>%
  ## make sure none have empty dispersal potential 
  filter(!is.na(dispersal_potential_kmY)) %>%
  ## make one column for annual dispersal potential, climate velo for easier plotting of lat x elev data together
  mutate(annual_dispersal_pot = ifelse(Gradient == "Elevation",
                                       dispersal_potential_mY,
                                       ifelse(Gradient == "Latitudinal",
                                              dispersal_potential_kmY,
                                              NA))) %>%
  ## reorder factors 
  mutate(Group = factor(Group, ordered = TRUE, levels = c("Birds", "Plants", "Mammals",
                                                          "Fish", "Amphibians", "Squamates"))) 

## convert units of latitudinal and elevation shifts & dispersal & climate velocity to km/y:
mod_data$ShiftKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$ShiftR / 1000, mod_data$ShiftR)
mod_data$ClimVeloTKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$EleVeloT / 1000, mod_data$LatVeloT)
mod_data$AnnualDispPotKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$annual_dispersal_pot / 1000, mod_data$annual_dispersal_pot)

mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)

mod_data <- select(mod_data, -c(SLDiff, CorrShift, PredSLShift, study_level_shift)) %>%
  distinct()

## make plots to summarize subset 
nrow(mod_data) # 1935 range shifts
length(unique(mod_data$Species)) # 409 species 

grad <- mod_data %>%
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

mod_data %>%
  group_by(Gradient) %>% tally()

groups <- mod_data %>%
  ggplot(aes(x = Group)) +
  geom_bar() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs (x = "", y = "Number of\nrange expansion\nobservations") 

ggsave(groups, path = "figures/sotm", filename = "barplot-groups.png", 
       device = "png", height = 2, width = 4)

mod_data %>%
  group_by(Group) %>% tally()


#################################################
####      fitting breakpoint regression      ####
#################################################
## filter to shifts > 0.0001
mod_data <- mod_data %>%
  filter(ShiftR > 0.0001) 

## split between elev and lat
lat = filter(mod_data, Gradient == "Latitudinal")
ele = filter(mod_data, Gradient == "Elevation")

## fit normal regression 
mod_lat <- lm(ShiftKmY ~ AnnualDispPotKmY, data = lat)

## fit breakpoint regression
mod_bp_lat <- segmented(mod_lat, 
                        seg.Z = ~ AnnualDispPotKmY,
                        psi = mean(lat$ClimVeloTKmY)) ## set break point to = mean climate velocity across sites

mod_bp_lat <- segmented(mod_lat, 
                        seg.Z = ~ AnnualDispPotKmY) ## do not set any starting value for break point

summary(mod_bp_lat)

## get estimated breakpoint and its standard error
mod_bp_lat$psi ## 4.78

## get the slopes
slope(mod_bp_lat) ## 0.615, -0.0003

## get the fitted data
my.fitted <- fitted(mod_bp_lat)
df <- data.frame(AnnualDispPotKmY = lat$AnnualDispPotKmY, ShiftKmY = my.fitted,
                 min = min(lat$ClimVeloTKmY), 
                 max = max(lat$ClimVeloTKmY),
                 mean = mean(lat$ClimVeloTKmY))

## plot the fitted model
df %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY)) + 
  geom_point(data = lat, aes(colour = ClimVeloTKmY), alpha = 0.5) + 
  geom_line() +
  theme_bw() +
  scale_x_log10() +
  stat_function(colour = "black", # add 1:1 line
                fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5)
# geom_vline(aes(xintercept = mean)) + # plot theoretical point 
# geom_rect(aes(xmin = min, xmax = max, ymin = 0, ymax =41), 
#           alpha = 0.003) 


## now bin by climate velocity and fit separate regression to each 
## choose 4 bins so that enough species with different dispersal abilities are sampled across different climate velocities
lat = filter(mod_data, Gradient == "Latitudinal")
hist(lat$ClimVeloTKmY)
q = quantile(lat$ClimVeloTKmY, probs = c(0,0.25, 0.5, 0.75,1))
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
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "") +
  scale_colour_manual(values = mycol) +
  stat_function(colour = "black", linetype = "dashed", fun = function(x){x}) 
## in some bins, probably not enough obs on either side of expected inflection point to find relationship

## now model
split <- split(lat, 
               f = lat$quant_mean)

mod_fits <- lapply(split, FUN = function(x) {
  mod_lat <- lm(ShiftKmY ~ AnnualDispPotKmY, data = x)
  
  fit = segmented(mod_lat,
                  seg.Z = ~ AnnualDispPotKmY)
  return(fit)
}
)

## extract break points  
results <- lapply(mod_fits, FUN = function(x) {
  x$psi
})

results <- as.data.frame(do.call(rbind, results))
results$quant_mean <- unique(lat$quant_mean)[order(unique(lat$quant_mean))]
results <- left_join(results, quants)

## plot theoretical versus real break points 
results %>% 
  ggplot(aes(y = quant_mean, x = ClimVeloTKmY)) +
  theme_bw() +
  labs(y = "Breakpoint", x = "Climate velocity") +
  geom_pointrange(aes(y = Est., x = ClimVeloTKmY, 
                      ymin = Est. - St.Err,
                      ymax = Est. + St.Err),
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
slopes <- left_join(slopes, quants)

## plot slopes 
slopes %>% 
  ggplot(aes(y = expected_slope, x = ClimVeloTKmY)) +
  theme_bw() +
  labs(y = "Slope", x = "Climate velocity", 
       colour = "Expected slope") +
  geom_pointrange(aes(y = Est., x = ClimVeloTKmY,
                      ymin = `CI(95%).l`,
                      ymax = `CI(95%).u`, 
                      colour = as.factor(expected_slope))) +
  scale_colour_discrete(c("red", "blue")) +
  geom_point() 


## get the fitted data
fitted <- lapply(mod_fits, FUN = function(x) {
  fitted(x)
})
real <- lapply(split, FUN = function(x) {
  x$AnnualDispPotKmY
})

fitted <- unlist(fitted)
real <- unlist(real)

reps = unlist(lapply(split, FUN = nrow))

df <- data.frame(AnnualDispPotKmY = real, ShiftKmY = fitted,
                 quant_mean = rep(unique(lat$quant_mean)[order(unique(lat$quant_mean))],
                                  reps),
                 quant_max = rep(unique(lat$quant_max)[order(unique(lat$quant_max))],
                                 reps),
                 quant_min = rep(unique(lat$quant_min)[order(unique(lat$quant_min))],
                                 reps),
                 ClimVeloTKmY = rep(unique(lat$ClimVeloTKmY)[order(unique(lat$ClimVeloTKmY))],
                                    reps))


## plot the fitted model
df %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY)) + 
  geom_point(data = lat, aes(colour = ClimVeloTKmY)) + 
  geom_line(aes(group = quant_max)) +
  theme_bw() +
  scale_x_log10() +
  facet_grid(~ClimVeloTKmY) +
  labs(x = "Annual dispersal potential (km/y)", y = "Observed range shift (km/y)",
       colour = "") +
  scale_colour_manual(values = mycol) +
  stat_function(colour = "black", # add 1:1 line
                fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
  geom_vline(aes(xintercept = quant_mean)) + # plot theoretical point 
  geom_rect(aes(xmin = quant_min, xmax = quant_max, ymin = 0, ymax =41), 
            alpha = 0.003)


