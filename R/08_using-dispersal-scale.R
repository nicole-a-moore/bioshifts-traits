## seeing whether empirical dispersal scale explains variation in range shifts
library(tidyverse)
library(PNWColors)
library(gridExtra)
library(grid)
library(cowplot)
library(lme4)
source("R/taxonomic-harmonization/clean_taxa_functions.R")


rs_data = read.table("data-raw/bioshifts-download/Lenoir_et_al/Analysis/Table_S1.csv",sep=";",h=T,dec=".",
                     stringsAsFactors = FALSE) 

## do some name fixes 
# rs_data$Species[which(!rs_data$Species %in% spv1$species)]

rs_data$Species = Clean_Names(rs_data$Species, return_gen_sps = F)

rs_data$Species[which(rs_data$Species == "Quercus x")] = "Quercus" 
rs_data$Species[which(rs_data$Species == "Mentha x")] = "Mentha" 
rs_data$Species[which(rs_data$Species == "Circaea x intermedia")] = "Circaea intermedia" 

rs_data$Species = Clean_Names(rs_data$Species, return_gen_sps = F)

#rs_data$Species[which(!rs_data$Species %in% spv1$species)]



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

#----------------------
#### get shift data for species with dispersal scale ####
## read in bioshifts v1
v1 = read.table("data-raw/bioshiftsv1/Shifts2018_checkedtaxo.txt",
                header = T,
                encoding="latin1") 

v1 %>% 
  filter(Type == "LAT")%>%
  group_by(Param) %>%
  summarise(plusplus = length(which(SHIFT > 0 & v.lat.mean> 0)),
            plusminus = length(which(SHIFT > 0 & v.lat.mean < 0)),
            minusplus = length(which(SHIFT < 0 & v.lat.mean> 0)),
            minusminus = length(which(SHIFT < 0 & v.lat.mean < 0))) 

v1 %>% 
  filter(Type == "ELE") %>%
  group_by(Param) %>%
  summarise(plusplus = length(which(SHIFT > 0 & v.ele.mean> 0)),
            plusminus = length(which(SHIFT > 0 & v.ele.mean < 0)),
            minusplus = length(which(SHIFT < 0 & v.ele.mean> 0)),
            minusminus = length(which(SHIFT < 0 & v.ele.mean < 0))) 


## clean names to make them match reported names in species list
v1$reported_name = Clean_Names(gsub("_", " ", v1$Publi_name), return_gen_sps = F)

spv1 <- filter(sp, v1 == 1) %>%
  select(reported_name, scientificName) %>%
  unique()

## yay! all are there
which(!v1$reported_name %in% spv1$reported_name)

## join to add scientific name column
v1 = left_join(v1, spv1)

## subset to species with dispersal scale 
v1 <- filter(v1, scientificName %in% dscale$scientificName)
length(unique(v1$scientificName)) #693 species 


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
length(unique(v2$scientificName)) #652 species 

## make sure all the species are there:
length(unique(c(v2$scientificName, v1$scientificName))) # 725 unique species 
length(unique(dscale$scientificName)) # 725 unique species 

#----------------------
### relate velocity of shift to dispersal scale ###
colnames(v1)
colnames(v2)

#----------------------
## clean bioshifts v1 

#--------------------------
#use insects orders to avoid mixing butterflies with others
dscale$class[dscale$class=="Insecta"] <- dscale$order[dscale$class=="Insecta"]
v1$Class[v1$Class=="Insecta"] <- v1$Order[v1$Class=="Insecta"]

#----------------------
#correct and reorganize methodological variables
v1$Data[v1$Data=="occurence-based"] <- "occurrence-based"
v1$Sampling <- ifelse(v1$Sampling == "TWO","TWO","MULTIPLE")
v1$Grain_size <- ifelse(v1$Grain_size %in% c("large","very_large"),"large",v1$Grain_size)
v1$Uncertainty_Distribution <- ifelse(v1$Uncertainty_Distribution %in% c("RESAMPLING","RESAMPLING(same)"),"RESAMPLING",
                                         ifelse(v1$Uncertainty_Distribution %in% c("MODEL","MODEL+RESAMPLING(same)","RESAMPLING+MODEL"),"MODEL",
                                                ifelse(v1$Uncertainty_Distribution %in% c("DETECTABILITY","RESAMPLING(same)+DETECTABILITY"),"DETECTABILITY",
                                                       v1$Uncertainty_Distribution)))
v1$Uncertainty_Distribution <- ifelse(v1$Uncertainty_Distribution == "RAW","OPPORTUNISTIC","PROCESSED")  

#transform study area
v1$ID.area <- log(v1$ID.area)


#----------------------
#remove freshwater fishes
#v1 <- v1[-which(v1$Class == "Actinopterygii" & v1$ECO=="T"),]

#remove marine birds
#v1 <- v1[-which(v1$Class == "Aves" & v1$ECO=="M"),]

#----------------------
## add dispersal scale 
## get rid of old taxonomy columns from v1 (they aren't right)
v1 <- select(v1, -c("Kingdom", "Phylum", "Class", "Order", "Family"))
v1_saved = v1

## after everything, there should be:
length(which(v1$Type == "ELE")) ## 1490 elevation shifts
length(which(v1$Type == "LAT")) ## 2312 latitude shifts

## get rid of columns that will cause duplication in dispersal scale database 
dscale <- select(dscale, -c("reported_name", "reported_name_fixed", "db", "db_code")) %>%
  unique()

v1 = left_join(v1, dscale) 

## check on merge
length(which(is.na(v1$DispersalDistanceKm))) #0 missing dispersal scale
length(unique(v1$scientificName)) #still 693 species


#----------------------
## calculate one dispersal distance value per species per metric
## ex. one mean, one max
# dscale %>%
#   ggplot(aes(x = Code)) + geom_bar() + coord_flip()

v1 = v1 %>% 
  group_by(scientificName, Code) %>%
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
  group_by(scientificName) %>%
  mutate(DispersalDistanceKm_max = max(DispersalDistanceKm)) %>%
  select(Code, scientificName, DispersalDistanceKm, DispersalDistanceKm_unique, DispersalDistanceKm_max, 
         everything()) %>%
  ungroup()

## look at the metrics that were larger than max 
max = v1 %>%
  filter(Code == "MaxDispersalDistance") %>%
  mutate(matches = ifelse(DispersalDistanceKm_unique == DispersalDistanceKm_max, "Y", "N")) %>%
  select(matches, scientificName) %>%
  unique() %>%
  filter(matches == "N")

length(unique(max$scientificName))
## 19 species have other dispersal distance metrics that are larger than their so-called maximum 

v1 %>%
  filter(scientificName %in% max$scientificName) %>%
  filter(DispersalDistanceKm == DispersalDistanceKm_max) %>% 
  select(scientificName, Code) %>% 
  unique() 
## most are mean values, some median
## mostly fish, 1 plant, 1 tree, some small mammals

## correct max dispersal distance for these species
v1 = v1 %>%
  group_by(scientificName, Code) %>% 
  mutate(DispersalDistanceKm_unique = ifelse(DispersalDistanceKm_unique != DispersalDistanceKm_max &
                                               Code == "MaxDispersalDistance", 
                                             DispersalDistanceKm_max,
                                             DispersalDistanceKm_unique)) %>%
  ungroup() %>%
  select(-DispersalDistanceKm_max)

v1 %>%
  select(Code, scientificName) %>%
  unique() %>%
  group_by(Code) %>%
  tally()

## make sure we don't drop any data
maxsp = unique(v1$scientificName[which(v1$Code == "MaxDispersalDistance")])
nomaxsp <- unique(v1$scientificName[which(!v1$scientificName %in% maxsp)])

## 589 species have max dispersal distance
## start with max!
v1 <- v1 %>%
  filter(Code == "MaxDispersalDistance") %>%
  select(-DispersalDistanceKm, -DispersalDistance, -Unit, -Field, -Source, -Database, 
         -Sex) %>%
  unique()

which(nomaxsp %in% v1$scientificName)
which(maxsp %in% v1$scientificName)

check = v1_saved %>%
  filter(scientificName %in% maxsp)
nrow(check) == nrow(v1)
## yay = all observations are there without any duplication

v1 <- rename(v1, "MaxDispersalDistanceKm" = DispersalDistanceKm_unique)

## save dataset 
#write.csv(v1, "data-processed/bioshiftsv1_max-dispersal-distance.csv", row.names = FALSE)
v1 = read.csv("data-processed/bioshiftsv1_max-dispersal-distance.csv")

### join age at maturity data 
am <- read.csv("data-processed/age-at-maturity.csv")

am %>% 
  ggplot(aes(x = AgeAtMaturity, fill = class)) + geom_histogram() +
  theme_bw()

## convert units 
unique(am$Unit)

am <- am %>%
  mutate(AgeAtMaturityDays = ifelse(Unit %in% c("yrs", "y", "years", "year"), 
                                    365*AgeAtMaturity,
                                    ifelse(Unit == "weeks", 
                                           7*AgeAtMaturity,
                                           ifelse(Unit == "months",
                                                  30*AgeAtMaturity, 
                                                  AgeAtMaturity)))) %>%
  mutate(GenerationLengthDays = ifelse(Unit %in% c("yrs", "y", "years", "year"), 
                                  365*GenerationLength,
                                  ifelse(Unit == "weeks", 
                                         7*GenerationLength,
                                         ifelse(Unit == "months",
                                                30*GenerationLength, 
                                                GenerationLength))))

am %>% 
  filter(!is.na(AgeAtMaturityDays)) %>%
  ggplot(aes(x = log(AgeAtMaturityDays), fill = class)) + geom_histogram() +
  theme_bw()


## if multiple estimates of age at maturity per species, keep the lowest 
am_join <- am %>%
  group_by(scientificName) %>%
  mutate(AgeAtMaturityDays = min(AgeAtMaturityDays)) %>% # select minimum per species 
  ungroup() %>%
  select(scientificName, AgeAtMaturityDays) %>%
  unique() %>%
  mutate(YearOfMaturity = ceiling(AgeAtMaturityDays/365)) ## make new column for a value that's rounded to the nearest year 

## join to dispersal data:
v1 <- left_join(v1, am_join)

length(unique(v1$scientificName)) ## still have all the species!
length(unique(v1$scientificName[which(is.na(v1$AgeAtMaturityDays))])) 
## 162 / 589 species do not have age at maturity data 

unique(v1$scientificName[which(is.na(v1$AgeAtMaturityDays))])

## calculate dispersal potential for species with age at maturity/longevity 
v1 = v1 %>%
  mutate(MaxDispersalPotentialKmY = ifelse(!is.na(YearOfMaturity), 
                                     MaxDispersalDistanceKm/YearOfMaturity,
                                     NA)) %>%
  mutate(MaxDispersalPotentialmY = ifelse(!is.na(YearOfMaturity), 
                                           (MaxDispersalDistanceKm*1000)/YearOfMaturity,
                                           NA)) %>%
  mutate(MaxDispersalDistancem = MaxDispersalDistanceKm*1000)

ggplot(v1, aes(x = log(MaxDispersalPotentialKmY), fill = class)) + geom_histogram()
ggplot(v1, aes(x = log(MaxDispersalPotentialmY), fill = class)) + geom_histogram()
ggplot(v1, aes(x = log(MaxDispersalDistanceKm), fill = class)) + geom_histogram()

## which taxa are missing age at maturity data?
v1 %>%
  filter(is.na(YearOfMaturity)) %>%
  ggplot(., aes(x = log(MaxDispersalDistanceKm), fill = class)) + geom_histogram()
## mostly magnoliopsida and mammals


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
ele <- v1[v1$Type == "ELE",]
ele <- ele[which(is.na(ele$v.ele.mean)==F),]
ele$lags <- ele$SHIFT - ele$v.ele.mean
hist(ele$lags)

#latitude
lat <- v1[v1$Type == "LAT",]
lat$lags <- lat$SHIFT - lat$v.lat.mean
hist(lat$lags)

## how many shifts for species with maximum dispersal distance in v1? 
nrow(ele) + nrow(lat) #4400

## how many leading edge shifts for species with maximum dispersal distance in v1? 
length(which(ele$Param == "LE")) + length(which(lat$Param == "LE")) #1367

lat = lat %>%
  filter(v.lat.mean >= 0) %>%
  ## get rid of cases where you expect contraction at leading edge (negative velocity)
  mutate(expect_tracking_dist = ifelse(MaxDispersalDistanceKm > log(v.lat.mean),"Yes", "No")) %>%
  mutate(expect_tracking_pot = ifelse(MaxDispersalPotentialKmY > log(v.lat.mean),"Yes", "No")) %>%
  ## get rid of trailing edge 
  filter(., Param != "TE") %>% 
  filter(!is.na(MaxDispersalPotentialKmY))

ele = ele %>%
  filter(v.ele.mean >= 0) %>%
  ## get rid of cases where you expect contraction at leading edge (negative velocity)
  mutate(expect_tracking_dist = ifelse(MaxDispersalDistancem > log(v.ele.mean),"Yes", "No")) %>%
  mutate(expect_tracking_pot = ifelse(MaxDispersalPotentialmY > log(v.ele.mean),"Yes", "No"))  %>%
  ## get rid of trailing edge 
  filter(., Param != "TE") %>%
  filter(!is.na(MaxDispersalPotentialKmY))


lat %>% rbind(., ele) %>%
  ggplot(aes(x = Param)) + geom_bar() + theme_bw() +
  labs(y = "Count", x = "Range shift parameter") + 
  facet_wrap(~Type)

lat %>% rbind(., ele) %>%
  group_by(Type, Param) %>%
  tally()

# Plot:
#----------------------
#prepare a dataset for each range shift parameter
ele_te <- ele[ele$Param=="TE",]
ele_le <- ele[ele$Param=="LE",]
ele_o <- ele[ele$Param=="O",]
lat_te <- lat[lat$Param=="TE",]
lat_le <- lat[lat$Param=="LE",]
lat_o <- lat[lat$Param=="O",]

pal = pnw_palette("Bay",7)


## how many species should be able to keep up with climate change across elev and lat?
pol1 <- data.frame(x = c(-10, 10, 10), y = c(-10, -10, 10))
pol2 <- data.frame(x = c(-10, -10, 10), y = c(-10, 10, 10))
pol3 <- data.frame(x = c(-4, 15, 15), y = c(-4, -4, 15))
pol4 <- data.frame(x = c(-4, -4, 15), y = c(-4, 15, 15))

## try with dispersal potential 
plot3 = lat %>% 
  ggplot(aes(y = MaxDispersalPotentialKmY, x = log(v.lat.mean))) + 
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
  ggplot(aes(y = MaxDispersalPotentialmY, x = log(v.ele.mean))) +
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
       filename = "maxpotential-versus-lag.png")

ele %>%
  group_by(expect_tracking_pot) %>%
  tally()

lat %>%
  group_by(expect_tracking_pot) %>%
  tally()

lat %>%
  rbind(., ele) %>%
  ggplot(aes(x = MaxDispersalDistanceKm, fill = class)) + geom_histogram() +
  theme_bw() + facet_wrap(~Type)

## take a closer look 
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
       filename = "ele-le_points_potential.png")

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
       filename = "ele-le_boxplot_potential.png")

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
       filename = "lat-le_points_potential.png")

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
       filename = "lat-le_boxplot_potential.png")

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
       filename = "lat-o_points_potential.png")

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
       filename = "lat-o_boxplot_potential.png")

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
       filename = "ele-o_points_potential.png")

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
       filename = "ele_o_boxplot_potential.png")

## so interesting!




### here is where things get messy





## model brainstorm:
## first want to ask if range shift lags are larger in species whose dispersal potential is less than the speed of climate change?


## see how to avoid singularity 
# key_lat = lat %>%
#   group_by(Param, class, expect_tracking) %>%
#   tally() %>%
#   ungroup() %>%
#   filter(n > 1) %>% ## get rid of classes with only 1 species per parameter per tracking group 
#   group_by(class, Param) %>%
#   filter(length(expect_tracking) == 2) %>%  ## get rid of classes that don't have species represented in each tracking group
#   select(-n)
# 
# lat = semi_join(lat, key_lat)

lat %>%
  group_by(Article_ID) %>%
  tally() %>%
  arrange(-n)
  


## filter out classes with fewer than 5 sp that have dispersal potential
## elevation:
dat.MaxDispPot.ele <- data.frame(ele[which(is.na(ele$MaxDispersalPotentialmY)==F),])

dsp <- unique(dat.MaxDispPot.ele[,c(which(colnames(dat.MaxDispPot.ele) == "class"),
                                     which(colnames(dat.MaxDispPot.ele) == "scientificName"))])
cl <- table(dsp$class)
nn <- names(which(cl >= 5))
dat.MaxDispPot.ele <- dat.MaxDispPot.ele[which(dat.MaxDispPot.ele$class %in% nn),]
dat.MaxDispPot.ele <-droplevels(dat.MaxDispPot.ele)
length(unique(dat.MaxDispPot.ele$class))

## latitude:
dat.MaxDispPot.lat <- data.frame(lat[which(is.na(lat$MaxDispersalPotentialKmY)==F),])

dsp <- unique(dat.MaxDispPot.lat[,c(which(colnames(dat.MaxDispPot.lat) == "class"),
                                    which(colnames(dat.MaxDispPot.lat) == "scientificName"))])
cl <- table(dsp$class)
nn <- names(which(cl >= 5))
dat.MaxDispPot.lat <- dat.MaxDispPot.lat[which(dat.MaxDispPot.lat$class %in% nn),]
dat.MaxDispPot.lat <-droplevels(dat.MaxDispPot.lat)
length(unique(dat.MaxDispPot.lat$class))

#scale variables
# elevation:
dat.MaxDispPot.ele$max_disp_pot <- scale(dat.MaxDispPot.ele$MaxDispersalPotentialmY)
dat.MaxDispPot.ele$n <- scale(dat.MaxDispPot.ele$N)
dat.MaxDispPot.ele$id.area <- scale(dat.MaxDispPot.ele$ID.area)
dat.MaxDispPot.ele$start <- scale(dat.MaxDispPot.ele$START)
dat.MaxDispPot.ele$dur <- scale(dat.MaxDispPot.ele$DUR)
dat.MaxDispPot.ele$shift <- scale(dat.MaxDispPot.ele$SHIFT)
dat.MaxDispPot.ele$lags_scaled <- scale(dat.MaxDispPot.ele$lags)
dat.MaxDispPot.ele$v.ele.mean_scaled <- scale(dat.MaxDispPot.ele$v.ele.mean)

# latitude:
dat.MaxDispPot.lat$max_disp_pot <- scale(dat.MaxDispPot.lat$MaxDispersalPotentialKmY)
dat.MaxDispPot.lat$n <- scale(dat.MaxDispPot.lat$N)
dat.MaxDispPot.lat$id.area <- scale(dat.MaxDispPot.lat$ID.area)
dat.MaxDispPot.lat$start <- scale(dat.MaxDispPot.lat$START)
dat.MaxDispPot.lat$dur <- scale(dat.MaxDispPot.lat$DUR)
dat.MaxDispPot.lat$shift <- scale(dat.MaxDispPot.lat$SHIFT)
dat.MaxDispPot.lat$lags_scaled <- scale(dat.MaxDispPot.lat$lags)
dat.MaxDispPot.lat$v.lat.mean_scaled <- scale(dat.MaxDispPot.lat$v.lat.mean)


## split by type of shift 
ele_te <- dat.MaxDispPot.ele[dat.MaxDispPot.ele$Param=="TE",]
ele_le <- dat.MaxDispPot.ele[dat.MaxDispPot.ele$Param=="LE",]
ele_o <- dat.MaxDispPot.ele[dat.MaxDispPot.ele$Param=="O",]
lat_te <- dat.MaxDispPot.lat[dat.MaxDispPot.lat$Param=="TE",]
lat_le <- dat.MaxDispPot.lat[dat.MaxDispPot.lat$Param=="LE",]
lat_o <- dat.MaxDispPot.lat[dat.MaxDispPot.lat$Param=="O",]

## first, model including only class as a random effect:
## latitude leading edge 
mod_lat_le = lmer(lags ~ expect_tracking_pot + (1|class), 
                  data = lat_le)

summary(mod_lat_le)
# fe = fixef(mod_lat_le)
# ci = confint(mod_lat_le, level = 0.95)
# 
# df = data.frame(expect_tracking = c("No", "Yes"),
#                 pred_lag = c(fe[1], (fe[1] + fe[2])), lower_ci = ci[3,1],
#                 upper_ci = c(ci[3,2], ci[4,2] + fe[1] )
# ## plot
# df %>%
#   ggplot(aes(x = expect_tracking, y = pred_lag)) + geom_point() +
#   geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci))

sjPlot::plot_model(mod_lat_le, type = "pred")


## elevation leading edge
mod_ele_le = lmer(lags ~ expect_tracking_pot + (1|class), 
                  data = ele_le)

summary(mod_ele_le)
sjPlot::plot_model(mod_ele_le, type = "pred")

## latitude optimum
mod_lat_o = lmer(lags ~ expect_tracking_pot + (1|class),
                 data = lat_o)

summary(mod_lat_o)
sjPlot::plot_model(mod_lat_o, type = "pred")

## elevation optimum
mod_ele_o = lmer(lags ~ expect_tracking_pot + (1|class), 
                 data = ele_o)

summary(mod_ele_o)
sjPlot::plot_model(mod_ele_o, type = "pred")




## next, add some important methodological covariates:
## Sampling = whether study was based on two time points or multiple
## Article_ID = control for study

## latitude leading edge 
mod_lat_le = lmer(lags ~ expect_tracking_pot + (1|Article_ID) + (1|class) + (1|Sampling), 
                  data = lat_le)

summary(mod_lat_le)
sjPlot::plot_model(mod_lat_le, type = "pred")


## elevation leading edge
mod_ele_le = lmer(lags ~ expect_tracking_pot + (1|Article_ID) + (1|class) + (1|Sampling), 
                  data = ele_le)

summary(mod_ele_le)
sjPlot::plot_model(mod_ele_le, type = "pred")

## latitude optimum
mod_lat_o = lmer(lags ~ expect_tracking_pot + (1|Article_ID) + (1|class) + (1|Sampling),
                 data = lat_o)

summary(mod_lat_o)
sjPlot::plot_model(mod_lat_o, type = "pred")

## elevation optimum
mod_ele_o = lmer(lags ~ expect_tracking_pot + (1|Article_ID) + (1|class) + (1|Sampling), 
                 data = ele_o)

summary(mod_ele_o)
sjPlot::plot_model(mod_ele_o, type = "pred")

## article id and sampling explain a lot of variation
## but even after accounting for them, species in which we expect climate tracking have a higher intercept (more positive lags)

## plot model predicted average with confidence intervals to see if confidence intervals overlap
new_data <- data.frame(expand_grid(expect_tracking_pot = c("Yes", "No"),
                                   Article_ID = unique(ele_le$Article_ID),
                                   class = unique(ele_le$class),
                                   Sampling = unique(ele_le$Sampling)))

pred <- predict(mod_ele_le)



## feedback:
## impute lifespan and age at maturity (try in a small subset (birds) and see what makes best imputation)
## see if v2 has velocity?
## see whether relationships do not exist at trailing edge
## add trailing edge + negative velocity shifts
## ask Sam for dispersal data
## see if we can add mean/mode other dispersal metrics
## try using taxonomic group instead of class (group plants together)
## co linearity between class and max dispersal distance?
## compare species with both elev and lat shifts to see if dispersal potential keeps up with climate velocity better across elev
## elev lags are in m/y!!!! convert


## modelling options:
## create model that predicts range shift values for each observation you need after accounting for methodological variables
## use these as your dependent variable
## easier to control for methodology when you can use ALL of the points 
## avoids signularities 
## OR
## include the methodological factors in your model 

## look within articles and see if the expected relationship is there 

## jonathon wants to use all dispersal data to look at the evolution of dispersal distance/potential

## Brown methodological affects

## high climate velicity and low range shift = high lag
## see 



## does range shift lag increase with dispersal potential in species we expect to have lags?
## latitude leading edge 
laggers_lat_le = filter(lat_le, expect_tracking_pot == "No")

mod_laggers_lat_le = lmer(lags_scaled ~ max_disp_pot + (1|class), 
                         data = laggers_lat_le)

summary(mod_laggers_lat_o)

laggers_lat_o = filter(lat_o, expect_tracking_pot == "No")

mod_laggers_lat_o = lmer(lags_scaled ~ max_disp_pot + (1|Article_ID) + (1|class) + (1|Sampling), 
                  data = laggers_lat_o)

summary(mod_laggers_lat_o)

laggers_ele_le = filter(ele_le, expect_tracking_pot == "No")

mod_laggers_ele_le = lmer(lags_scaled ~ max_disp_pot + (1|Article_ID) + (1|class) + (1|Sampling), 
                          data = laggers_ele_le)

summary(mod_laggers_ele_le)
## no 


#----------------------
# scrap code 

#----------------------
##select only classes with >5 species
## elevation
dat.MaxDispDist.ele <- data.frame(ele[which(is.na(ele$MaxDispersalDistanceKm)==F),])

dsp <- unique(dat.MaxDispDist.ele[,c(which(colnames(dat.MaxDispDist.ele) == "class"),
                                     which(colnames(dat.MaxDispDist.ele) == "scientificName"))])
cl <- table(dsp$class)
nn <- names(which(cl >10))
dat.MaxDispDist.ele <- dat.MaxDispDist.ele[which(dat.MaxDispDist.ele$class %in% nn),]
dat.MaxDispDist.ele <-droplevels(dat.MaxDispDist.ele)
length(unique(dat.MaxDispDist.ele$class))

#scale variables
dat.MaxDispDist.ele$max_disp_dist <- scale(dat.MaxDispDist.ele$MaxDispersalDistanceKm)
dat.MaxDispDist.ele$n <- scale(dat.MaxDispDist.ele$N)
dat.MaxDispDist.ele$id.area <- scale(dat.MaxDispDist.ele$ID.area)
dat.MaxDispDist.ele$start <- scale(dat.MaxDispDist.ele$START)
dat.MaxDispDist.ele$dur <- scale(dat.MaxDispDist.ele$DUR)
dat.MaxDispDist.ele$shift <- scale(dat.MaxDispDist.ele$SHIFT)
dat.MaxDispDist.ele$v.ele.mean_scaled <- scale(dat.MaxDispDist.ele$v.ele.mean)

## latitude 
dat.MaxDispDist.lat <- data.frame(lat[which(is.na(lat$MaxDispersalDistanceKm)==F),])

dsp <- unique(dat.MaxDispDist.lat[,c(which(colnames(dat.MaxDispDist.lat) == "class"),
                                     which(colnames(dat.MaxDispDist.lat) == "scientificName"))])
cl <- table(dsp$class)
nn <- names(which(cl >10))
dat.MaxDispDist.lat <- dat.MaxDispDist.lat[which(dat.MaxDispDist.lat$class %in% nn),]
dat.MaxDispDist.lat <-droplevels(dat.MaxDispDist.lat)
length(unique(dat.MaxDispDist.lat$class))

#scale variables
dat.MaxDispDist.lat$max_disp_dist <- scale(dat.MaxDispDist.lat$MaxDispersalDistanceKm)
dat.MaxDispDist.lat$n <- scale(dat.MaxDispDist.lat$N)
dat.MaxDispDist.lat$id.area <- scale(dat.MaxDispDist.lat$ID.area)
dat.MaxDispDist.lat$start <- scale(dat.MaxDispDist.lat$START)
dat.MaxDispDist.lat$dur <- scale(dat.MaxDispDist.lat$DUR)
dat.MaxDispDist.lat$shift <- scale(dat.MaxDispDist.lat$SHIFT)
dat.MaxDispDist.lat$v.lat.mean_scaled <- scale(dat.MaxDispDist.lat$v.lat.mean)


#----------------------
#prepare a dataset for each range shift parameter
ele_te <- dat.MaxDispDist.ele[dat.MaxDispDist.ele$Param=="TE",]
ele_le <- dat.MaxDispDist.ele[dat.MaxDispDist.ele$Param=="LE",]
ele_o <- dat.MaxDispDist.ele[dat.MaxDispDist.ele$Param=="O",]
lat_te <- dat.MaxDispDist.lat[dat.MaxDispDist.lat$Param=="TE",]
lat_le <- dat.MaxDispDist.lat[dat.MaxDispDist.lat$Param=="LE",]
lat_o <- dat.MaxDispDist.lat[dat.MaxDispDist.lat$Param=="O",]

## simple simple model to warm up:
lm(shift ~ max_disp_dist, data = lat_le)
lm(shift ~ max_disp_dist, data = ele_le)

## average across study for two time point data ?
## think about negative two time point studies 
## think about adding study duration interaction with dispersal 
## what was climate velocity per study duration?
## facet by low/med/heigh climate velocity 





m_lat = lmer(shift ~ v.lat.mean_scaled*max_disp_dist + (1|Article_ID) + (max_disp_dist|class), 
             data = lat_le)
summary(m_lat)

m_lat = lmer(shift ~ v.lat.mean_scaled*max_disp_dist + dur + id.area + n + (1|Article_ID) + (max_disp_dist|class), 
             data = lat_le)
summary(m_lat)





## things to ask: 
## is dispersal scale more important at leading edge or trailing edge?
## how much variation does dispersal scale explain at leading edge?



### thoughts:
## reproductive/dispersal frequency also matters!!!! need to know this too

## we have multiple observations of the same species, and the traits for each species are the same - include as a random effect?
## (Lise did not)
## count species within same article 
lat_le %>% 
  group_by(Article_ID, scientificName) %>%
  arrange(-n) %>%
  tally() %>% View




#----------------------
## model building notes from Lise:
# shifts ~ velocity x bodySize + StudyArea + StudyDuration + dataType + 
#   Nsamplings + dataProcessing + ( velocity x bodySize | Class) + (1|ArticleID)
#  To avoid singular fit (and shrinkage), I only selected taxonomic groups with > 10 or 15 species


## StudyArea = ID.area 
## StudyDuration = DUR
## dataType = Type
## Nsamplings = N
## dataProcessing = Uncertainty_Distribution
## Type = LE, TE, O
