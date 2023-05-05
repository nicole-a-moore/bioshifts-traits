## trying Lise & Jenn's new way of dealing with methodological variation
library(tidyverse)
library(lme4)
library(car)

## the idea: 
## remove variation that is due to differences in methodology between studies
## 1. calculate the average shift per study 
## 2. fit a model to the average shift per study as a function of methodological variables 
## (make sure none covary with each other or with variables of interest in our final model)
## 3. choose levels of each methodological factor and predict what mean range shift for each study would be if methods were consistent across studies 
## 4. calculate corrected shift by adding study-level predicted range shift to each range shift observation within studies 
## 5. model corrected shift as a function of variables of interest (for me: dispersal potential, climate velocity)

#############################
####   data preparation  ####
#############################
rs_data = read.table("data-raw/bioshifts-download/Lenoir_et_al/Analysis/Table_S1.csv",sep=";",h=T,dec=".",
                     stringsAsFactors = FALSE) 

## get start and end / duration data from other version of bioshifts 
other_rs_data <- read.delim('data-raw/bioshiftsv1/Shifts2018_checkedtaxo.txt')

other_rs_data <- select(other_rs_data, c(ID, START, END, DUR)) %>%
  distinct()

length(unique(other_rs_data$ID)) # 325

## if there are multiple durations per study, choose shortest 
other_rs_data = other_rs_data %>%
  group_by(ID) %>%
  mutate(DUR = min(DUR), 
         START = min(START),
         END = min(END)) %>%
  ungroup() %>%
  distinct() %>%
  rename("Source" = ID)

#it's the data used in Lenoir et al. 2020; so it's a subset of Bioshift v1 available here: https://figshare.com/articles/dataset/BioShifts_a_global_geodatabase_of_climate-induced_species_redistribution_over_land_and_sea/7413365?file=22057815

#Kingdom=	Taxonomic kingdom of the species for which the range shift has been estimated (according to NCBI; September 2019) 
#Phylum=	Taxonomic phylum of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Order=	Taxonomic order of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Class=	Taxonomic class of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Family=	Taxonomic family of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Genus=	Taxonomic genus of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Species=	Validated species name of the species for which the range shift has been estimated (according to NCBI; September 2019)
#Hemisphere=	Hemisphere of the study: North or South
#Ecosystem=	Ecosystem of the study: Marine or terrestrial (note that freshwater species have been classified as terrestrial)
#Gradient=	Spatial gradient along which the species range shift has been estimated: Elevation or Latitude
#Position=	"Range shift parameter for which the range shift has been estimated: Leading edge, Trailing edge, Centroid"
#ShiftR=	Range shift estimate standardized by the length of the study period over which the range shift has been estimated
#Unit=	Unit of the range shift estimate: m/year or km/year
#EleVeloT= elevational temperature velocity (m/year)
#LatVeloT= latitudinal temperature velocity (km/year)
#HFI = human footprint index (see Lenoir et al. 2020)
#BaselineT = baseline temperature (°C)
#LifeForm = lifeform of the species (endo, ecto, crypto, plant)
#Start = year of the beginning of the range shift monitoring
#Area = area of the polygon of the study (km2)
#Ntaxa = number of taxa monitored in the study
#PrAb=	Type of the data used to estimate the range shift: Occurrence (OCCUR) or Abundance (ABUND)
#Sampling=	"Characteristics of the time periods over which the range shift has been estimated: Continuous (CONT; yearly data or less), Irregular (IRR; multiple periods irregularly distributed), Multiple (MULT; multiple periods regularly distributed), Two (TWO; two periods)"
#Grain=	"Spatial grain of the data used to estimate the range shift: fine (FINE; spatial resolution lower than 10 km), coarse (COARSE; spatial resolution greater than 100 km), or medium (MEDIUM; intermediate situations)"
#Signif=	Whether the significance of the range shift estimate has been tested in the original publication (note that we do not report the result of such test in the database): YES or NO
#Quality=	"The quality of the approach used to estimate the range shift: LOW (no pre-processing of the raw data), BALANCED (data cleaning or resampling procedures were carried out to quantify the range shift on a balanced dataset), MODELED (range shifts were quantified species ranges modelled using  species distribution models calibrated independently for each time period); RESURVEYED (range shifts were quantified from paired designs such as permanent plots)"
#Start=	The first year of the temporal period over which the range shift has been estimated
#Reference=	The reference to the original publication
#Source=	Unique identifier of the spatial polygons delinating the study areas (provided in the geodatabase)

rs_data = rs_data[order(rs_data$n_cl),]

## group irregular sampling with multiple sampling  
rs_data$Sampling = ifelse(rs_data$Sampling %in% c("IRR","MULT"),"MULT", rs_data$Sampling)

## group by study ID and class and calculate a mean range shift value 
rs_study_level <- rs_data %>%
  group_by(Gradient, Source, Class) %>%
  mutate(study_level_shift = mean(ShiftR)) %>%
  select(-ShiftR, -Species, -Genus, -Kingdom, -Phylum, -Order, -Family, -Species,
         -n_cl, -n_sp, -Position) %>% ## make it such that there is one row per study
  distinct() %>%
  ungroup()

length(unique(paste(rs_data$Source, rs_data$Gradient, rs_data$Class, sep = "_"))) # 621 studies should remain

unique(rs_study_level$Source)[which(!unique(rs_data$Source) %in% unique(other_rs_data$Source))] ## make sure all are there

## add columns for start/end and duration 
rs_study_level <- left_join(rs_study_level, other_rs_data)

## log transform area and duration
rs_study_level$LogArea = log(rs_study_level$Area)
rs_study_level$LogDuration = log(rs_study_level$DUR)

## split data into latitudinal and elevational observations
rs_lat <- filter(rs_study_level, rs_study_level$Gradient == "Latitudinal")
rs_ele <- filter(rs_study_level, rs_study_level$Gradient == "Elevation")

## get complete cases of variables we want to include 
rs_lat <- select(rs_lat, Source,
                 study_level_shift, LogArea, Grain, Quality, Sampling, PrAb, Signif, LogDuration, LatVeloT,
                 Start, Start, Class)
# rs_lat <- select(rs_lat,  Source,
#                  study_level_shift, LogArea, Grain, Quality, Sampling, PrAb, Signif, LogDuration)
rs_lat <- rs_lat[complete.cases(rs_lat),]

rs_ele <- select(rs_ele, Source,
                 study_level_shift, LogArea, Grain, Quality, Sampling, PrAb, Signif, LogDuration, EleVeloT, 
                 Class, Start)
# rs_ele <- select(rs_ele, Source,
#                  study_level_shift, LogArea, Grain, Quality, Sampling, PrAb, Signif, LogDuration)
rs_ele <- rs_ele[complete.cases(rs_ele),]

rs_ele$Gradient = "Elevation"
rs_lat$Gradient = "Latitude"
combined <- select(rs_lat, -LatVeloT) %>%
  rbind(., select(rs_ele, -EleVeloT))

#####################################
####        model fitting        ####
#####################################
## fit one model per gradient (latitudinal, elevational) 

############################
####      latitude      ####
############################
## look at data 
plot1 = rs_lat %>%
  ggplot(., aes(x = Sampling, fill = Class)) +
  geom_bar() +
  #geom_point(position = position_jitter()) +
  labs(y = "Mean study-level shift (km/y)") 

legend = cowplot::get_legend(plot1)

ggsave(legend, filename = "lat_class_legend.png", path = "figures/methodo", width = 6, height = 6, 
       device = "png")

plot1 <- plot1 + theme(legend.position = "none")

ggsave(plot1, filename = "lat_class_sampling.png", path = "figures/methodo", width = 3, height = 3, 
       device = "png")


combined %>%
  ggplot(., aes(x = Sampling, fill = Class)) +
  geom_bar() +
  facet_wrap(~Gradient) +
  #geom_point(position = position_jitter()) +
  labs(y = "Mean study-level shift (km/y)")  +
  theme(legend.position = "none")

ggsave(filename = "class_sampling.png", path = "figures/methodo", width = 5, height = 3, 
       device = "png")

combined %>%
  group_by(Gradient, Class, Sampling) %>%
  tally() %>%
  group_by(Gradient, Class) %>%
  tally() %>%
  group_by(Gradient, n) %>%
  tally()

combined %>%
  ggplot(., aes(x = PrAb, fill = Class)) +
  geom_bar() +
  facet_wrap(~Gradient) +
  #geom_point(position = position_jitter()) +
  labs(y = "Number of study/class combinations")  +
  theme(legend.position = "none")

ggsave(filename = "class_prab.png", path = "figures/methodo", width = 5, height = 3, 
       device = "png")

combined %>%
  group_by(Gradient, Class, PrAb) %>%
  tally() %>%
  group_by(Gradient, Class) %>%
  tally() %>%
  group_by(Gradient, n) %>%
  tally()


combined %>%
  ggplot(., aes(x = Grain, fill = Class)) +
  geom_bar() +
  facet_wrap(~Gradient) +
  #geom_point(position = position_jitter()) +
  labs(y = "Number of study/class combinations")  +
  theme(legend.position = "none")

ggsave(filename = "class_grain.png", path = "figures/methodo", width = 5, height = 3, 
       device = "png")

combined %>%
  group_by(Gradient, Class, Grain) %>%
  tally() %>%
  group_by(Gradient, Class) %>%
  tally() %>%
  group_by(Gradient, n) %>%
  tally()





## how many classes have been studied using all types of sampling at least once?  
length(unique(rs_lat$Class))

rs_lat %>%
  group_by(Class, Sampling) %>%
  tally() %>%
  group_by(Class) %>%
  tally() %>%
  group_by(n) %>%
  tally()






rs_lat %>%
  ggplot(., aes(x = Grain, y = study_level_shift)) +
  geom_violin() +
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "lat_classMSLS_x_grain.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")

rs_lat %>%
  ggplot(., aes(x = Quality, y = study_level_shift)) +
  geom_violin()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "lat_classMSLS_x_quality.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")

rs_lat %>%
  ggplot(., aes(x = PrAb, y = study_level_shift)) +
  geom_violin()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "lat_classMSLS_x_prab.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")

rs_lat %>%
  ggplot(., aes(x = Signif, y = study_level_shift)) +
  geom_violin()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "lat_classMSLS_x_signif.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")

rs_lat %>%
  ggplot(., aes(x = LogArea, y = study_level_shift)) +
  geom_point()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "lat_classMSLS_x_logarea.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")

rs_lat %>%
  ggplot(., aes(x = LogDuration, y = study_level_shift)) +
  geom_point()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "lat_classMSLS_x_logduration.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")


rs_lat %>%
  ggplot(., aes(x = Start, y = study_level_shift)) +
  geom_point()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "lat_classMSLS_x_start.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")


## now elevation
rs_ele %>%
  ggplot(., aes(x = Grain, y = study_level_shift)) +
  geom_violin() +
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "ele_classMSLS_x_grain.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")

rs_ele %>%
  ggplot(., aes(x = Quality, y = study_level_shift)) +
  geom_violin()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "ele_classMSLS_x_quality.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")

rs_ele %>%
  ggplot(., aes(x = PrAb, y = study_level_shift)) +
  geom_violin()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "ele_classMSLS_x_prab.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")

rs_ele %>%
  ggplot(., aes(x = Signif, y = study_level_shift)) +
  geom_violin()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "ele_classMSLS_x_signif.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")

rs_ele %>%
  ggplot(., aes(x = LogArea, y = study_level_shift)) +
  geom_point()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "ele_classMSLS_x_logarea.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")

rs_ele %>%
  ggplot(., aes(x = LogDuration, y = study_level_shift)) +
  geom_point()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "ele_classMSLS_x_logduration.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")


rs_ele %>%
  ggplot(., aes(x = Start, y = study_level_shift)) +
  geom_point()+
  geom_point(position = position_jitter(), aes(colour = Class))+
  labs(y = "Mean study-level shift (km/y)") +
  theme(legend.position = "none")

ggsave(filename = "ele_classMSLS_x_start.png", path = "figures/methodo", width = 4, height = 4, 
       device = "png")




## see if things vary with climate velocity:
rs_lat %>%
  ggplot(., aes(x = Sampling, y = LatVeloT)) +
  geom_violin()

rs_lat %>%
  ggplot(., aes(x = Grain, y = LatVeloT)) +
  geom_violin()

rs_lat %>%
  ggplot(., aes(x = Quality, y = LatVeloT)) +
  geom_violin()

rs_lat %>%
  ggplot(., aes(x = PrAb, y = LatVeloT)) +
  geom_violin()

rs_lat %>%
  ggplot(., aes(x = Signif, y = LatVeloT)) +
  geom_violin()

rs_lat %>%
  ggplot(., aes(x = LogArea, y = LatVeloT)) +
  geom_point() 

rs_lat %>%
  ggplot(., aes(x = LogDuration, y = LatVeloT)) +
  geom_point() 

rs_lat %>%
  ggplot(., aes(x = Start, y = LatVeloT)) +
  geom_point() 

## order factors in a logical way
rs_lat$Grain <- factor(rs_lat$Grain, levels = c("FINE", "MEDIUM", "COARSE"))
## change to continuous
## keep only area, grain, duration

## fit model for latitude 
mod_lat <- lm(study_level_shift ~ LogArea + Grain + Quality + Sampling + PrAb + Signif + LogDuration,
              data = rs_lat) # velocity not colinear with methodological variables
# mod_lat <- lm(study_level_shift ~ LogArea + Grain + LogDuration, 
#    data = rs_lat)
summary(mod_lat)

## check for colinearity using VIFs:
vif(mod_lat)

## get r2 
summary(mod_lat)$r.squared
# 0.2219915

## plot model predictions 
## continuous variables
new_data <- data.frame(expand_grid(LogArea = seq(min(rs_lat$LogArea), max(rs_lat$LogArea), by = 0.01),
                                   Grain = unique(rs_lat$Grain)[1],
                                   Quality = unique(rs_lat$Quality)[1],
                                   Sampling = unique(rs_lat$Sampling)[1],
                                   PrAb = unique(rs_lat$PrAb)[1],
                                   Signif = unique(rs_lat$Signif)[1],
                                   LogDuration = seq(min(rs_lat$LogDuration), max(rs_lat$LogDuration), by = 0.01)))

pred_lat <- predict(mod_lat, new_data, level = 0, se.fit = T, re.form = NA)

pred_lat <- new_data %>%
  mutate(pred_val = pred_lat$fit) %>%
  mutate(pred_val_SE = pred_lat$se.fit)

## AREA
pred_lat %>%
  filter(LogDuration == first(.$LogDuration)) %>%
  ggplot(., aes(x = LogArea, y = pred_val)) + 
  geom_point() +
  theme_bw() +
  labs(y = "Predicted study-level range shift")
## shift increases with study area

## DURATION
pred_lat %>%
  filter(LogArea == first(.$LogArea)) %>%
  ggplot(., aes(x = LogDuration, y = pred_val)) + 
  geom_point() +
  theme_bw() +
  labs(y = "Predicted study-level range shift", x = "Study duration") +
  geom_point(data = rs_lat, aes(x  = LogDuration, y = study_level_shift), inherit.aes = FALSE)
## shift decreases with study duration

## categorical variables 
new_data2 <- data.frame(expand_grid(LogArea = median(rs_lat$LogArea),
                                   Grain = unique(rs_lat$Grain),
                                   Quality = unique(rs_lat$Quality),
                                   Sampling = unique(rs_lat$Sampling),
                                   PrAb = unique(rs_lat$PrAb),
                                   Signif = unique(rs_lat$Signif),
                                   LogDuration = median(rs_lat$LogDuration)))

pred_lat2 <- predict(mod_lat, new_data2, level = 0, se.fit = T, re.form = NA)

pred_lat2 <- new_data2 %>%
  mutate(pred_val = pred_lat2$fit) %>%
  mutate(pred_val_SE = pred_lat2$se.fit)

## GRAIN
pred_lat2 %>%
  ggplot(., aes(x = Grain, y = pred_val)) + 
  geom_boxplot() +
  theme_bw() +
  labs(y = "Predicted study-level range shift")
## course grain shifts are smaller than fine/medium

## QUALITY
pred_lat2 %>%
  ggplot(., aes(x = Quality, y = pred_val)) + 
  geom_boxplot() +
  theme_bw() +
  labs(y = "Predicted study-level range shift")
## low-quality shifts are larger, modeled shifts are smallest

## SIGNIF
pred_lat2 %>%
  ggplot(., aes(x = Signif, y = pred_val)) + 
  geom_boxplot() +
  theme_bw() +
  labs(y = "Predicted study-level range shift")
## shifts where authors test significance are smaller

## SAMPLING
pred_lat2 %>%
  ggplot(., aes(x = Sampling, y = pred_val)) + 
  geom_boxplot() +
  theme_bw() +
  labs(y = "Predicted study-level range shift")
## two-time point shifts are smallest 

## PRAB
pred_lat2 %>%
  ggplot(., aes(x = PrAb, y = pred_val)) + 
  geom_boxplot() +
  theme_bw() +
  labs(y = "Predicted study-level range shift",
       x = "Observation type")
## shifts measuring abundance are smaller 


#############################
####      elevation      ####
#############################
## look at data 
rs_ele %>%
  ggplot(., aes(x = Sampling, y = study_level_shift)) +
  geom_violin()

rs_ele %>%
  ggplot(., aes(x = Grain, y = study_level_shift)) +
  geom_violin()

rs_ele %>%
  ggplot(., aes(x = Quality, y = study_level_shift)) +
  geom_violin()

rs_ele %>%
  ggplot(., aes(x = PrAb, y = study_level_shift)) +
  geom_violin()

rs_ele %>%
  ggplot(., aes(x = Signif, y = study_level_shift)) +
  geom_violin()

rs_ele %>%
  ggplot(., aes(x = LogArea, y = study_level_shift)) +
  geom_point()

rs_ele %>%
  ggplot(., aes(x = LogDuration, y = study_level_shift)) +
  geom_point() 

## order factors in a logical way
rs_ele$Grain <- factor(rs_ele$Grain, levels = c("FINE", "MEDIUM", "COARSE"))


## fit model for latitude 
mod_ele <- lm(study_level_shift ~ LogArea + Grain + Quality + Sampling + PrAb + Signif + LogDuration,
              data = rs_ele) # velocity not colinear with methodological variables
# mod_ele <- lm(study_level_shift ~ LogArea + Grain + LogDuration, 
#               data = rs_ele)
summary(mod_ele)

## check for colinearity using VIFs:
vif(mod_ele)

## get r2 
summary(mod_ele)$r.squared
# 0.1327721

## plot model predictions 
## continuous variables
new_data <- data.frame(expand_grid(LogArea = seq(min(rs_ele$LogArea), max(rs_ele$LogArea), by = 0.01),
                                   Grain = unique(rs_ele$Grain)[1],
                                   Quality = unique(rs_ele$Quality)[1],
                                   Sampling = unique(rs_ele$Sampling)[1],
                                   PrAb = unique(rs_ele$PrAb)[1],
                                   Signif = unique(rs_ele$Signif)[1],
                                   LogDuration = seq(min(rs_ele$LogDuration), max(rs_ele$LogDuration), by = 0.01)))

pred_ele <- predict(mod_ele, new_data, level = 0, se.fit = T, re.form = NA)

pred_ele <- new_data %>%
  mutate(pred_val = pred_ele$fit) %>%
  mutate(pred_val_SE = pred_ele$se.fit)

## AREA
pred_ele %>%
  filter(LogDuration == first(.$LogDuration)) %>%
  ggplot(., aes(x = LogArea, y = pred_val)) + 
  geom_point() +
  theme_bw() +
  labs(y = "Predicted study-level range shift")
## shift increases with study area

## DURATION
pred_ele %>%
  filter(LogArea == first(.$LogArea)) %>%
  ggplot(., aes(x = LogDuration, y = pred_val)) + 
  geom_point() +
  theme_bw() +
  labs(y = "Predicted study-level range shift", x = "Study duration")
## shift decreases with study duration

## categorical variables 
new_data2 <- data.frame(expand_grid(LogArea = median(rs_ele$LogArea),
                                    Grain = unique(rs_ele$Grain),
                                    Quality = unique(rs_ele$Quality),
                                    Sampling = unique(rs_ele$Sampling),
                                    PrAb = unique(rs_ele$PrAb),
                                    Signif = unique(rs_ele$Signif),
                                    LogDuration = median(rs_ele$LogDuration)))

pred_ele2 <- predict(mod_ele, new_data2, level = 0, se.fit = T, re.form = NA)

pred_ele2 <- new_data2 %>%
  mutate(pred_val = pred_ele2$fit) %>%
  mutate(pred_val_SE = pred_ele2$se.fit)

## GRAIN
pred_ele2 %>%
  ggplot(., aes(x = Grain, y = pred_val)) + 
  geom_boxplot() +
  theme_bw() +
  labs(y = "Predicted study-level range shift")
## course grain shifts are larger than fine/medium

## QUALITY
pred_ele2 %>%
  ggplot(., aes(x = Quality, y = pred_val)) + 
  geom_boxplot() +
  theme_bw() +
  labs(y = "Predicted study-level range shift")
## low-quality shifts are larger, modeled shifts are smallest

## SIGNIF
pred_ele2 %>%
  ggplot(., aes(x = Signif, y = pred_val)) + 
  geom_boxplot() +
  theme_bw() +
  labs(y = "Predicted study-level range shift")
## shifts where authors test significance are smaller

## SAMPLING
pred_ele2 %>%
  ggplot(., aes(x = Sampling, y = pred_val)) + 
  geom_boxplot() +
  theme_bw() +
  labs(y = "Predicted study-level range shift")
## two-time point shifts are smallest 

## PRAB
pred_ele2 %>%
  ggplot(., aes(x = PrAb, y = pred_val)) + 
  geom_boxplot() +
  theme_bw() +
  labs(y = "Predicted study-level range shift",
       x = "Observation type")
## about the same 


## make sure duration isn't confounded by taxonomy
# rs_study_level <- rs_data %>%
#   group_by(Gradient, Source) %>%
#   mutate(study_level_shift = mean(ShiftR)) %>%
#   select(-ShiftR, -Species, -Genus, -Order,-Family, -Species,
#          -n_cl, -n_sp, -Position) %>% ## make it such that there is one row per study
#   distinct() %>%
#   ungroup()
# 
# rs_study_level <- left_join(rs_study_level, other_rs_data)
# 
# rs_study_level %>%
#   ggplot(., aes(x = DUR, y = study_level_shift, colour = Kingdom)) +
#   geom_point() +
#   scale_x_log10() +
#   facet_wrap(~Gradient)


############################################
####    applying correction to shifts   ####
############################################
## calculate difference between observed study-level mean and predicted study-level mean, holding methods constant across all studies

## DATA FRAME STRUCTURE:
## Source (study ID)
## Observed mean study-level shift 
## Predicted study-level shift, if all variables held to mode value 

# function to find mode of each predictor
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## get mode of each categorical predictor across studies
rs_lat$Grain <- as.character(rs_lat$Grain)
rs_ele$Grain <- as.character(rs_ele$Grain)
modes_lat <- sapply(select(rs_lat, c("Grain", "Quality", "Sampling", 
                                     "PrAb", "Signif")), FUN = getmode)
modes_ele <- sapply(select(rs_ele, c("Grain", "Quality", "Sampling", 
                                     "PrAb", "Signif")), FUN = getmode)

## get median of each continuous predictor across studies 
medians_lat <- sapply(select(rs_lat, c("LogArea", "LogDuration")), FUN = median)
medians_ele <- sapply(select(rs_ele, c("LogArea", "LogDuration")), FUN = median)

## predict mean study level range shift, holding all methodological variables at their median/mode value
new_data_lat <- cbind(data.frame(as.list(modes_lat)), data.frame(as.list(medians_lat)))
new_data_ele <- cbind(data.frame(as.list(modes_ele)), data.frame(as.list(medians_ele)))

corr_rs_lat <- predict(mod_lat, new_data_lat)
corr_rs_ele <- predict(mod_ele, new_data_ele)

corr_rs_lat <- select(rs_lat, Source, study_level_shift) %>%
  mutate(PredSLShift = corr_rs_lat, 
         Gradient = "Latitudinal")
corr_rs_ele <- select(rs_ele, Source, study_level_shift) %>%
  mutate(PredSLShift = corr_rs_ele,
         Gradient = "Elevation")

## bind ele and lat corr shift data 
corr_rs <- rbind(corr_rs_ele, corr_rs_lat) 

## calculate difference between observed study-level mean shift and predicted study-level mean shift 
corr_rs$SLDiff <- corr_rs$study_level_shift - corr_rs$PredSLShift

## now, calculate corrected range shift
## subtract study level difference from raw range shift observations within each study
rs_data <- left_join(rs_data, corr_rs) %>%
  distinct()

rs_data$CorrShift <- rs_data$ShiftR - rs_data$SLDiff



#######################
####    plotting   ####
#######################
hist(rs_data$CorrShift)

## now plot raw versus corrected shift 
o_vs_c <- rs_data %>%
  ggplot(., aes(x = ShiftR, y = CorrShift, colour = Source)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Observed range shift", y = "Corrected range shift") +
  geom_abline(intercept = 0, slope = 1)

ggsave(o_vs_c, path = "figures/corrected-shifts", filename = "observed-versus-corrected-points.png", 
       "png", width = 5, height = 4)

hists <- rs_data %>%
  gather(key = "shift_type", value = "shift_value", 
         c(ShiftR, CorrShift)) %>% 
  ggplot(., aes(x = shift_value, fill = shift_type)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(Gradient~Position) +
  scale_fill_discrete(labels = c("Corrected shift", "Observed shift")) +
  labs(x = "Range shift", y = "Frequency", fill = "")

ggsave(hists, path = "figures/corrected-shifts", filename = "observed-versus-corrected-histograms.png", 
       "png", width = 7, height = 4)

diffs <- rs_data %>%
  ggplot(., aes(x = SLDiff)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(Gradient~Position) +
  labs(x = "Study-level difference", y = "Frequency")

ggsave(diffs, path = "figures/corrected-shifts", filename = "observed-versus-corrected-diff-histograms.png", 
       "png", width = 7, height = 4)

## save data set 
write.csv(rs_data, "data-processed/corrected-bioshifts.csv", row.names = FALSE)



## plot example:
lat_corr <- filter(rs_data, Gradient == "Latitudinal")

lat_corr %>%
  ggplot(aes(x = log(Area), y = ShiftR)) +
  geom_point() +
  geom_point(data = rs_lat, aes(x = LogArea, y = study_level_shift), colour = "red")

lat_corr %>%
  ggplot(aes(x = log(Area), y = ShiftR)) +
  geom_point() +
  geom_point(data = rs_lat, aes(x = LogArea, y = study_level_shift), colour = "red") +
  geom_point(aes(y = CorrShift), colour = "blue")

