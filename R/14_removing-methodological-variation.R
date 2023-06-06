## trying new way of dealing with methodological variation
## BIOSHIFTS meeting April 27 
library(tidyverse)
library(nlme)
library(car)

theme_set(theme_bw())

## The idea: 

# 1. First, fit a class- and study-level model of methods (separate for lat and ele):
#   mean shift class within study ~ sampling + abun/occur + grain_size + modeled/raw + area, 1| Class
#   -(actual variables will vary by version of bioshifts)
#   -weight = num_species OR remove low-species number studies e.g. <5 OR weight by study-level variance (i.e. precision around mean)
#   -remove variables that are co-linear with climate velocity (e.g. start year)
#   -random effect of Class on the intercept because methods might differ across classes
#   -do not include parameter (upper, lower, o) because methods don’t differ across params, assume effect of method is the same across params
#   -only need to include methods thought to influence range shift direction (not variance)

# 2. Then, standardize to a given method level
#   -calculate the study-level residual from model
#   -choose standardized level for every method
#   -predict model at standardized levels

# 3. Add study-level residual to the model-prediction, and add this same adjustment to every data point within that study
#   -predicted value + study residual – study value


#############################
####   data preparation  ####
#############################
## read in data:
rs_data = read.table("data-raw/bioshifts-download/Lenoir_et_al/Analysis/Table_S1.csv",sep=";",h=T,dec=".",
                     stringsAsFactors = FALSE) 
#it's the data used in Lenoir et al. 2020; so it's a subset of Bioshift v1 available here: https://figshare.com/articles/dataset/BioShifts_a_global_geodatabase_of_climate-induced_species_redistribution_over_land_and_sea/7413365?file=22057815

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

rs_data = rs_data[order(rs_data$n_cl),]

## log transform area 
rs_data$LogArea = log(rs_data$Area)

## group irregular sampling with multiple sampling  
rs_data$Sampling = ifelse(rs_data$Sampling %in% c("IRR","MULT"),"MULT", rs_data$Sampling)

## order factors in a logical way
rs_data$Grain <- factor(rs_data$Grain, levels = c("FINE", "MEDIUM", "COARSE"))
## change to continuous?
## comes from log scale 

rs_data$ContinuousGrain <- as.numeric(as.character(ifelse(rs_data$Grain == "FINE", 
                                                                 1, 
                                                                 ifelse(rs_data$Grain == "MEDIUM",
                                                                        10,
                                                                        ifelse(rs_data$Grain == "COARSE",
                                                                               100,
                                                                               NA))))) 
unique(rs_data$ContinuousGrain)

## make new var for unique study id + class 
rs_data$ID = paste(rs_data$Source, rs_data$Gradient, rs_data$Class, sep = "_")

## get rid of study with insanely large area
rs_data <- filter(rs_data, LogArea < 9)

## split raw data by gradient
rs_data_lat <- filter(rs_data, Gradient == "Latitudinal")
rs_data_ele <- filter(rs_data, Gradient == "Elevation")
  

## group by study ID and class and calculate a mean range shift value and variance for each class within studies
rs_study_level <- rs_data %>%
  group_by(Gradient, Source, Class) %>%
  mutate(study_level_shift = mean(ShiftR),
         study_level_var = var(ShiftR)) %>%
  select(-ShiftR, -Species, -Genus, -Kingdom, -Phylum, -Order, -Family, -Species,
         -n_cl, -n_sp, -Position) %>% ## make it such that there is one row per study
  distinct() %>%
  ungroup()

## filter to studies with at least 5 species 
## justification - average shift of studies with few species might be biased because of species traits, or publication bias
## (might want to test the sensitivity to this)
rs_study_level <- filter(rs_study_level, Ntaxa >= 5)

length(unique(paste(rs_data$Source, rs_data$Gradient, rs_data$Class, sep = "_"))) # 621 studies should remain

rs_study_level <- select(rs_study_level, Gradient,
       Source, study_level_shift, study_level_var, LogArea, ContinuousGrain, Quality, PrAb, Class, Ntaxa,
       Sampling, ID)
rs_study_level <- rs_study_level[complete.cases(rs_study_level),]

## split data into latitudinal and elevational observations
rs_lat <- filter(rs_study_level, rs_study_level$Gradient == "Latitudinal")
rs_ele <- filter(rs_study_level, rs_study_level$Gradient == "Elevation")


#####################################
####        model fitting        ####
#####################################
## fit one model per gradient (latitudinal, elevational) 

############################
####      latitude      ####
############################
## fit model for latitude 
mod_lat <- lme(study_level_shift ~ LogArea + ContinuousGrain + Quality + PrAb,
                random = (~1|Class), 
               weights = ~I(1/study_level_var),
              data = rs_lat) 
summary(mod_lat)

## get residuals 
rs_lat$residuals <- resid(mod_lat)
hist(rs_lat$residuals)
plot(mod_lat)

## check for colinearity using VIFs:
vif(mod_lat)

###############################################
####      plot predictions - latitude      ####
###############################################
## plot model predictions 
new_data <- data.frame(expand_grid(LogArea = seq(min(rs_lat$LogArea), 
                                                 max(rs_lat$LogArea), by = 0.1),
                                   ContinuousGrain = c(1, 10, 100),
                                   PrAb = unique(rs_lat$PrAb),
                                   Quality = unique(rs_lat$Quality),
                                   Class = unique(rs_lat$Class)))

pred_lat <- predict(mod_lat, new_data, level = 0:1, se.fit = T, re.form = NA)

new_data$pred_val <- pred_lat$predict.fixed
new_data$pred_val_class <- pred_lat$predict.Class

# function to find the mode of each predictor
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## get mode of each predictor
mode_grain = getmode(rs_lat$ContinuousGrain)
mode_prab = getmode(rs_lat$PrAb)
mode_quality = getmode(rs_lat$Quality)
area = mean(rs_lat$LogArea)
mode_class = getmode(rs_lat$Class)

## plot predictions and raw data together:
## AREA
new_data %>%
  filter(PrAb == mode_prab, ## plot predictions for most common levels
         ContinuousGrain == mode_grain,
         Quality == mode_quality,
         Class == mode_class
         ) %>%
  ggplot(aes(x = LogArea, y = pred_val_class)) +
  geom_line(aes(group = Class)) +
  geom_point(data = rs_data_lat, aes(x = LogArea, y = ShiftR, colour = ID), 
             size = 0.25, alpha = 0.75) + ## small circles are raw shifts within studies
  geom_point(data = rs_lat, aes(x = LogArea, y = study_level_shift, fill = ID),
             shape = 21, colour = "black", alpha = 0.75) + ## large circles are mean study-level shift per class
  theme(legend.position = "none") +
  labs(y = "Mean range shift rate (class within study)")
## study level shift decreases with study area

## GRAIN and QUALITY
new_data %>%
  filter(LogArea == new_data$LogArea[which(new_data$LogArea - area == min(new_data$LogArea - area))],
         PrAb == mode_prab,
         Class == mode_class,
         Quality == mode_quality
         ) %>%
  ggplot(aes(x = ContinuousGrain, y = pred_val_class)) +
  geom_line()+
  theme(legend.position = "none") +
  geom_point(data = rs_data_lat, aes(x = ContinuousGrain, y = ShiftR, colour = ID), 
             size = 0.25, alpha = 0.75,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  geom_point(data = rs_lat, aes(x = ContinuousGrain, y = study_level_shift, fill = ID),
             shape = 21, colour = "black", alpha = 0.75,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  labs(y = "Mean range shift rate (class within study)")
## no obvious trend with grain size
## note: we may need to model as categorical variable?

## PRAB 
new_data %>%
  filter(LogArea <= 7.859033, ## choose large area
    #LogArea == new_data$LogArea[which(new_data$LogArea - area == min(new_data$LogArea - area))],
         ContinuousGrain == mode_grain,
         Class == mode_class,
         ) %>%
  ggplot(aes(x = PrAb, y = pred_val_class)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  geom_point(data = rs_data_lat, aes(x = PrAb, y = ShiftR, colour = ID), 
             size = 0.25, alpha = 0.75,
             position = position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0)) +
  geom_point(data = rs_lat, aes(x = PrAb, y = study_level_shift, fill = ID),
             shape = 21, colour = "black", alpha = 0.75,
             position = position_jitterdodge(jitter.width = 0.3,
                                             jitter.height = 0)) +
  labs(y = "Mean range shift rate (class within study)")
## note: predictions look funky here and don't match the data but I can't figure out why
## study level shift is slightly greater for abundance data 

#############################
####      elevation      ####
#############################
## fit model for elevation 
## get rid of study with zero variance
rs_ele <- filter(rs_ele, study_level_var != 0)
mod_ele <- lme(study_level_shift ~ LogArea + ContinuousGrain + Quality + PrAb,
               random = (~1|Class), 
               weights = ~I(1/study_level_var),
               data = rs_ele) 
summary(mod_ele)

## get residuals 
rs_ele$residuals <- resid(mod_ele)
hist(rs_ele$residuals)
plot(mod_ele)

## check for colinearity using VIFs:
vif(mod_ele)

###############################################
####      plot predictions - elevation      ####
###############################################
## plot model predictions 
new_data <- data.frame(expand_grid(LogArea = seq(min(rs_ele$LogArea), 
                                                 max(rs_ele$LogArea), by = 0.1),
                                   ContinuousGrain = c(1, 10, 100),
                                   PrAb = unique(rs_ele$PrAb),
                                   Quality = unique(rs_ele$Quality),
                                   Class = unique(rs_ele$Class)))

pred_ele <- predict(mod_ele, new_data, level = 0:1, se.fit = T, re.form = NA)

new_data$pred_val <- pred_ele$predict.fixed
new_data$pred_val_class <- pred_ele$predict.Class

## get mode of each predictor
mode_grain = getmode(rs_ele$ContinuousGrain)
mode_prab = getmode(rs_ele$PrAb)
mode_quality = getmode(rs_ele$Quality)
area = mean(rs_ele$LogArea)
mode_class = getmode(rs_ele$Class)

## plot predictions and raw data together:
## AREA
new_data %>%
  filter(PrAb == mode_prab, ## plot predictions for most common levels
         ContinuousGrain == mode_grain,
         Quality == mode_quality,
         Class == mode_class
  ) %>%
  ggplot(aes(x = LogArea, y = pred_val_class)) +
  geom_line(aes(group = Class)) +
  geom_point(data = rs_data_ele, aes(x = LogArea, y = ShiftR, colour = ID), 
             size = 0.25, alpha = 0.75) + ## small circles are raw shifts within studies
  geom_point(data = rs_ele, aes(x = LogArea, y = study_level_shift, fill = ID),
             shape = 21, colour = "black", alpha = 0.75) + ## large circles are mean study-level shift per class
  theme(legend.position = "none") +
  labs(y = "Mean range shift rate (class within study)")
## no obvious trend between study level shift and study area

## GRAIN and QUALITY
new_data %>%
  filter(LogArea == new_data$LogArea[which(new_data$LogArea - area == min(new_data$LogArea - area))],
         PrAb == mode_prab,
         Class == mode_class,
         Quality == mode_quality
  ) %>%
  ggplot(aes(x = ContinuousGrain, y = pred_val_class)) +
  geom_line()+
  theme(legend.position = "none") +
  geom_point(data = rs_data_ele, aes(x = ContinuousGrain, y = ShiftR, colour = ID), 
             size = 0.25, alpha = 0.75,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  geom_point(data = rs_ele, aes(x = ContinuousGrain, y = study_level_shift, fill = ID),
             shape = 21, colour = "black", alpha = 0.75,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  labs(y = "Mean range shift rate (class within study)")
## no obvious trend with grain size
## note: we may need to model as categorical variable?

## PRAB 
new_data %>%
  filter(LogArea <= 7.818, ## choose large area
         #LogArea == new_data$LogArea[which(new_data$LogArea - area == min(new_data$LogArea - area))],
         ContinuousGrain == mode_grain,
         Class == mode_class,
  ) %>%
  ggplot(aes(x = PrAb, y = pred_val_class)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  geom_point(data = rs_data_lat, aes(x = PrAb, y = ShiftR, colour = ID), 
             size = 0.25, alpha = 0.75,
             position = position_jitterdodge(jitter.width = 0.3, 
                                             jitter.height = 0)) +
  geom_point(data = rs_lat, aes(x = PrAb, y = study_level_shift, fill = ID),
             shape = 21, colour = "black", alpha = 0.75,
             position = position_jitterdodge(jitter.width = 0.3,
                                             jitter.height = 0)) +
  labs(y = "Mean range shift rate (class within study)")
## study level shift is slightly greater for abundance data 



############################################
####    applying correction to shifts   ####
############################################
## next step is to calculate the predicted study-level shift for each study, holding methods constant across all studies
## we want to leave variation due to class in, so we will ignore the random effect

## predict mean study level range shift per class, holding all methodological variables at their median/mode value
new_data_lat <- data.frame(expand_grid(LogArea = seq(min(rs_lat$LogArea), 
                                                                 max(rs_lat$LogArea), by = 0.1),
                                                   ContinuousGrain = c(1, 10, 100),
                                                   PrAb = unique(rs_lat$PrAb),
                                                   Quality = unique(rs_lat$Quality)))

corr_rs_lat <- predict(mod_lat, new_data_lat, level = 0, se.fit = T, re.form = NA)
new_data_lat$pred_val <- corr_rs_lat

new_data_ele <- data.frame(expand_grid(LogArea = seq(min(rs_ele$LogArea), 
                                                     max(rs_ele$LogArea), by = 0.1),
                                       ContinuousGrain = c(1, 10, 100),
                                       PrAb = unique(rs_ele$PrAb),
                                       Quality = unique(rs_ele$Quality)))

corr_rs_ele <- predict(mod_ele, new_data_ele, level = 0, se.fit = T, re.form = NA)
new_data_ele$pred_val <- corr_rs_ele

## filter to only predictions for the mode of each level + the smallest area value 
new_data_lat <- filter(new_data_lat, 
       ContinuousGrain == getmode(rs_lat$ContinuousGrain),
       PrAb == getmode(rs_lat$PrAb),
       LogArea == first(new_data_lat$LogArea),
       Quality == getmode(rs_lat$Quality))

new_data_ele <- filter(new_data_ele, 
                       ContinuousGrain == getmode(rs_ele$ContinuousGrain),
                       PrAb == getmode(rs_ele$PrAb),
                       LogArea == first(new_data_ele$LogArea),
                       Quality == getmode(rs_ele$Quality))


## make a data frame!
## DATA FRAME STRUCTURE:
## Source (study ID)
## Observed mean study-level shift 
## Predicted study-level shift 
corr_rs_lat <- select(rs_lat, Source, study_level_shift) %>%
  mutate(PredSLShift = new_data_lat$pred_val, 
         Gradient = "Latitudinal")
corr_rs_ele <- select(rs_ele, Source, study_level_shift) %>%
  mutate(PredSLShift =  new_data_ele$pred_val,
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

## logic check:
rs_data %>%
  select(ShiftR, PredSLShift, SLDiff, CorrShift) %>% View()


#######################
####    plotting   ####
#######################
hist(rs_data$CorrShift)

## plot raw versus corrected shift 
rs_data %>%
  ggplot(., aes(x = ShiftR, y = CorrShift, colour = Source)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Observed range shift", y = "Corrected range shift") +
  geom_abline(intercept = 0, slope = 1)

rs_data %>%
  gather(key = "shift_type", value = "shift_value", 
         c(ShiftR, CorrShift)) %>% 
  ggplot(., aes(x = shift_value, fill = shift_type)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(Gradient~Position) +
  scale_fill_discrete(labels = c("Corrected shift", "Observed shift")) +
  labs(x = "Range shift", y = "Frequency", fill = "")

rs_data %>%
  ggplot(., aes(x = SLDiff)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(Gradient~Position) +
  labs(x = "Study-level difference", y = "Frequency")

