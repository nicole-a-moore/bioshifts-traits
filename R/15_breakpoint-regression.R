## applying breakpoint regression 
library(tidyverse)
library(segmented)

select <- dplyr::select

theme_set(theme_bw())

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

## for now, calculate the mean and max dispersal potential within species 
v1 = v1 %>%
  group_by(Species) %>%
  mutate(MeanDispersalPotentialKmY = mean(DispersalPotentialKmY)) %>%
  mutate(MeanDispersalPotentialmY = mean(DispersalPotentialmY)) %>%
  mutate(MaxDispersalPotentialKmY = max(DispersalPotentialKmY)) %>%
  mutate(MaxDispersalPotentialmY = max(DispersalPotentialmY)) %>%
  ungroup() %>%
  mutate(CodeOfMax = MaxDispersalPotentialKmY == DispersalPotentialKmY) %>%
  filter(CodeOfMax == TRUE) %>%
  filter(Field != "ArithmeticMeanBreedingDispersal") %>% ## remove breeding dispersal estimates 
  select(-CodeOfMax)

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
         -DispersalDistancem, -Code, 
         #-Field, 
         -Sex, -DispersalUnit, -Database, -DispersalSource,
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

v1 %>%
  filter(Group == "Birds") %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, fill = Field)) +
  geom_histogram() +
  scale_x_log10()

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
mod_data$ShiftKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$ShiftR / 1000, mod_data$ShiftR)
mod_data$ClimVeloTKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$EleVeloT / 1000, mod_data$LatVeloT)
mod_data$AnnualDispPotKmY <- ifelse(mod_data$Gradient == "Elevation", mod_data$annual_dispersal_pot / 1000, mod_data$annual_dispersal_pot)

mycolours <- colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(10)

mod_data <- select(mod_data, -c(SLDiff, CorrShift, PredSLShift, study_level_shift)) %>%
  distinct()

## make plots to summarize subset 
nrow(mod_data) # 1919 range shifts
length(unique(mod_data$Species)) # 404 species 

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


#########################################################
####     fit segmented linear mixed effect model     ####
#########################################################
#----------------------
## allow breakpoint to vary between climate velocity bins
## hold slopes left and right of the breakpoint constant across climate velocity
## test whether:
## - left slope is 1 
## - right slope is 0
## - intercept is 0 
## - breakpoint increases with climate velocity (if possible?)

## get only latitudinal observations 
lat = filter(mod_data, Gradient == "Latitudinal")
write.csv(lat, "data-processed/model-data_lat.csv", row.names = FALSE)

## choose 4 climate velocity bins so that enough species with different dispersal abilities are sampled across different climate velocities
hist(lat$ClimVeloTKmY) ## right skewed
q = quantile(lat$ClimVeloTKmY, probs = c(0,0.25, 0.5, 0.75,1))
lat$ClimVeloTKmY_cont <- lat$ClimVeloTKmY # save original climate velocity as new variable
lat$ClimVeloTKmY = cut(lat$ClimVeloTKmY,
                       breaks = q,
                       include.lowest = T)
plot(lat$ClimVeloTKmY)

lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\[", "(") 
lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\]", ")") 

## calculate max, min and mean climate velocity of each quantile
lat <- lat %>%
  mutate(quant_max = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,2], "\\)", " "))) %>%
  mutate(quant_min = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,1], "\\(", " "))) %>%
  mutate(quant_mean = (quant_min + quant_max)/2) %>%
  group_by(ClimVeloTKmY) %>%
  mutate(quant_median = median(ClimVeloTKmY_cont)) %>%
  ungroup()

## plot distribution of climate velocities within each quantile
lat %>%
  ggplot(aes(x = ClimVeloTKmY_cont)) + geom_histogram() +
  facet_grid(~ClimVeloTKmY) +
  geom_vline(aes(xintercept = quant_mean), colour = "red") +
  geom_vline(aes(xintercept = quant_median), colour = "blue")

most_spec = filter(lat, Source == "A19_P2")

most_spec %>%
  ggplot(aes(x = AnnualDispPotKmY, y = ShiftKmY)) +
  geom_point() +
  # scale_x_log10() +
  # coord_trans(x = "log") + 
  geom_vline(aes(xintercept = ClimVeloTKmY_cont)) + 
  geom_hline(aes(yintercept = ClimVeloTKmY_cont)) + 
  geom_abline(slope = 1, intercept = 0)

## get quantiles 
quants <- lat %>%
  ungroup() %>%
  dplyr::select(quant_mean, quant_min, quant_max, ClimVeloTKmY) %>%
  distinct()

mycol <- rev(colorRampPalette(RColorBrewer::brewer.pal(6, "RdBu"))(4))

## plot by binned climate velocity 
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

##  fit linear mixed effect model:
mod_lat <- lme(ShiftKmY ~ AnnualDispPotKmY, 
               random = ~1|ClimVeloTKmY,
               data = lat)

## and feed to breakpoint regression model:
fit = segmented(mod_lat, 
                seg.Z = ~ AnnualDispPotKmY,
                z.psi = ~ ClimVeloTKmY_cont, ## add continuous climate velocity as a covariate to breakpoint param
                ## allow breakpoint, but not slopes, to vary by climate velocity bin:
                random = list(ClimVeloTKmY = pdDiag(~1 + AnnualDispPotKmY + G0)),
                ## set starting breakpoint values as mean climate velocity across sites:
                psi = mean(lat$ClimVeloTKmY_cont)) 

fit ## high variance of G0 as a fixed effect = lots of variation in change point 
summary(fit$lme.fit)
## intercept: 0.5475239

#plot(fit)
## low variance among slope U (when slope is also allowed to vary)
## high variance among breakpoints G0
fit$psi.i
mean(fit$psi.i)
## all are reliable (within covariate range) 

## get slopes
left = fit$lme.fit$coefficients$fixed[2]
right = fit$lme.fit$coefficients$fixed[2] +  fit$lme.fit$coefficients$fixed[3]
## left = 0.793884
## right = 0.0001327778

AIC(mod_lat, fit$lme.fit)
## random bp model is best fit 

## plot the residuals 
hist(fit$lme.fit$residuals)
## not normal - right skewed 
qqnorm(fit$lme.fit$residuals)
qqline(fit$lme.fit$residuals, col = "steelblue", lwd = 2)

## plot residuals versus independent var. and make sure there is no structure
df <- data.frame(resid = as.numeric(fit$lme.fit$residuals),
                 disp_pot = lat$AnnualDispPotKmY,
                 clim_velo = lat$ClimVeloTKmY,
                 fitted = as.numeric(fit$lme.fit$fitted))


df %>%
  ggplot(aes(x = disp_pot, y = resid)) + geom_point() +
  scale_x_log10() +
  geom_hline(yintercept = 0, colour = "red")

df %>%
  ggplot(aes(y = resid, x = clim_velo)) + geom_boxplot() +
  geom_hline(yintercept = 0, colour = "red")

## fitted vs. residuals 
df %>%
  ggplot(aes(x = fitted, y = resid)) + geom_point() +
  scale_x_log10() +
  geom_hline(yintercept = 0, colour = "red")
## heteroscedastic - higher residual error for high fitted values 

## get confidence intervals 
ci = intervals(zfit$lme.fit)
ci_leftslope = ci$fixed[2,]
ci_rightslope = ci$fixed[3,] + ci$fixed[2,]
ci_intercept = ci$fixed[1,] 
ci_breakpoints = ci$fixed[4,]

low = ci_breakpoints[2] - ci_breakpoints[1]
up = ci_breakpoints[3] - ci_breakpoints[2]

## calculate the y coordinates of the breakpoints 
breakpoints <- as.numeric(fit$psi.i)
slope_left <- slope(fit)[1,1]
slope_right <- slope(fit)[2,1]
intercept <- fit$lme.fit$coefficients$fixed[1]
intercepts <- fit$lme.fit$coefficients$random$id[,1] + intercept

y_coords <- slope_left*breakpoints + intercepts

## plot 
df <- data.frame(theoretical_bp = sort(unique(lat$quant_mean)),
                 model_fitted_bp = y_coords) %>%
  left_join(lat, ., by = c("quant_mean" = "theoretical_bp"))

## median
df %>%
  ggplot(aes(y = model_fitted_bp, x = quant_median)) +
  geom_point() +
  geom_linerange(aes(ymax = model_fitted_bp + up, ymin = model_fitted_bp - low)) +
  geom_point(data = lat, aes(x = quant_median, y = ClimVeloTKmY_cont), 
             size = 0.25, colour = "red") +  ## small red points = plot all climate velocities
  geom_point(aes(y = quant_median), colour = "red") + ## large red point =  median climate velocity within quantile bins 
  facet_wrap(~ClimVeloTKmY, nrow = 1) +
  labs(x = 'Climate velocity', y = "Breakpoint") 
## range of theoretical intercepts (red) is lower than detected breakpoint 
## this means: species are shifting faster than expected 

#----------------------
## plotting theoretical predictions against data and model predictions 
## note: there is no prediction function for segmented lme yet, so I'm on my own here

## make prediction data frame 
pred <- data.frame(expand_grid(AnnualDispPotKmY = seq(min(lat$AnnualDispPotKmY), max(lat$AnnualDispPotKmY), by = 0.01), 
                   ClimVeloTKmY = unique(lat$ClimVeloTKmY)))

## attach breakpoints for each quantile 
breakpoints <- data.frame(fitted_breakpoint = as.numeric(fit$psi.i),
                          ClimVeloTKmY = sort(unique(lat$ClimVeloTKmY)),
                          theoretical_breakpoint = sort(unique(lat$quant_median)),
                          intercept = intercepts)
pred = left_join(pred, breakpoints)

## calculate predictions 
pred = pred %>%
  mutate(predShiftKmY = ifelse(AnnualDispPotKmY < fitted_breakpoint,
                               slope_left*AnnualDispPotKmY + intercept,
                               ifelse(AnnualDispPotKmY > fitted_breakpoint,
                                      slope_right*AnnualDispPotKmY + intercept + slope_left*fitted_breakpoint,
                                      NA))) ## y = mx+ b, where m is left slope if AnnualDispPotKmY < breakpoint, m is right slope if AnnualDispPotKmY > breakpoint

## plot 
pred %>%
  ggplot(aes(x = AnnualDispPotKmY, y = predShiftKmY)) +
  geom_line() +
  facet_grid(~ClimVeloTKmY) 

## on log x axis:
pred %>%
  ggplot(aes(x = AnnualDispPotKmY, y = predShiftKmY)) +
  geom_line() +
  facet_grid(~ClimVeloTKmY) +
  scale_x_log10()

## alongside theoretical prediction:
pred = pred %>%
  mutate(theoreticalShiftKmY = ifelse(AnnualDispPotKmY < theoretical_breakpoint,
                               1*AnnualDispPotKmY + 0,
                               ifelse(AnnualDispPotKmY > theoretical_breakpoint,
                                      0*AnnualDispPotKmY + 0 + 1*theoretical_breakpoint,
                                      NA))) 
  
pred %>%
  ggplot(aes(x = AnnualDispPotKmY, y = predShiftKmY)) +
  geom_line() +
  facet_grid(~ClimVeloTKmY) +
  scale_x_log10() +
  geom_line(aes(x = AnnualDispPotKmY, y = theoreticalShiftKmY), colour = "red")

## with raw data: 
pred %>%
  ggplot(aes(x = AnnualDispPotKmY, y = predShiftKmY)) +
  geom_point(data = lat, 
             aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY_cont), 
             alpha = 0.5) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_line(linewidth = 1) + ## solid line = model fitted relationship
  facet_grid(~ClimVeloTKmY) +
  scale_x_log10() +
  geom_line(aes(x = AnnualDispPotKmY, y = theoreticalShiftKmY), colour = "red",
            linewidth = 1, linetype = "dotted") +  ## dotted line = theoretical prediction
  labs(x = "Potential dispersal rate (km/y)", y = "Range expansion rate (km/y)",
       colour = "Climate velocity (km/y)")


###############################################################################
####      analyzing influence of outliers / points with high leverage      ####
###############################################################################
## assess model fit
## https://online.stat.psu.edu/stat462/node/171/
## reset data
lat = filter(mod_data, Gradient == "Latitudinal")

## choose 4 climate velocity bins so that enough species with different dispersal abilities are sampled across different climate velocities
hist(lat$ClimVeloTKmY) ## right skewed
q = quantile(lat$ClimVeloTKmY, probs = c(0,0.25, 0.5, 0.75,1))
lat$ClimVeloTKmY_cont <- lat$ClimVeloTKmY # save original climate velocity as new variable
lat$ClimVeloTKmY = cut(lat$ClimVeloTKmY,
                       breaks = q,
                       include.lowest = T)
plot(lat$ClimVeloTKmY)

lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\[", "(") 
lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\]", ")") 

## calculate max, min and mean climate velocity of each quantile
lat <- lat %>%
  mutate(quant_max = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,2], "\\)", " "))) %>%
  mutate(quant_min = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,1], "\\(", " "))) %>%
  mutate(quant_mean = (quant_min + quant_max)/2) %>%
  group_by(ClimVeloTKmY) %>%
  mutate(quant_median = median(ClimVeloTKmY_cont)) %>%
  ungroup()

## which range shift points have high leverage in the model?
## unfortuantely there is no framework for checking this in for a segmented.lme object :(
## but, there is a good description of how it can be done in Nieuwenhuis et al.:
## DOI: 10.32614/RJ-2012-011

## the procedure:
## remove points one-by-one from data, each time refitting the same model 
## calculate the influence of each point on the model parameters by comparing full model to model with each point removed 
## Cook's distance can be used to calculate influence on all parameters at once 
## DFBETAS can be used to calculate influence on a particular parameter 

## we care how each point affects:
## - the right / left slopes
## - the breakpoint
## - the intercept

## DFBETAS = difference in the magnitude of the parameter estimate between the model including and the model excluding the case, divided by the standard error of the parameter estimate excluding the case

##  fit linear mixed effect model to full data:
mod_lat <- lme(ShiftKmY ~ AnnualDispPotKmY, 
               random = ~1|ClimVeloTKmY,
               data = lat)

summary(mod_lat)

## and feed to breakpoint regression model:
fit = segmented(mod_lat, 
                seg.Z = ~ AnnualDispPotKmY,
                z.psi = ~ ClimVeloTKmY_cont, ## add continuous climate velocity as a covariate to breakpoint param
                ## allow breakpoint, but not slopes, to vary by climate velocity bin:
                random = list(ClimVeloTKmY = pdDiag(~1 + AnnualDispPotKmY + G0)),
                ## set starting breakpoint values as mean climate velocity across sites:
                psi = mean(lat$ClimVeloTKmY_cont)) 

## get params of full model
params <- fit$lme.fit$coefficients$fixed

## now: loop through points, removing one by one and fitting model

## note: must assign new data frame to old 'lat' object, or else weird error occurs 
## so must refresh data each time 
results <- c()
x = 1
while(x < nrow(lat)) {
  lat = filter(mod_data, Gradient == "Latitudinal")
  
  ## choose 4 climate velocity bins so that enough species with different dispersal abilities are sampled across different climate velocities
  q = quantile(lat$ClimVeloTKmY, probs = c(0,0.25, 0.5, 0.75,1))
  lat$ClimVeloTKmY_cont <- lat$ClimVeloTKmY # save original climate velocity as new variable
  lat$ClimVeloTKmY = cut(lat$ClimVeloTKmY,
                         breaks = q,
                         include.lowest = T)
  lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\[", "(") 
  lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\]", ")") 
  
  ## remove point x
  lat <- lat[-x,]
  
  ##  fit linear mixed effect model to data without the point:
  mod_lat_sub <- lme(ShiftKmY ~ AnnualDispPotKmY, 
                     random = ~1|ClimVeloTKmY,
                     data = lat)
  
  ## and feed to breakpoint regression model:
  fit_sub = segmented(mod_lat, 
                  seg.Z = ~ AnnualDispPotKmY,
                  z.psi = ~ ClimVeloTKmY_cont, ## add continuous climate velocity as a covariate to breakpoint param
                  ## allow breakpoint, but not slopes, to vary by climate velocity bin:
                  random = list(ClimVeloTKmY = pdDiag(~1 + AnnualDispPotKmY + G0)),
                  ## set starting breakpoint values as mean climate velocity across sites:
                  psi = mean(lat$ClimVeloTKmY_cont)) 

  ## get magnitude of parameter estimates
  params_sub <- fit_sub$lme.fit$coefficients$fixed
  
  ## get standard error of parameter estimate excluding case 
  se_sub <- sqrt(diag(vcov(fit_sub)))
  
  ## calculate DFBETAS
  dfbetas = abs(params - params_sub)/se_sub
  
  ## get the parameter estimates 
  breakpoints_sub <- as.numeric(fit_sub$psi.i)
  slope_left_sub <- slope(fit_sub)[1,1]
  slope_right_sub <- slope(fit_sub)[2,1]
  intercepts_sub <- fit_sub$lme.fit$coefficients$fixed[1] + fit_sub$lme.fit$coefficients$random$id[,1]
  subreg = fit_sub$lme.fit$coefficients$fixed[5]
 
  param_ests = data.frame(breakpoint = breakpoints_sub,
                          left_slope = rep(slope_left_sub, 4),
                          left_right = rep(slope_right_sub, 4),
                          intercept = intercepts_sub, 
                          subregression_slope = subreg,
                          ClimVeloTKmY = sort(unique(lat$ClimVeloTKmY)),
                          point_excluded = x)
  
  results_mini <- full_join(param_ests, data.frame(point_excluded = x, 
                                              param = names(params),
                                              param_full = params,
                                              param_sub = params_sub,
                                              se_sub = se_sub,
                                              dfbetas = dfbetas),
                       by = "point_excluded")
  
  ## save results:
  results <- rbind(results, results_mini)
  
  print(paste0("On obs. number: ", x))
  
  x = x + 1
}

## save: 
# write.csv(results, "data-processed/breakpoint-regression-leverage-analysis-results.csv", row.names = FALSE)
results <- read.csv("data-processed/breakpoint-regression-leverage-analysis-results.csv")

## calculate cutoff: 2 / root(n)
cutoff <- 2/sqrt(nrow(lat))

results$cutoff = cutoff

## how many are bigger than cutoff?
results %>%
  mutate(is_bigger = cutoff < dfbetas) %>%
  #filter(is_bigger == TRUE) 
  ggplot(aes(x = is_bigger)) + geom_bar()

## how big are the actual parameter differences?
results %>%
  mutate(param_diff = abs(param_full - param_sub)) %>%
  ggplot(aes(x = param_diff)) + geom_histogram()
## most are small 
## some are large 

## which parameters are affected most?
results %>%
  mutate(param_diff = abs(param_full - param_sub)) %>%
  ggplot(aes(x = param_diff)) + geom_histogram() + 
  facet_grid(~param) 
## mostly the breakpoint param is affected 

## how different are the actual breakpoints from the full model?
results %>%
  select(breakpoint, point_excluded, ClimVeloTKmY) %>%
  distinct() %>%
  ggplot(aes(x = breakpoint)) + geom_histogram() +
  facet_grid(~ClimVeloTKmY) +
  geom_vline(data = real_bps, aes(xintercept = breakpoint), colour = "red")

## and to theoretical breakpoints?
results %>%
  select(breakpoint, ClimVeloTKmY, point_excluded) %>%
  distinct()
  
## plot breakpoints from models with 1 point removed alongside theoretical breakpoints (climate velocity)
results %>% 
  ggplot(aes(y = breakpoint, x = ClimVeloTKmY)) +
  geom_point() + ## large black dots = breakpoints from models with 1 point removed 
  facet_wrap(~ClimVeloTKmY, nrow = 1) +
  geom_point(data = lat, aes(x = ClimVeloTKmY, y = ClimVeloTKmY_cont), 
             size = 0.25, colour = "red") + ## small red points = plot all climate velocities
  geom_point(data = lat, aes(y = quant_median), colour = "red") ## large red point =  median climate velocity within quantile bins 

## plot variation in slope of subregression (breakpoint location ~ climate velocity)
results %>% 
  filter(param == "G.ClimVeloTKmY_cont") %>%
  select(subregression_slope, param_full, param_sub) %>%
  distinct() %>%
  ggplot(aes(x = param_sub)) +
  geom_histogram() +
  geom_vline(aes(xintercept = param_full), colour = "red")


#############################################################################
####     fit segmented linear mixed effect model - leading edge only     ####
#############################################################################
#----------------------
## fit to only leading edge data, not centroid 

## get only latitudinal observations 
lat = filter(mod_data, Gradient == "Latitudinal" & Position == "Leading edge")

## choose 4 climate velocity bins so that enough species with different dispersal abilities are sampled across different climate velocities
hist(lat$ClimVeloTKmY) ## right skewed
q = quantile(lat$ClimVeloTKmY, probs = c(0,0.25, 0.5, 0.75,1))
lat$ClimVeloTKmY_cont <- lat$ClimVeloTKmY # save original climate velocity as new variable
lat$ClimVeloTKmY = cut(lat$ClimVeloTKmY,
                       breaks = q,
                       include.lowest = T)
plot(lat$ClimVeloTKmY)

lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\[", "(") 
lat$ClimVeloTKmY <- str_replace_all(lat$ClimVeloTKmY, "\\]", ")") 

## calculate max, min and mean climate velocity of each quantile
lat <- lat %>%
  mutate(quant_max = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,2], "\\)", " "))) %>%
  mutate(quant_min = as.numeric(str_replace_all(str_split_fixed(ClimVeloTKmY, "\\,", 2)[,1], "\\(", " "))) %>%
  mutate(quant_mean = (quant_min + quant_max)/2) 

## plot distribution of climate velocities within each quantile
lat %>%
  ggplot(aes(x = ClimVeloTKmY_cont)) + geom_histogram() +
  facet_grid(~ClimVeloTKmY) +
  geom_vline(aes(xintercept = quant_mean), colour = "red")

## get quantiles 
quants <- lat %>%
  ungroup() %>%
  dplyr::select(quant_mean, quant_min, quant_max, ClimVeloTKmY) %>%
  distinct()

mycol <- rev(colorRampPalette(RColorBrewer::brewer.pal(6, "RdBu"))(4))

## plot by binned climate velocity 
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

##  fit linear mixed effect model:
mod_lat <- lme(ShiftKmY ~ AnnualDispPotKmY, 
               random = ~1|ClimVeloTKmY,
               data = lat)

summary(mod_lat)

## and feed to breakpoint regression model:
fit = segmented(mod_lat, 
                seg.Z = ~ AnnualDispPotKmY,
                ## allow breakpoint, but not slopes, to vary by climate velocity bin:
                random = list(ClimVeloTKmY = pdDiag(~1 + G0)),
                ## set starting breakpoint values as mean climate velocity across sites:
                psi = mean(lat$ClimVeloTKmY_cont)) 

fit ## less variance in G0 as a fixed effect = lots of variation in change point 
summary(fit$lme.fit)
## intercept: 1.083429


plot(fit)
## low variance among slope U (when slope is also allowed to vary)
## high variance among breakpoints G0
fit$psi.i
mean(fit$psi.i)
## all are reliable (within covariate range) 

## get slopes
left = fit$lme.fit$coefficients$fixed[2]
right =  fit$lme.fit$coefficients$fixed[2] +  fit$lme.fit$coefficients$fixed[3]
## left = 0.1544394
## right = -0.002850327 

AIC(mod_lat, fit$lme.fit)
## random bp model is best fit 

## plot the residuals 
hist(fit$lme.fit$residuals)
## more normal than model fit to full dataset  
qqnorm(fit$lme.fit$residuals)
qqline(fit$lme.fit$residuals, col = "steelblue", lwd = 2)
## distribution still has long tails

## plot residuals versus independent var. and make sure there is no structure
df <- data.frame(resid = as.numeric(fit$lme.fit$residuals),
                 disp_pot = lat$AnnualDispPotKmY,
                 clim_velo = lat$ClimVeloTKmY,
                 fitted = as.numeric(fit$lme.fit$fitted))

df %>%
  ggplot(aes(x = disp_pot, y = resid)) + geom_point() +
  scale_x_log10() +
  geom_hline(yintercept = 0, colour = "red")

df %>%
  ggplot(aes(y = resid, x = clim_velo)) + geom_boxplot() +
  geom_hline(yintercept = 0, colour = "red")

## fitted vs. residuals 
df %>%
  ggplot(aes(x = fitted, y = resid)) + geom_point() +
  scale_x_log10() +
  geom_hline(yintercept = 0, colour = "red")
## heteroscedastic - higher residual error for high fitted values 

## get confidence intervals 
ci = intervals(fit$lme.fit)
ci_leftslope = ci$fixed[2,]
ci_rightslope = ci$fixed[3,] + ci$fixed[2,]
ci_intercept = ci$fixed[1,] 
ci_breakpoints = ci$fixed[4,]

low = ci_breakpoints[2] - ci_breakpoints[1]
up = ci_breakpoints[3] - ci_breakpoints[2]

## calculate the y coordinates of the breakpoints 
breakpoints <- as.numeric(fit$psi.i)
slope_left <- slope(fit)[1,1]
slope_right <- slope(fit)[2,1]
intercept <- fit$lme.fit$coefficients$fixed[1]
intercepts <- fit$lme.fit$coefficients$random$id[,1] + intercept

y_coords <- slope_left*breakpoints + intercepts

## plot 
df <- data.frame(theoretical_bp = sort(unique(lat$quant_mean)),
                 model_fitted_bp = y_coords) %>%
  left_join(lat, ., by = c("quant_mean" = "theoretical_bp"))


df %>%
  ggplot(aes(y = model_fitted_bp, x = quant_mean)) +
  geom_point() +
  geom_linerange(aes(ymax = model_fitted_bp + up, ymin = model_fitted_bp - low)) +
  geom_point(data = lat, aes(x = quant_mean, y = ClimVeloTKmY_cont), 
             size = 0.25, colour = "red") +  ## small red points = plot all climate velocities
  geom_point(aes(y = quant_mean), colour = "red") + ## large red point =  mean climate velocity within quantile bind 
  facet_wrap(~ClimVeloTKmY) +
  labs(x = 'Climate velocity', y = "Breakpoint") 
## range of theoretical intercepts (red) is within detected breakpoint confidence intervals

#----------------------
## plotting theoretical predictions against data and model predictions 
## note: there is no prediction function for segmented lme yet, so I'm on my own here

## make prediction data frame 
pred <- data.frame(expand_grid(AnnualDispPotKmY = seq(min(lat$AnnualDispPotKmY), max(lat$AnnualDispPotKmY), by = 0.01), 
                               ClimVeloTKmY = unique(lat$ClimVeloTKmY)))

## attach breakpoints for each quantile 
breakpoints <- data.frame(fitted_breakpoint = as.numeric(fit$psi.i),
                          ClimVeloTKmY = sort(unique(lat$ClimVeloTKmY)),
                          theoretical_breakpoint = sort(unique(lat$quant_mean)),
                          intercept = intercepts)
pred = left_join(pred, breakpoints)

## calculate predictions 
pred = pred %>%
  mutate(predShiftKmY = ifelse(AnnualDispPotKmY < fitted_breakpoint,
                               slope_left*AnnualDispPotKmY + intercept,
                               ifelse(AnnualDispPotKmY > fitted_breakpoint,
                                      slope_right*AnnualDispPotKmY + intercept + slope_left*fitted_breakpoint,
                                      NA))) ## y = mx+ b, where m is left slope if AnnualDispPotKmY < breakpoint, m is right slope if AnnualDispPotKmY > breakpoint

## plot 
pred %>%
  ggplot(aes(x = AnnualDispPotKmY, y = predShiftKmY)) +
  geom_line() +
  facet_grid(~ClimVeloTKmY) 

## on log x axis:
pred %>%
  ggplot(aes(x = AnnualDispPotKmY, y = predShiftKmY)) +
  geom_line() +
  facet_grid(~ClimVeloTKmY) +
  scale_x_log10()

## alongside theoretical prediction:
pred = pred %>%
  mutate(theoreticalShiftKmY = ifelse(AnnualDispPotKmY < theoretical_breakpoint,
                                      1*AnnualDispPotKmY + 0,
                                      ifelse(AnnualDispPotKmY > theoretical_breakpoint,
                                             0*AnnualDispPotKmY + 0 + 1*theoretical_breakpoint,
                                             NA))) 

pred %>%
  ggplot(aes(x = AnnualDispPotKmY, y = predShiftKmY)) +
  geom_line() +
  facet_grid(~ClimVeloTKmY) +
  scale_x_log10() +
  geom_line(aes(x = AnnualDispPotKmY, y = theoreticalShiftKmY), colour = "red")

## with raw data: 
pred %>%
  ggplot(aes(x = AnnualDispPotKmY, y = predShiftKmY)) +
  geom_point(data = lat, 
             aes(x = AnnualDispPotKmY, y = ShiftKmY, colour = ClimVeloTKmY_cont), 
             alpha = 0.5) +
  scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
  geom_line(linewidth = 1) + ## solid line = model fitted relationship
  facet_grid(~ClimVeloTKmY) +
  scale_x_log10() +
  geom_line(aes(x = AnnualDispPotKmY, y = theoreticalShiftKmY), colour = "red",
            linewidth = 1, linetype = "dotted") +  ## dotted line = theoretical prediction
  labs(x = "Potential dispersal rate (km/y)", y = "Range expansion rate (km/y)",
       colour = "Climate velocity (km/y)")


## make up some data with breakpoints to see how the model is being fit 
## breakpoint increases with level, intercept stays at 0
data1 <- data.frame(x = seq(1:60),
                   y = c(seq(1:20), rep(0+20, 40)))
data2 <- data.frame(x = seq(1:60),
                    y = c(seq(1:30), rep(0+30, 30)))
data3 <- data.frame(x = seq(1:60),
                    y = c(seq(1:40), rep(0+40, 20)))
data4 <- data.frame(x = seq(1:60),
                    y = c(seq(1:50), rep(0+50, 10)))

data = rbind(data1, data2) %>%
  rbind(., data3) %>%
  rbind(., data4)

data$level = rep(1:4, each = 60)

## add noise to avoid singular fit 
data$y = data$y + rnorm(240, mean = 0, sd = 0.1)

ggplot(data, aes(x = x, y = y, colour = level)) + geom_point()

##  fit linear mixed effect model:
mod<- lme(y ~ x, 
          random = ~1|level,
          data = data)

summary(mod)

## and feed to breakpoint regression model:
fit = segmented(mod, 
                seg.Z = ~ x,
                ## allow breakpoint, but not slopes, to vary by level:
                random = list(level = pdDiag(~1 + G0)),
                ## set starting breakpoint values as mean climate velocity across sites:
                psi = mean(data$y),
                data = data) 

fit
summary(fit)

fit$psi.i
## break points are: 20, 30, 40, 50

slope(fit)
## left slope 1, right slope 0

fit$lme.fit$coefficients$fixed[1] + fit$lme.fit$coefficients$random$id[,1]
## intercepts: 0, 0, 0, 0

mean(fit$lme.fit$coefficients$fixed[1] + fit$lme.fit$coefficients$random$id[,1])
## overall intercept: 0

## what happens if we have imperfect sampling across x?
## remove points where x < 10 from level 1
## prediction: intercept for this level will increase, overall intercept will increase

data <- data %>%
  filter(!(x < 20 & level == 1))

##  fit linear mixed effect model:
mod<- lme(y ~ x, 
          random = ~1|level,
          data = data)

summary(mod)

## and feed to breakpoint regression model:
fit = segmented(mod, 
                seg.Z = ~ x,
                ## allow breakpoint, but not slopes, to vary by level:
                random = list(level = pdDiag(~1 + G0)),
                ## set starting breakpoint values as mean climate velocity across sites:
                psi = mean(data$y),
                data = data) 

summary(fit)

fit$psi.i
mean(fit$psi.i)
## break points are: 20, 30, 40, 50
## mean is 34.98

slope(fit)
## left slope 1, right slope 0

fit$lme.fit$coefficients$fixed[1] + fit$lme.fit$coefficients$random$id[,1]
## intercepts: 0, 0, 0, 0

mean(fit$lme.fit$coefficients$fixed[1] + fit$lme.fit$coefficients$random$id[,1])
## overall intercept: 0








## garbage 
breakpoints <- as.numeric(fit$psi.i)
slope_left <- slope(fit)[1,1]
slope_right <- slope(fit)[2,1]
intercepts <- fit$lme.fit$coefficients$fixed[1] + fit$lme.fit$coefficients$random$id[,1]

breakpoints_sub <- as.numeric(fit_sub$psi.i)
slope_left_sub <- slope(fit_sub)[1,1]
slope_right_sub <- slope(fit_sub)[2,1]
intercepts_sub <- fit_sub$lme.fit$coefficients$fixed[1] + fit_sub$lme.fit$coefficients$random$id[,1]

diff_bps = abs(breakpoints - breakpoints_sub)
diff_lslope = abs(slope_left - slope_left_sub)
diff_rslope = abs(slope_right - slope_right_sub)
diff_ints = abs(intercepts - intercepts_sub)



## notes:
## then test for independence and normality of residuals and ensure variance is evenly distributed 
## bootstrap to get confidence intervals 

### plot the breakpoints 
## calculate y value of intercepts 
## try to see if we can put climate velocity in as a continuous variable 
## slope less than one interpretation = even when species have dispersal capacity to keep up, they are lagging for other reasons 
## try and figure out why this version's slope is 0.5
## previous appraoch but binned would probably give more similar answer to this?
## change original model to deal with zero boundedness if we have a problem with the residuals 

## decide what to do about non-normal residuals 



  


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
length(xlev1) ## 60 points identified 
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
# Residual standard error: 3.186 on 1274 degrees of freedom
# Multiple R-squared:  0.03951,	Adjusted R-squared:  0.03876 

## fit breakpoint regression
mod_bp_lev <- segmented(mod_lat, 
                        seg.Z = ~ AnnualDispPotKmY) ## do not set any starting value for break point

summary(mod_bp_lev)
# Residual standard error: 2.94 on 1272 degrees of freedom
# Multiple R-Squared: 0.1834,  Adjusted R-squared: 0.1815 

## get estimated breakpoint and its standard error
mod_bp_lev$psi ## 5.65

## get the slopes
slope(mod_bp_lev) ## 0.512, 0.0018

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
# Residual standard error: 2.086 on 1301 degrees of freedom
# Multiple R-squared:  0.02247,	Adjusted R-squared:  0.02172 

## fit breakpoint regression
mod_bp_out <- segmented(mod_lat, 
                        seg.Z = ~ AnnualDispPotKmY) ## do not set any starting value for break point

summary(mod_bp_out)
# Residual standard error: 1.814 on 1299 degrees of freedom
# Multiple R-Squared: 0.2621,  Adjusted R-squared: 0.2604 

## get estimated breakpoint and its standard error
mod_bp_out$psi ## 3.455

## get the slopes
slope(mod_bp_out) ## 0.473, 0.00

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











## garbage 
########################################################
####      fitting global breakpoint regression      ####
########################################################
## filter to shifts > 0.0001
mod_data <- mod_data %>%
  filter(ShiftR > 0.0001) 

## split between elev and lat
lat = filter(mod_data, Gradient == "Latitudinal")
ele = filter(mod_data, Gradient == "Elevation")

## remove outlier
# hist(lat$ShiftKmY)
# lat <- filter(lat, ShiftKmY < 20)

## fit normal regression 
hist(lat$ShiftKmY)

## must log response variable because otherwise non-normal residuals
## but do we need to? do we care if residuals are non-normal? 
## - if we take the log, we aren't really testing our hypothesis - no reason to believe a non-linear relationship exists 
## - structure in residuals affects predictions and error estimates, but we only care about slope estimates and whether they align with our predictions
## - plus - we already *know* this model will be lacking explanatory variables that could explain this leftover structure in the residuals (all the other things affecting range shifts)

lat$logShiftKmY = log(lat$ShiftKmY)

mod_lat <- lm(ShiftKmY ~ AnnualDispPotKmY, data = lat)
summary(mod_lat)
# Residual standard error: 3.229 on 1334 degrees of freedom
# Multiple R-squared:  0.01402,	Adjusted R-squared:  0.01328 

## fit breakpoint regression
mod_bp_lat <- segmented(mod_lat, 
                        seg.Z = ~ AnnualDispPotKmY) ## do not set any starting value for break point

summary(mod_bp_lat)
# Residual standard error: 2.946 on 1332 degrees of freedom
# Multiple R-Squared: 0.1806,  Adjusted R-squared: 0.1787 

## get estimated breakpoint and its standard error
mod_bp_lat$psi ## 5.03

## get the slopes
slope(mod_bp_lat) ## 0.591, -0.0002
exp(0.6927)
exp(0.00019)
exp(1)^mod_bp_lat$coefficients[1] ## intercept of 0.171
# for every one unit increase in dispersal rate, range shift increases by a factor of 1.99

## plot the residuals
hist(mod_bp_lat$residuals)

## plot the predictions
pred <- data.frame(AnnualDispPotKmY = lat$AnnualDispPotKmY)
pred <- data.frame(predict(mod_bp_lat, newdata = pred, interval = "confidence"))
ci <- confint(mod_bp_lat)

bp_y <- predict(mod_bp_lat, newdata = data.frame(AnnualDispPotKmY = mod_bp_lat$psi[,2]))

df <- data.frame(AnnualDispPotKmY = lat$AnnualDispPotKmY,
                 logShiftKmY = lat$logShiftKmY,
                 pred_shift = pred$fit,
                 min = min(lat$ClimVeloTKmY), 
                 max = max(lat$ClimVeloTKmY),
                 mean = mean(lat$ClimVeloTKmY), 
                 breakpoint_x = mod_bp_lat$psi[,2],
                 breakpoint_y = bp_y,
                 ci_lower = pred$lwr,
                 ci_upper = pred$upr,
                 ci_bp_lower = ci[,2],
                 ci_bp_upper = ci[,3])

## plot the fitted model
df %>%
  ggplot(aes(x = AnnualDispPotKmY, y = exp(1)^pred_shift)) +
  #geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey70", alpha = 0.7) +
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
  geom_pointrange(aes(xmin = ci_bp_lower, xmax = ci_bp_upper, x = breakpoint_x, y = exp(1)^breakpoint_y),
                  size = 0.3, colour = "black") +
  geom_vline(aes(xintercept = mean))  # plot theoretical point
# geom_rect(aes(xmin = min, xmax = max, ymin = 0, ymax = 41),
#           alpha = 0.002) 

# df %>%
#   ggplot(aes(x = AnnualDispPotKmY, y = exp(1)^pred_shift)) +
#   #geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey70", alpha = 0.7) +
#   geom_point(data = lat, aes(colour = ClimVeloTKmY, y = ShiftKmY), alpha = 0.5) + 
#   geom_line() +
#   theme_bw() +
#   stat_function(colour = "black", # add 1:1 line
#                 fun = function(x){x},
#                 linetype = "dashed") +
#   scale_y_continuous(limits = c(0, 41), 
#                      expand = c(0.1, 0.1)) +
#   scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
#                 labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
#   scale_colour_gradient2(high = "#B2182B", low = "#2166AC", mid = "#F8DCCB", midpoint = 3.5) +
#   geom_pointrange(aes(xmin = ci_bp_lower, xmax = ci_bp_upper, x = breakpoint_x, y = exp(1)^breakpoint_y),
#             size = 0.3, colour = "black") +
#   geom_vline(aes(xintercept = mean))  # plot theoretical point 
#   # geom_rect(aes(xmin = min, xmax = max, ymin = 0, ymax = 41),
#   #           alpha = 0.002) 

########################################################
####      fitting global breakpoint regression      ####
####          (only leading edge obs.)              ####
########################################################
## reset data 
lat = filter(mod_data, Gradient == "Latitudinal")

## remove outlier
# lat <- filter(lat, ShiftKmY < 20)

## remove centroid obs 
lat <- filter(lat, Position != "Centroid")

## fit normal regression 
mod_lat <- lm(ShiftKmY ~ AnnualDispPotKmY, data = lat)
summary(mod_lat)
# Residual standard error: 3.382 on 249 degrees of freedom
# Multiple R-squared:  0.01715,	Adjusted R-squared:  0.0132 

## fit breakpoint regression
mod_bp_lat_le <- segmented(mod_lat, 
                           seg.Z = ~ AnnualDispPotKmY) ## do not set any starting value for break point

summary(mod_bp_lat_le)
# Residual standard error: 2.913 on 247 degrees of freedom
# Multiple R-Squared: 0.2767,  Adjusted R-squared: 0.2679 

# ## test for normality of residuals 
# res <- residuals(mod_bp_lat)
# hist(res) ## a few large residuals
# 
# fitted = fitted(mod_bp_lat)
# res_first <- res[which(lat$AnnualDispPotKmY <= mod_bp_lat$psi[,2])]
# res_sec <- res[which(lat$AnnualDispPotKmY > mod_bp_lat$psi[,2])]
# fitted_first <- fitted[which(lat$AnnualDispPotKmY <= mod_bp_lat$psi[,2])]
# fitted_sec <- fitted[which(lat$AnnualDispPotKmY > mod_bp_lat$psi[,2])]

# qqnorm(res_first)
# qqline(res_first)
# qqnorm(res_sec)
# qqline(res_sec)
# 
# shapiro.test(res_first)
# shapiro.test(res_sec)
# ## not normal :(
# 
# data.frame(res = res, fitted = fitted(mod_bp_lat)) %>%
#   ggplot(., aes(x = fitted, y = res)) +
#   geom_point() +
#   scale_x_log10()

# ## test for independence of residuals 
# # plot residuals versus disp potential
# data.frame(res = res, dispersal_pot = lat$AnnualDispPotKmY) %>%
#   ggplot(., aes(x = dispersal_pot, y = res)) +
#   geom_point() +
#   scale_x_log10()
# ## i think it is because bird range shifts are more variable than plants
# ## you can also see the effect of the lower bound on range shift rate here (0)

# ## check variance of residuals 
# copy <- lat
# copy$res2 <- res^2
# 
# mod_res <- lm(res2 ~ AnnualDispPotKmY, data = copy) #ANOVA of the squared residuals
# anova(mod_res) #displays the results


## get estimated breakpoint and its standard error
mod_bp_lat_le$psi ## 5.73

## get the slopes
slope(mod_bp_lat_le) ## 0.691, -0.0017

## plot the predictions
pred <- data.frame(AnnualDispPotKmY = lat$AnnualDispPotKmY)
pred <- data.frame(predict(mod_bp_lat_le, newdata = pred, interval = "confidence"))
ci <- confint(mod_bp_lat_le)

bp_y <- predict(mod_bp_lat_le, newdata = data.frame(AnnualDispPotKmY = mod_bp_lat_le$psi[,2]))

df <- data.frame(AnnualDispPotKmY = lat$AnnualDispPotKmY,
                 ShiftKmY = lat$ShiftKmY,
                 pred_shift = pred$fit,
                 min = min(lat$ClimVeloTKmY), 
                 max = max(lat$ClimVeloTKmY),
                 mean = mean(lat$ClimVeloTKmY), 
                 breakpoint_x = mod_bp_lat_le$psi[,2],
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
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c('0.0001', "0.001", "0.01","0.1", "1","10", "100", "1000")) +
  stat_function(colour = "black", # add 1:1 line
                fun = function(x){x},
                linetype = "dashed") +
  scale_y_continuous(limits = c(0, 41), 
                     expand = c(0.1, 0.1)) +
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
# Residual standard error: 0.8736 on 404 degrees of freedom
# Multiple R-squared:  0.03039,	Adjusted R-squared:  0.02799 
summary(mod_fits[[1]])
# Residual standard error: 0.7004 on 402 degrees of freedom
# Multiple R-Squared: 0.3799, Adjusted R-squared: 0.3753 

summary(mod_lms[[2]])
# Residual standard error: 4.838 on 272 degrees of freedom
# Multiple R-squared:  0.01643,	Adjusted R-squared:  0.01281 
summary(mod_fits[[2]])
# Residual standard error: 4.562 on 270 degrees of freedom
# Multiple R-Squared: 0.1318,  Adjusted R-squared: 0.1221 

summary(mod_lms[[3]])
# Residual standard error: 3.249 on 333 degrees of freedom
# Multiple R-squared:  0.006419,	Adjusted R-squared:  0.003435 
summary(mod_fits[[3]])
# Residual standard error: 3.149 on 331 degrees of freedom
# Multiple R-Squared: 0.07206,  Adjusted R-squared: 0.06365 

summary(mod_lms[[4]])
# Residual standard error: 2.657 on 319 degrees of freedom
# Multiple R-squared:  0.001132,	Adjusted R-squared:  -0.001999 
summary(mod_fits[[4]])
# Residual standard error: 2.514 on 317 degrees of freedom
# Multiple R-Squared: 0.1116,  Adjusted R-squared: 0.1032 

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
# Residual standard error: 1.293 on 83 degrees of freedom
# Multiple R-squared:  0.0314,	Adjusted R-squared:  0.01973 
summary(mod_fits[[1]])
# Residual standard error: 0.9052 on 81 degrees of freedom
# Multiple R-Squared: 0.537,  Adjusted R-squared: 0.5198 

summary(mod_lms[[2]])
# Residual standard error: 4.586 on 86 degrees of freedom
# Multiple R-squared:  0.01926,	Adjusted R-squared:  0.00786 
summary(mod_fits[[2]])
# Residual standard error: 4.251 on 84 degrees of freedom
# Multiple R-Squared: 0.1766,  Adjusted R-squared: 0.1472 

summary(mod_lms[[3]])
# Residual standard error: 2.235 on 76 degrees of freedom
# Multiple R-squared:  0.01036,	Adjusted R-squared:  -0.002666 
summary(mod_fits[[3]])
# Residual standard error: 1.802 on 74 degrees of freedom
# Multiple R-Squared: 0.373,  Adjusted R-squared: 0.3476


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



