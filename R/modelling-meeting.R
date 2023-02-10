## script to fit models using the two different approaches
library(tidyverse)
library(ggplot2)
library(lme4)
library(PNWColors)
library(MuMIn)
library(effects)
library(rr2)
pal = pnw_palette("Bay",7)

## read in dispersal data 
data <- read.csv( "data-processed/dispersal-data_modelling-meeting.csv")
## reorder factor 
data$expect_tracking = factor(data$expect_tracking, levels = c("No", "Yes"), ordered = TRUE)

data %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, y = abs(ShiftR))) + geom_point() +
  scale_x_log10() + scale_y_log10() +
  geom_abline(intercept = 0, slope = 1)

data %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = MaxDispersalPotentialmY, y = abs(ShiftR))) + geom_point() +
  scale_x_log10() + scale_y_log10() +
  geom_abline(intercept = 0, slope = 1)

data %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, y = abs(corrected_shift))) + geom_point() +
  scale_x_log10() + scale_y_log10() +
  geom_abline(intercept = 0, slope = 1)

data %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = MaxDispersalDistanceKm, y = abs(corrected_shift))) + geom_point() +
  scale_x_log10() + scale_y_log10() +
  geom_abline(intercept = 0, slope = 1)


data %>%
  mutate(negshift = ifelse(ShiftR < 0, "neg", "pos")) %>%
  group_by(negshift) %>%
  tally()

data %>%
  mutate(negshift = ifelse(ShiftR < 0, "neg", "pos")) %>%
  select(Species, negshift) %>%
  unique() %>%
  group_by(negshift) %>%
  tally()


## remove species with negative lags 
## see of results change 
## fit new models:
## shift ~ disp*cv + rand eff
## lag ~ disp + rand eff
## (should explain the same thing - does it?)

## lag ~ disp*HFI + rand eff

##############
## data vis ##
##############

## How many species should vs. shouldn't be able to keep up with temperature change?
## Is dispersal potential less limiting across elevation than latitude?
data %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, fill = expect_tracking)) + 
  geom_histogram() +
  theme_classic() +
  scale_fill_manual(values = c(pal[6], pal[4])) + 
  facet_grid(~Gradient) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), 
                labels = c("0.0001","0.001", "0.01", "0.1", "1", "10", "100", "1000")) +
  labs(x = "Annual dispersal potential (km/y)",
       y = "Number of range shift observations", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?")
## Dispersal is more limiting across latitude than elevation

data %>%
  ggplot(aes(x = MaxDispersalPotentialKmY, fill = expect_tracking)) + 
  geom_histogram() +
  theme_bw() +
  scale_fill_manual(values = c(pal[6], pal[4])) + 
  facet_grid(Gradient~Group) +
  labs(x = "Annual dispersal potential (km/y)",
       y = "Number of range shift observations", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000), 
                labels = c("1e-4    ","1e-3", " 0.01", "0.1", "1", "10", "100", "1000"))  +
  theme(panel.grid = element_blank(), strip.background = element_rect(colour="black", fill="white"),
        legend.position = "none")
## All animals have an annual dispersal potential greater than the velocity of temperature across elevation, but not latitude
## In other words, all animals should be able to disperse fast enough to keep up with climate change across elevation 

## Is only the speed of climate change important? Or is it both climate change and dispersal capacity?

data %>%
  ggplot(aes(x = Group, fill = expect_tracking)) + 
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = c(pal[6], pal[4])) + 
  coord_flip() +
  facet_wrap(~Gradient) +
  labs(x = "Taxonomic group", 
       y = "Number of range shift observations", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?")

## Do species that are expect to be able to keep up with climate change have less lags?
## Plot raw lags (raw shift - climate velocity), then plot inferred lags (residual shift - climate velocity)
## Here, more negative value = greater lag behind climate velocity 
lags %>%
  filter(Gradient == "Latitudinal") %>%
  ggplot(aes(x = annual_dispersal_pot, y = raw_lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (km/y)", 
       y = "Inferred range shift lag (km/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Position~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

lags %>%
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
  facet_grid(Position~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

lags %>%
  filter(Gradient == "Elevation") %>%
  ggplot(aes(x = annual_dispersal_pot, y = raw_lag, fill = expect_tracking)) + 
  geom_boxplot() +
  theme_bw() +
  scale_x_log10() +
  scale_fill_manual(values = c(pal[4], pal[6])) + 
  labs(x = "Annual dispersal potential (m/y)", 
       y = "Inferred range shift lag (m/y)", 
       fill = "Expect dispersal\npotential to allow\nspecies to keep up\nwith temp change?") +
  geom_abline(intercept = 0, slope = 0) +
  facet_grid(Position~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 

lags %>%
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
  facet_grid(Position~Group) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) 
## Within groups, species that are expected to be able to keep up with temperature change have greater range expansion lags than species who are not
## Patterns stronger in residual shifts 



#############################
## prep data for modelling ##
#############################
## get rid of groups that have only species who can or only species who can't keep up with climate change 
## separate data into latitudinal and elevation shifts
data$expect_tracking = factor(data$expect_tracking, levels = unique(data$expect_tracking), ordered = F)

lat <- filter(data, 
              Group %in% c("Birds", "Plants", "Fish") & Gradient == "Latitudinal") %>%
  mutate(Group = factor(.$Group, ordered = F))
ele <- filter(data, 
              Group %in% c("Plants") & Gradient == "Elevation") 

#####################################
## approach 1: modelling residuals ##
#####################################
## column "lag" is the residual shift - velocity 

## start with latitude
## model inferred lag as a function of taxonomic group, range shift position, and expect tracking
lat_mod <- lm(lag ~ expect_tracking*Group*Position, data = lat)
summary(lat_mod)

## inspect model fit
plot(lat_mod)
hist(lat_mod$residuals)
lat$residuals = lat_mod$residuals
ggplot(data = lat, aes(x = expect_tracking, y = residuals)) + geom_boxplot()
ggplot(data = lat, aes(x = Group, y = residuals)) + geom_boxplot()
ggplot(data = lat, aes(x = Position, y = residuals)) + geom_boxplot()

## problem: variance of residuals differs by taxonomic group and position
## solution: glm with VarIdent variance structure
## https://link-springer-com.proxy3.library.mcgill.ca/chapter/10.1007/978-0-387-87458-6_4

## compete models with different variance structures
vf <- varIdent(form = ~ 1 | Group)
vf2 <- varIdent(form = ~ 1 | Position)
vf3 <- varIdent(form = ~ 1 | Position*Group)
mod_a <- gls(lag ~ expect_tracking*Group*Position, data = lat)
mod_b <- gls(lag ~ expect_tracking*Group*Position, data = lat, weights = vf)
mod_c <- gls(lag ~ expect_tracking*Group*Position, data = lat, weights = vf2)
mod_d <- gls(lag ~ expect_tracking*Group*Position, data = lat, weights = vf3)
AIC(mod_a, mod_b, mod_c, mod_d)
## the best model is the one that allows variance structure to differ by position AND taxonomic group\
summary(mod_d)

## inspect model fit
plot(mod_d)
hist(mod_d$residuals)
lat$residuals = mod_d$residuals
ggplot(data = lat, aes(x = expect_tracking, y = residuals)) + geom_boxplot()
ggplot(data = lat, aes(x = Group, y = residuals)) + geom_boxplot()
ggplot(data = lat, aes(x = Position, y = residuals)) + geom_boxplot()


## plot the model predictions:
pred <- predict(mod_d, type = "response")

lat_pred = lat 
lat_pred$predicted_val = pred

lat_pred = lat_pred %>%
  select(Group, expect_tracking, Position, predicted_val) %>%
  distinct()

lat_pred$expect_tracking = factor(lat_pred$expect_tracking, levels = c("No", "Yes"),
                                  ordered = T)

## plot
ggplot(data = lat_pred, aes(x = expect_tracking, colour = expect_tracking, y = predicted_val,
                            shape = Position)) + 
  geom_point() +
  facet_wrap(~Group) +
  theme_classic() +
  labs(x = "Expect dispersal potential to allow species to keep up with temp change?",
       y = "Model-predicted inferred range shift lag (km/y)",
       title = "Latitude") +
  scale_colour_manual(values = c(pal[6], pal[4])) +
  theme(legend.position = "none")
## don't know how to get se/conf ints yet for this type of model


## do the same for elevation
## fit simple model and inspect residuals 
ele_mod <- lm(lag ~ expect_tracking*Position, data = ele)
summary(ele_mod)
plot(ele_mod)
hist(ele_mod$residuals)
ele$residuals = ele_mod$residuals
ggplot(data = ele, aes(x = expect_tracking, y = residuals)) + geom_boxplot()
ggplot(data = ele, aes(x = Position, y = residuals)) + geom_boxplot()

## inspect model fit
plot(ele_mod)
hist(ele_mod$residuals)
ele$residuals = ele_mod$residuals
ggplot(data = ele, aes(x = expect_tracking, y = residuals)) + geom_boxplot()
ggplot(data = ele, aes(x = Position, y = residuals)) + geom_boxplot()

## same problem; so compete models with different variance structure
vf2 <- varIdent(form = ~ 1 | Position)
mod_a <- gls(lag ~ expect_tracking*Position, data = ele)
mod_c <- gls(lag ~ expect_tracking*Position, data = ele, weights = vf2)
AIC(mod_a, mod_c)
## the best model is the one that allows variance structure to differ by position
summary(mod_c)

## inspect model fit
plot(mod_c)
hist(mod_c$residuals)
ele$residuals = mod_c$residuals
ggplot(data = ele, aes(x = expect_tracking, y = residuals)) + geom_boxplot()
ggplot(data = ele, aes(x = Position, y = residuals)) + geom_boxplot()


## plot the model predictions:
pred <- predict(mod_c, type = "response")

ele_pred = ele 
ele_pred$predicted_val = pred

ele_pred = ele_pred %>%
  select(Group, expect_tracking, Position, predicted_val) %>%
  distinct()

ele_pred$expect_tracking = factor(ele_pred$expect_tracking, levels = c("No", "Yes"),
                                  ordered = T)

## plot
ggplot(data = ele_pred, aes(x = expect_tracking, colour = expect_tracking, y = predicted_val,
                            shape = Position)) + 
  geom_point() +
  facet_wrap(~Group) +
  theme_classic() +
  labs(x = "Expect dispersal potential to allow species to keep up with temp change?",
       y = "Model-predicted inferred range shift lag (km/y)",
       title = "Elevation") +
  scale_colour_manual(values = c(pal[6], pal[4])) +
  theme(legend.position = "none")



####################################
## approach 2: modelling raw lags ##
####################################
## column "raw_lag" is the raw range shift value - climate velocity 
## here, we will fit one model to each combination of position x gradient

## start with leading edge x latitude:
lat_le <- filter(lat, Position == "Leading edge" & Gradient == "Latitudinal")

## fit model with fixed effects of interest + methods
mod_lat_le = lmer(raw_lag ~ expect_tracking*Group + AreaF + PrAb + Grain +
              Sampling + Quality + Signif + (1|Source),
            data = lat_le,
            REML = TRUE,
            na.action = "na.fail")
summary(mod_lat_le)

## fit model with only fixed effects of interest but same random structure
mod_simp = lmer(raw_lag ~ expect_tracking*Group + (1|Source),
                  data = lat_le,
                  REML = TRUE,
                  na.action = "na.fail")
summary(mod_simp)

## compare conditional and marginal r2
r.squaredGLMM(mod_lat_le)
r.squaredGLMM(mod_simp) ## model without fixed effects methodology 

## compare AIC
AIC(mod_lat_le, mod_simp)
## model with methods has lower AIC

## plot model parameters
ef <- effect("expect_tracking*Group", mod_lat_le)
summary(ef)

plt_dat <- data.frame(ef)

ggplot(plt_dat, aes(x = expect_tracking, y = fit, color = Group)) + geom_point() + 
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), width=0.4) + 
  theme_bw() +
  facet_wrap(~Group) +
  labs(title = "Latitude - leading edge")


## now do the same for centroid x latitude 
lat_o <- filter(lat, Position == "Centroid" & Gradient == "Latitudinal")

mod_lat_o = lmer(raw_lag ~ expect_tracking*Group + AreaF + PrAb + Grain +
              Sampling + Quality + Signif + (1|Source),
            data = lat_o,
            REML = TRUE,
            na.action = "na.fail")
summary(mod_lat_o)

## fit model with only fixed effects of interest but same random structure
mod_simp = lmer(raw_lag ~ expect_tracking*Group + (1|Source),
                data = lat_o,
                REML = TRUE,
                na.action = "na.fail")
summary(mod_simp)

## compare conditional and marginal r2
r.squaredGLMM(mod_lat_o)
r.squaredGLMM(mod_simp) ## model without fixed effects methodology 

## compare AIC
AIC(mod_lat_o, mod_simp)
## model with methods has lower AIC

## plot model parameters
ef <- effect("expect_tracking*Group", mod_lat_o)
summary(ef)

plt_dat <- data.frame(ef)

ggplot(plt_dat, aes(x = expect_tracking, y = fit, color = Group)) + geom_point() + 
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), width=0.4) + 
  theme_bw() +
  facet_wrap(~Group) +
  labs(title = "Latitude - centroid")

## now do the same for leading edge x elevation 
ele_le <- filter(ele, Position == "Leading edge" & Gradient == "Elevation")

mod_ele_le = lmer(raw_lag ~ expect_tracking + AreaF + PrAb + Grain +
             Sampling + Quality + Signif + (1|Source),
           data = ele_le,
           REML = TRUE,
           na.action = "na.fail")
summary(mod_ele_le)

## fit model with only fixed effects of interest but same random structure
mod_simp = lmer(raw_lag ~ expect_tracking + (1|Source),
                data = ele_le,
                REML = TRUE,
                na.action = "na.fail")
summary(mod_simp)

## compare conditional and marginal r2
r.squaredGLMM(mod_ele_le)
r.squaredGLMM(mod_simp) ## model without fixed effects methodology 

## compare AIC
AIC(mod_ele_le, mod_simp)
## model with methods has lower AIC

## plot model parameters
ef <- effect("expect_tracking", mod_ele_le)
summary(ef)

plt_dat <- data.frame(ef)

ggplot(plt_dat, aes(x = expect_tracking, y = fit)) + geom_point() + 
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), width=0.4) + 
  theme_bw() + 
  labs(title = "Elevation - leading edge")

## now do the same for centroid x elevation 
ele_o <- filter(ele, Position == "Centroid" & Gradient == "Elevation")

mod_ele_o = lmer(raw_lag ~ expect_tracking + AreaF + PrAb + Grain +
             Sampling + Quality + Signif + (1|Source),
           data = ele_o,
           REML = TRUE,
           na.action = "na.fail")
summary(mod_ele_o)

## fit model with only fixed effects of interest but same random structure
mod_simp = lmer(raw_lag ~ expect_tracking + (1|Source),
                data = ele_o,
                REML = TRUE,
                na.action = "na.fail")
summary(mod_simp)

## compare conditional and marginal r2
r.squaredGLMM(mod_ele_o)
r.squaredGLMM(mod_simp) ## model without fixed effects methodology 

## compare AIC
AIC(mod_ele_o, mod_simp)
## model with methods has lower AIC

## plot model parameters
ef <- effect("expect_tracking", mod_ele_o)
summary(ef)

plt_dat <- data.frame(ef)

ggplot(plt_dat, aes(x = expect_tracking, y = fit)) + geom_point() + 
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), width=0.4) + 
  theme_bw() + 
  labs(title = "Elevation - leading edge")


## Concerns: approach 1
## 1. Issues raised by Freckleton 2002 - problems arise if fixed effects are colinear with methodological variables used in first model
## 2. Potential issue that a different methdological model is fit for each position x gradient 

## Concerns: approach 2
## 1. how to check for overfitting?

## Concerns: general
## 1. results from approach 1 are not comparable to results from approach 2. Why? Approach 1 accounts for variation due to methodology in the raw shift data (methodological variation accounted for *before* calculating lag), while approach 2 accounts for variation due to methodology in the lag data
## 2. in both models, we sometimes have multiple range shift observations per species and are not accounting for potential non-independence 



## check whether lm results of residuals are different from lme results with methodological factors 

## one big model instead of Romain's many models
## no model selection - but check for covariance
## include position
## then use predict and hold all methodological variables to same level, make separate predictions for different levels of position, marine v terr, etc. (similar approach to modelling residuals)

## when velocity is positive and shift is negative, lag should be value of velocity (so that max lag is climate velocity)

## test for non-linear relationship between annual dispersal potential and lag 

## check dispersal distances and see if they represent maximum dispersal potential in lifetime, check units 






