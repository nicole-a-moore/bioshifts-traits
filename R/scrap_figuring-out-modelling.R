mod.full.lat.LE <- lmer(shift ~ v.lat.mean * BodyLengthMean + dur + id.area + data + 
                          sampling + uncertainty_parameter + uncertainty_distribution + 
                          (v.lat.mean * BodyLengthMean|class)+(1|article_id),
                        data = lat_le,
                        REML = TRUE,
                        na.action = "na.fail")

mod.full.lat.LE <- lmer(lag ~ expect_tracking*Group + Duration + Area + PrAb + 
                          Sampling + Quality + Signif + (1|article_id),
                        data = lat_le,
                        REML = TRUE,
                        na.action = "na.fail")

lags

f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus) + (1|StartF) + (1|Signif)"

lat <- filter(lags, 
              Group %in% c("Birds", "Plants", "Fish") & Gradient == "Latitudinal") %>%
  mutate(Group = factor(.$Group, ordered = F))
ele <- filter(lags, 
              Group %in% c("Plants") & Gradient == "Elevation") 

lat_le <- filter(lat, Position == "Leading edge")

lat_le$Reference = as.factor(lat_le$Reference)

no_fish <- filter(lat_le, Group !=)

mod.full.lat.LE <- lmer(raw_lag ~ expect_tracking*Group + (1|Sampling)  + (1|PrAb) 
                        + (1|AreaF) + (1|Quality) ,
                        data = lat_le,
                        REML = TRUE,
                        na.action = "na.fail")
summary(mod.full.lat.LE)

data = expand.grid(expect_tracking = unique(lat_le$expect_tracking), 
                   Group = unique(lat_le$Group),
                   Sampling = unique(lat_le$Sampling),
                   PrAb = unique(lat_le$PrAb),
                   Grain = unique(lat_le$Grain)) 

pred = predict(mod.full.lat.LE, data)

data$pred = pred

data %>%
  ggplot(aes(x = expect_tracking, y = pred, colour = Sampling)) +
  geom_point() +
  facet_grid(Group~PrAb~Grain)

lat_le %>%
  group_by(expect_tracking, Group, Sampling, PrAb, AreaF, Quality, Grain) %>%
  tally() %>% 
  View()




## after accounting for methods, how much residual variation is explained by expect_tracking?
## ugh so many singular fits


## look at methdological vars 
lags %>%
  ggplot(aes(x = PrAb, y = raw_lag)) +
  geom_boxplot() +
  facet_grid(Position~Gradient)

lags %>%
  ggplot(aes(x = Sampling, y = raw_lag)) +
  geom_boxplot() +
  facet_grid(Position~Gradient)

lags %>%
  ggplot(aes(x = Grain, y = raw_lag)) +
  geom_boxplot() +
  facet_grid(Position~Gradient)

lags %>%
  ggplot(aes(x = Quality, y = raw_lag)) +
  geom_boxplot() +
  facet_grid(Position~Gradient)

lags %>%
  ggplot(aes(x = AreaF, y = raw_lag)) +
  geom_boxplot() +
  facet_grid(Position~Gradient)

lags %>%
  ggplot(aes(x = NtaxaF, y = raw_lag)) +
  geom_boxplot() +
  facet_grid(Position~Gradient)

lags %>%
  ggplot(aes(x = StartF, y = raw_lag)) +
  geom_boxplot() +
  facet_grid(Position~Gradient)

lags %>%
  ggplot(aes(x = Signif, y = raw_lag)) +
  geom_boxplot() +
  facet_grid(Position~Gradient)

lags %>%
  ggplot(aes(x = as.factor(Reference), y = raw_lag)) +
  geom_boxplot() +
  facet_grid(Position~Gradient)



## see if expect_tracking is correlated with methodological variables:
f="ShiftR ~ PrAb + Sampling + Grain + Quality + AreaF + NtaxaF + (1|Family/Genus) + (1|StartF) + (1|Signif)"

vars = select(lat, expect_tracking,
              PrAb, Sampling, Grain, Quality, AreaF, NtaxaF, StartF, Signif)

chi2 = sapply(vars, FUN = function(x) {
  chi = chisq.test(vars$expect_tracking, x, correct=FALSE)
  return(c(chi$statistic, chi$p.value))
})




chi = chisq.test(lat_le$expect_tracking, lat_le$PrAb, correct=FALSE)
c(chi$statistic, chi$p.value)
chi$p.value < 0.05

chi = chisq.test(lat_le$expect_tracking, lat_le$Sampling, correct=FALSE)
c(chi$statistic, chi$p.value)
chi$p.value < 0.05

chi = chisq.test(lat_le$expect_tracking, lat_le$Grain, correct=FALSE)
c(chi$statistic, chi$p.value)
chi$p.value < 0.05

chi = chisq.test(lat_le$expect_tracking, lat_le$Quality, correct=FALSE)
c(chi$statistic, chi$p.value)
chi$p.value < 0.05

chi = chisq.test(lat_le$expect_tracking, lat_le$AreaF, correct=FALSE)
c(chi$statistic, chi$p.value)
chi$p.value < 0.05

chi = chisq.test(lat_le$expect_tracking, lat_le$NtaxaF, correct=FALSE)
c(chi$statistic, chi$p.value)
chi$p.value < 0.05

chi = chisq.test(lat_le$expect_tracking, lat_le$StartF, correct=FALSE)
c(chi$statistic, chi$p.value)
chi$p.value < 0.05

chi = chisq.test(lat_le$expect_tracking, lat_le$Signif, correct=FALSE)
c(chi$statistic, chi$p.value)
chi$p.value < 0.05
