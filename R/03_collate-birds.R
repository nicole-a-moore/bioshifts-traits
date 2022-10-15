## script to check availability of traits of interest for birds
## renames columns in primary trait databases to common names
library(tidyverse)

## read in key so we know which columns = which common names 
key <- read.csv("data-processed/traits-by-group/birds_filled.csv")

## filter out ones without common names 
## NOTE: excluding MarineSpeciesTraits database for now because weird and confusing  
## NOTE: also excluding trait bank - ask lise about it
key <- key %>%
  filter(!is.na(trait_common_name)) %>% 
  filter(!Database_name == "Marine Species traits") %>%
  filter(!Database_name == "TraitBank")

unique(key$trait_common_name)


## read in each bioshifts trait database subset one by one
db_names <- unique(key$Database_name)
filenames <- c("amniota", "AnAge", "AVONET", "BIOTIME", "EltonTraits", "GARD_appendix_1", 
               "GARD_novosolov_et_al_2017",
               "Heinen_et_al_2017", "lislevand_et_al_2007",
               #"MarineSpeciesTraits",
                "Pigot_et_al_2020", "Storchov_and_Horak_2018", "Sutherland_2000_birds"
               #, "TraitBank"
               )
all_dbs <- c()
n <- 1
while(n <= length(unique(key$Database_name))) {
  ## get name of current database 
  cur_db_name <- db_names[n]
  
  ## read in current database 
  db <- read_csv(paste("data-processed/primary-traits-bioshifts/", filenames[n], ".csv", sep = ""))
  
  ## filter key to cur db traits 
  db_key <- filter(key, Database_name == cur_db_name)
  
  ## make sure all trait names in key match those in database 
  nmissing <- length(which(!db_key$trait %in% colnames(db))) 
  
  ## if they are not, print a warning message
  if(nmissing != 0) {
    print(paste("Warning! ", cur_db_name, " has ", nmissing, " columns that do not match the key!", sep = ""))
  }
  
  ## get rid of non-common columns, rename columns to common trait name  
  order <- sapply(1:nrow(db_key), 
                  FUN = function(x) {return(which(db_key$trait[x] == colnames(db)))})
  
  cols_to_keep <- c("genus_species", "genus", "family", "class", 
                    "order", "kingdom", "db", "db_code")
  
  ## getting rid of non-common columns except for taxonomy
  db <- db %>%
    select(cols_to_keep, colnames(db)[order])
  
  ## renaming columns to common names 
  colnames(db)[(length(cols_to_keep)+1):ncol(db)] <- db_key$trait_common_name
  
  ## add empty columns for missing common traits to allow rbinding later 
  missing_cols <- unique(key$trait_common_name)[which(!unique(key$trait_common_name) %in% colnames(db))]
  
  db[,ncol(db):(ncol(db)+length(missing_cols))] <- NA
  
  colnames(db)[(ncol(db) - length(missing_cols)+1):ncol(db)] <- missing_cols
  
  ## add primary_source column to keep track of where each trait comes from
  db$primary_source <- cur_db_name
  
  ## bind to others 
  all_dbs <- rbind(all_dbs, db)
  
  n = n + 1
}

write.csv(all_dbs, "data-processed/common-traits/birds_prelim.csv", row.names = FALSE)


all_dbs <- read_csv("data-processed/common-traits/birds_prelim.csv") %>%
  filter(kingdom != "Plantae")

## resolve taxonomy issues 
## species that were not found by taxize have some taxonomic issues to sort out 
tax_data <- select(all_dbs, genus_species, order, family, class, kingdom) %>%
  unique(.) %>%
  filter(genus_species %in% .$genus_species[which(duplicated(.$genus_species))]) 

sub <- all_dbs[which(all_dbs$genus_species %in% tax_data$genus_species),]

## change "Fringillidae." to "Fringillidae"
sub <- sub %>%
  mutate(family = ifelse(family == "Fringillidae.",
                         "Fringillidae",
                         family))

## resolve species that have different values for family:
multifam <- sub %>%
  group_by(genus_species) %>%
  select(genus_species, family) %>%
  filter(!is.na(family)) %>%
  unique(.) %>%
  tally(.) %>%
  filter(n > 1)

multifam <- sub[which(sub$genus_species %in% multifam$genus_species),] 
multifam <- select(multifam, genus_species, family) %>%
  arrange(genus_species)

multifam2 <- multifam %>%
  mutate(actual_genus_species = ifelse(genus_species == "Amphispiza belli",
                                       "Artemisiospiza belli",
                                       NA)) %>%
  mutate(actual_genus_species = ifelse(genus_species == "Aimophila cassinii", ## ITIS
                                       "Peucaea cassinii",
                                       actual_genus_species)) %>%
  mutate(actual_genus_species = ifelse(genus_species == "Monarcha trivirgatus", ## ITIS
                                       "Symposiachrus trivirgatus",
                                       actual_genus_species)) %>%
  mutate(actual_genus_species = ifelse(genus_species == "Pitohui nigrescens", ## ITIS
                                       "Melanorectes nigrescens",
                                       actual_genus_species)) %>%
  mutate(actual_genus_species = ifelse(genus_species == "Spizella arborea", ## ITIS
                                       "Spizelloides arborea",
                                       actual_genus_species)) 

## search for taxonomy
multifam_key <- harmonize_taxnomy(unique(multifam2$actual_genus_species))[[1]]
multifam_key$genus_species_old <- unique(multifam2$genus_species)


## fill the rest by group
therest <- sub %>%
  group_by(genus_species) %>%
  filter(!genus_species %in% multifam$genus_species) %>%
  fill(family, kingdom, class, order, .direction = "updown") %>%
  select(genus_species, family, order, class, kingdom, db, db_code) %>%
  rename("genus_species_accepted" = genus_species) %>%
  mutate(genus_species_old = NA)

## combine database of new taxonomy
multifam <- multifam_key %>%
  select(genus_species, family, order, class, kingdom, db, db_code, genus_species_old) %>%
  rename("genus_species_accepted" = genus_species) %>%
  rbind(., therest) %>%
  mutate(genus_species_old = ifelse(is.na(genus_species_old),
                                    genus_species_accepted,
                                    genus_species_old))

## combine with original data 
without_sp <- filter(all_dbs, !genus_species %in% multifam$genus_species_old & 
                       !genus_species %in% multifam$genus_species_accepted)
sp_mf <- filter(all_dbs, genus_species %in% multifam$genus_species_old)

sp_mf <- sp_mf %>%
  select(-family, -order, -class, -kingdom, -db, -db_code) %>%
  left_join(., multifam, by = c("genus_species" = "genus_species_old")) %>%
  mutate(genus_species = genus_species_accepted) %>%
  select(-genus_species_accepted)

all_dbs <- rbind(sp_mf, without_sp)


## group by species and resolve issues 
## result: database where every species has one row, so that missing value actually means we have no value for that species
all_dbs %>%
  group_by(genus_species) %>%
  tally() %>%
  View()

## if there are more than one body mass/length measurements per species, keep the maximum
all_dbs %>%
  group_by(genus_species, kingdom) %>%
  summarize(adult_body_mass_g = paste(adult_body_mass_g,  collapse = ", ")) %>%
  View()

all_dbs %>%
  group_by(genus_species, kingdom) %>%
  summarize(female_body_mass_g = paste(female_body_mass_g,  collapse = ", ")) %>%
  View()


## gather all body mass measurements
## if there are multiple per species, keep the maximum one
resolved <- all_dbs %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(female_body_mass_g = max(female_body_mass_g, na.rm = TRUE),
         male_body_mass_g = max(male_body_mass_g, na.rm = TRUE),
         adult_body_mass_g = max(adult_body_mass_g, na.rm = TRUE),
         unsexed_body_mass_g = max(unsexed_body_mass_g, na.rm = TRUE),
         unk_body_mass_g = max(unk_body_mass_g, na.rm = TRUE),
         female_body_mass_at_maturity_g = max(female_body_mass_at_maturity_g, na.rm = TRUE)) %>%
  ungroup() %>%
  unique(.)

resolved2 <- resolved %>%
  gather(key = "measurement_type", value = "body_mass",
         c(female_body_mass_g, male_body_mass_g, adult_body_mass_g,
           unk_body_mass_g, unsexed_body_mass_g, female_body_mass_at_maturity_g)) %>%
  group_by(genus_species,  kingdom, class, order, family) %>%
  mutate(body_mass_g = max(body_mass, na.rm = TRUE)) %>%
  mutate(body_mass_g_source = ifelse(!is.infinite(body_mass_g),
                                     paste(unique(primary_source[which(body_mass == body_mass_g)]), 
                                        collapse = ", "),
                                     NA)) %>%
  mutate(body_mass_g_measurement_type = ifelse(!is.infinite(body_mass_g),
                                               paste(unique(measurement_type[which(body_mass == body_mass_g)]), 
                                        collapse = ", "),
                                        NA)) %>%
  select(-body_mass, -measurement_type) %>%
  ungroup() %>%
  mutate(body_mass_g = ifelse(is.infinite(body_mass_g), NA, body_mass_g)) %>%
  unique(.)

## do the same with length measurements 
which(str_detect(colnames(resolved2), "length"))

resolved3 <- resolved2 %>%
  mutate(unsexed_body_length_mm = unsexed_body_length_cm*1000) %>%
  gather(key = "measurement_type", value = "body_length",
         c(unsexed_body_length_mm, 
           body_length_mm, wing_length_mm,
           mean_unsexed_wing_length_mm,
           mean_male_wing_length_mm,
           mean_female_wing_length_mm,
           unsexed_wing_length_mm,
           male_wing_length_mm,
           female_wing_length_mm,
           qualitative_body_size)) %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(body_length_mm = max(body_length, na.rm = TRUE)) %>%
  mutate(body_length_mm_source = ifelse(!is.infinite(body_length_mm),
                                        paste(unique(primary_source[which(body_length == body_length_mm)]), 
                                        collapse = ", "),
                                        NA)) %>%
  mutate(body_length_mm_measurement_type = ifelse(!is.infinite(body_length_mm),
                                                  paste(unique(measurement_type[which(body_length == body_length_mm)]), 
                                                  collapse = ", "),
                                                  NA)) %>%
  select(-body_length, -measurement_type, -unsexed_body_length_cm) %>%
  ungroup() %>%
  mutate(body_length_mm = ifelse(is.infinite(body_length_mm), NA, body_length_mm)) %>%
  unique(.)


## do the same with litter/clutch size measurements 
which(str_detect(colnames(resolved3), "litter"))

resolved4 <- resolved3 %>%
  gather(key = "measurement_type", value = "litter_clutch_size",
         c(litter_or_clutch_size_n, litter_or_clutch_size_n_min, 
           litter_or_clutch_size_n_max, litter_or_clutch_size_n_mean)) %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(litter_or_clutch_size_n = max(litter_clutch_size, na.rm = TRUE)) %>%
  mutate(litter_or_clutch_size_n_source = ifelse(!is.infinite(litter_or_clutch_size_n),
                                                 paste(unique(primary_source[which(litter_clutch_size == 
                                                                               litter_or_clutch_size_n)]), 
                                           collapse = ", "),
                                           NA)) %>%
  mutate(litter_or_clutch_size_n_measurement_type = ifelse(!is.infinite(litter_or_clutch_size_n),
                                                           paste(unique(measurement_type[which(litter_clutch_size == 
                                                                                     litter_or_clutch_size_n)]), 
                                                     collapse = ", "),
                                                     NA)) %>%
  select(-litter_clutch_size, -measurement_type) %>%
  ungroup() %>%
  mutate(litter_or_clutch_size_n = ifelse(is.infinite(litter_or_clutch_size_n), NA, litter_or_clutch_size_n)) %>%
  unique(.)

## do the same for dispersal distance
resolved5 <- resolved4 %>%
  gather(key = "measurement_type", value = "dd_km",
         c(median_dispersal_distance_km, maximum_dispersal_distance_km)) %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(dispersal_distance_km = max(dd_km, na.rm = TRUE)) %>%
  mutate(dispersal_distance_km_source = ifelse(!is.infinite(dispersal_distance_km),
                                               paste(unique(primary_source[which(dd_km == 
                                                                            dispersal_distance_km)]), 
                                                    collapse = ", "),
                                               NA)) %>%
  mutate(dispersal_distance_km_measurement_type = ifelse(!is.infinite(dispersal_distance_km),
                                                         paste(unique(measurement_type[which(dd_km ==
                                                                                        dispersal_distance_km)]), 
                                                              collapse = ", "),
                                                         NA)) %>%
  select(-dd_km, -measurement_type) %>%
  ungroup() %>%
  mutate(dispersal_distance_km = ifelse(is.infinite(dispersal_distance_km), NA, dispersal_distance_km))  %>%
  unique(.)

## do the same for maturity
resolved6 <- resolved5 %>%
  gather(key = "measurement_type", value = "maturity",
         c(female_maturity_d, male_maturity_d, maturity_d)) %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(maturity_d = max(maturity, na.rm = TRUE)) %>%
  mutate(maturity_d_source = ifelse(!is.infinite(maturity_d),
                                    paste(unique(primary_source[which(maturity == maturity_d)]), 
                                              collapse = ", "),
                                    NA)) %>%
  mutate(dispersal_distance_km_measurement_type = ifelse(!is.infinite(maturity_d),
           paste(unique(measurement_type[which(maturity == maturity_d)]), 
                                                        collapse = ", "),
           NA)) %>%
  select(-maturity, -measurement_type) %>%
  ungroup() %>%
  mutate(maturity_d = ifelse(is.infinite(maturity_d), NA, maturity_d)) %>%
  unique(.)

## do the same for longevity
resolved7 <- resolved6 %>%
  gather(key = "measurement_type", value = "longevity",
         c(longevity_y, maximum_longevity_y)) %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(longevity_y = max(longevity, na.rm = TRUE)) %>%
  mutate(longevity_y_source = ifelse(!is.infinite(longevity_y),
                                     paste(unique(primary_source[which(longevity == longevity_y)]), 
                                   collapse = ", "),
                                   NA)) %>%
  mutate(longevity_y_measurement_type = ifelse(!is.infinite(longevity_y),
           paste(unique(measurement_type[which(longevity == longevity_y)]), 
                 collapse = ", "),
           NA)) %>%
  select(-longevity, -measurement_type) %>%
  ungroup() %>%
  mutate(longevity_y = ifelse(is.infinite(longevity_y), NA, longevity_y)) %>%
  unique(.)

##gather trophic info
resolved8 <- resolved7 %>%
  mutate(trophic_guild = ifelse(!is.na(trophic_guild), 
                                paste(trophic_guild, " (", primary_source, ")", sep = ""),
                                NA),
         trophic_niche = ifelse(!is.na(trophic_niche), 
                                paste(trophic_niche, " (", primary_source, ")", sep = ""),
                                NA),
         trophic_level = ifelse(!is.na(trophic_level), 
                                paste(trophic_level, " (", primary_source, ")", sep = ""),
                                NA)) %>%
  gather(key = "measurement_type", value = "trophic",
         c(trophic_guild, trophic_niche, trophic_level)) %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(trophic_group = paste(unique(trophic)[which(!is.na(unique(trophic)))], collapse = ", ")) %>%
  select(-trophic, -measurement_type) %>%
  ungroup() %>%
  unique(.)
  
colnames(resolved8)

## for the rest of the columns, choose max and change source
resolved9 <- resolved8 %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(max_litters_or_clutches_per_y = max(litters_or_clutches_per_y, na.rm = TRUE)) %>%
  mutate(litters_or_clutches_per_y_source = ifelse(!is.infinite(max_litters_or_clutches_per_y),
                                                   paste(unique(primary_source[which(litters_or_clutches_per_y 
                                                                              == max_litters_or_clutches_per_y)]), 
                                    collapse = ", "),
                                    NA)) %>%
  select(-litters_or_clutches_per_y) %>%
  rename("litters_or_clutches_per_y" = max_litters_or_clutches_per_y) %>%
  ungroup() %>%
  mutate(litters_or_clutches_per_y = ifelse(is.infinite(litters_or_clutches_per_y), NA, litters_or_clutches_per_y)) %>%
  unique(.)

resolved10 <- resolved9 %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(max_gestation_d = max(gestation_d, na.rm = TRUE)) %>%
  mutate(gestation_d_source = ifelse(!is.infinite(max_gestation_d),
                                     paste(unique(primary_source[which(gestation_d == max_gestation_d)]), 
                                                  collapse = ", "),
                                     NA)) %>%
  select(-gestation_d) %>%
  rename("gestation_d" = max_gestation_d) %>%
  ungroup() %>%
  mutate(gestation_d = ifelse(is.infinite(gestation_d), NA, gestation_d)) %>%
  unique(.)

resolved11 <- resolved10 %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(max_weaning_d = max(weaning_d, na.rm = TRUE)) %>%
  mutate(weaning_d_source = ifelse(!is.infinite(max_weaning_d),
         paste(unique(primary_source[which(weaning_d == max_weaning_d)]), 
                                    collapse = ", "),
         NA)) %>%
  select(-weaning_d) %>%
  rename("weaning_d" = max_weaning_d) %>%
  ungroup() %>%
  mutate(weaning_d = ifelse(is.infinite(weaning_d), NA, weaning_d)) %>%
  unique(.)

resolved12 <- resolved11 %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(max_growth_rate = max(growth_rate, na.rm = TRUE)) %>%
  mutate(growth_rate_source = ifelse(!is.infinite(max_growth_rate),
                                     paste(unique(primary_source[which(growth_rate == max_growth_rate)]), 
                                  collapse = ", "),
                                  NA)) %>%
  select(-growth_rate) %>%
  rename("growth_rate" = max_growth_rate) %>%
  ungroup() %>%
  mutate(growth_rate = ifelse(is.infinite(growth_rate), NA, growth_rate)) %>%
  unique(.)

resolved13 <- resolved12 %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(min_min_latitude_deg = min(min_latitude_deg, na.rm = TRUE)) %>%
  mutate(min_latitude_source = ifelse(!is.infinite(min_min_latitude_deg),
                                      paste(unique(primary_source[which(min_latitude_deg == min_min_latitude_deg)]), 
                                    collapse = ", "),
                                    NA)) %>%
  select(-min_latitude_deg) %>%
  rename("min_latitude_deg" = min_min_latitude_deg) %>%
  ungroup() %>%
  mutate(min_latitude_deg = ifelse(is.infinite(min_latitude_deg), NA, min_latitude_deg)) %>%
  unique(.)

resolved14 <- resolved13 %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(max_max_latitude_deg = max(max_latitude_deg, na.rm = TRUE)) %>%
  mutate(max_latitude_deg_source = ifelse(!is.infinite(max_max_latitude_deg),
                                          paste(unique(primary_source[which(max_latitude_deg == max_max_latitude_deg)]), 
                                    collapse = ", "),
                                    NA)) %>%
  select(-max_latitude_deg) %>%
  rename("max_latitude_deg" = max_max_latitude_deg) %>%
  ungroup() %>%
  mutate(max_latitude_deg = ifelse(is.infinite(max_latitude_deg), NA, max_latitude_deg)) %>%
  unique(.)

resolved15 <- resolved14 %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(max_range_size_km2 = max(range_size_km2, na.rm = TRUE)) %>%
  mutate(max_latitude_deg_source =  ifelse(!is.infinite(max_range_size_km2),
                                           paste(unique(primary_source[which(range_size_km2 == max_range_size_km2)]), 
                                         collapse = ", "),
                                         NA)) %>%
  select(-range_size_km2) %>%
  rename("range_size_km2" = max_range_size_km2) %>%
  ungroup() %>%
  mutate(range_size_km2 = ifelse(is.infinite(range_size_km2), NA, range_size_km2)) %>%
  unique(.)

resolved16 <- resolved15 %>%
  group_by(genus_species, kingdom, class, order, family) %>%
  mutate(max_population_density_n_per_ha = max(population_density_n_per_ha, na.rm = TRUE)) %>%
  mutate(population_density_n_per_ha_source = ifelse(!is.infinite(max_population_density_n_per_ha),
                                                     paste(unique(primary_source[which(population_density_n_per_ha == 
                                                                       max_population_density_n_per_ha)]), 
                                         collapse = ", "),
                                         NA)) %>%
  select(-population_density_n_per_ha) %>%
  rename("population_density_n_per_ha" = max_population_density_n_per_ha) %>%
  ungroup() %>%
  mutate(population_density_n_per_ha = ifelse(is.infinite(population_density_n_per_ha), 
                                              NA, population_density_n_per_ha)) %>%
  unique(.)

## get rid of primary source column and non-unique rows
data_final <- resolved16 %>%
  select(-primary_source) %>%
  unique(.) %>%
  mutate(trophic_group =  ifelse(trophic_group == "", NA,
                                 trophic_group)) ## make sure empty trophic position is NA

## now if I did everything right there should be one unique row per species 
## check:
length(unique(data_final$genus_species))
nrow(data_final) ## 2507 species 

## next: see how many species are missing each trait
trait_cols <-  colnames(data_final)[9:41]
trait_cols <- trait_cols[-c(which(str_detect(trait_cols, "source")))]
trait_cols <- trait_cols[-c(which(str_detect(trait_cols, "measurement_type")))]

data_final %>%
  gather(key = "trait", value = "measurement", trait_cols) %>%
  mutate(is_na = ifelse(is.na(measurement), 
                        "Missing", 
                        "Not missing")) %>%
  ggplot(aes(y = trait, fill = is_na)) + geom_bar() 

## would be good to also visualize the measurement types for each trait 


## and next look at and clean all of the different values of each trait

