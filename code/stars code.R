#code for honours project, sorry Claire, it won't be as clean as you'd like

#because of renv will need to redownload some packages


#load in packages----
library(tidyverse)
library(ggforce)
library(glmmTMB)
library(DHARMa)
library(ggeffects)

#load in data----

#belt data + wrangling

belt <- read_csv("data/stars_belt - Sheet1 (1).csv") %>%
  add_count(site, species, name = "abundance")

#get belt site data

belt_site <- read_csv("data/stars_belt - Sheet1 (1).csv") %>%
  group_by(site, date) %>%
  mutate(average_size = mean(size, na.rm = TRUE)) %>%
  group_by(site, date, species, average_size) %>%
  summarise(species_abundance = n(), .groups = "drop") %>%
  group_by(site, date) %>%
  mutate(total_abundance = sum(species_abundance, na.rm = TRUE)) %>%
  ungroup()
            

belt_site_wide <- belt_site %>%
  pivot_wider(
    names_from = species,
    values_from = species_abundance,
    values_fill = 0
  )

#this dataset now gives me the 0's needed for absence

belt_site_long <- belt_site_wide %>%
  pivot_longer(
    cols = c(leather, mottled, ochre, blood, pink, sunflower),
    names_to = "species",
    values_to = "species_abundance"
  ) %>%
  mutate(species = as.factor(species),
         )
#stars rove

stars_rove <- read_csv("data/stars_rove - Sheet1.csv")

#debris data


debris <- read_csv("data/debris_rove - Sheet1 (2).csv") %>%
  mutate(debris_score = as.factor(debris_score),
         mussel_score = as.factor(mussel_score))  %>%
  add_count(site, debris_type, name = "debris_count")
  

#site data

site <- read_csv("data/stars_site - Sheet1 (2).csv") %>%
  mutate(type = factor(type, levels = c("dock", "inter-dock", "no-dock")),
         depth_change = end_depth - start_depth,
         slope = asin(depth_change/transect_length)*180/pi,
         cluster = as.factor(cluster),
         pair = as.factor(pair),
         site = as.factor(site))



#realistically, the data should exist at two levels, within site and at site


#leftjoin/fulljoin time----

#within sites df (UNCOMPLETED YET)

within <- belt %>%
  full_join(site, by = c("site","date"), relationship = "many-to-one")

#among sites df

among <- site %>%
  full_join(stars_rove, by = c("site","date"), relationship = "many-to-many")

among <- among %>%
  full_join(belt_site, by = c("site", "date"), relationship = "many-to-many")

among <- among %>%
  full_join(debris, by = c("site","date"), relationship = "many-to-many")

#get it so only 1 row per observation
among_site <- among %>%
  select(site ,date ,total_abundance, slope, cluster, type, pair, average_size,
         richness, transect_length, debris_score, mussel_score, debris_count) %>%
  distinct()

among_site <- among_site %>%
  mutate(density = total_abundance/(transect_length*2))

#trying it again so I have it with a dataset that includes 0's

trial_zero <- belt_site_long %>%
  full_join(site, by = c("site","date"), relationship = "many-to-many")

trial_zero <- trial_zero %>%
  full_join(stars_rove, by = c("site", "date"), relationship = "many-to-many")

trial_zero <- trial_zero %>%
  full_join(debris, by = c("site", "date"), relationship = "many-to-many")

#get only one row per site

among_zero <- trial_zero %>%
  select(!average_size) %>%
  distinct()

#justifying variable transect length

#quickly do the lm for transect length

a <- lm(total_abundance ~ transect_length, data = among_site)
summary(a)

ggplot(data = among_site, aes(x = transect_length, y = total_abundance)) + geom_point()

#good to have as a back up

ggplot(data = among_site, (aes(x = transect_length, y = total_abundance))) +
  geom_point() + geom_smooth(method = "lm")

simulateResiduals(a, plot = TRUE)

#ok let's do more work to justify it, what substrate are the stars hanging out on ----

belt_j <- read_csv("data/stars_belt - Sheet1.csv (1)") %>%
  drop_na() %>%
 count(site,date,substrate, name = "substrate_abundance") %>%
  group_by(site,date) %>%
  mutate(total_abundance = sum(substrate_abundance),
         substrate_percent = substrate_abundance/total_abundance * 100) %>%
  ungroup()

#summarize some substrate

belt_j <- belt_j %>%
  mutate(
    substrate = case_when(
      substrate %in% c("sand/boulder", "boulder/sand") ~
        "boulder/sand",
      substrate %in% c("cobble/boulder", "boulder/cobble") ~
        "boulder/cobble",
      substrate %in% c("sand/cobble", "cobble/sand") ~
        "cobble/sand",
      TRUE ~ substrate)
  )

#pivot belt j wider so I can have 0's in categories that they are not in

belt_j_wider <- belt_j %>%
  select(site, date, substrate, substrate_abundance) %>%
  group_by(site, date, substrate) %>%
  summarise(
    substrate_abundance = sum(substrate_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = c(site, date),
    names_from = substrate,
    values_from = substrate_abundance,
    values_fill = list(substrate_abundance = 0)
  )

#percent time

belt_j_long <- belt_j_long %>%
  group_by(site, date) %>%
  mutate(
    total_abundance = sum(substrate_abundance, na.rm = TRUE),
    substrate_percent = substrate_abundance /
      total_abundance * 100
  ) %>%
  ungroup()

#pivot it back to normal

belt_j_long <- belt_j_wider %>%
  pivot_longer(
    cols = c(boulder, sand, debris, cobble, `boulder/sand`, `boulder/cobble`, `cobble/sand`),
    names_to = "substrate",
    values_to = "substrate_abundance"
  ) %>%
  mutate(substrate = as.factor(substrate),
  )


#without care for species

ggplot(data = belt_j_long, aes(x = substrate, y = substrate_percent)) + 
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.1,
    dodge.width = 0.3),
    alpha = 0.4,
    size = 2) +
  stat_summary(fun = mean,
               geom = "point",
               position = position_dodge(width = 0.3),
               size = 3) +
  stat_summary(fun.data = mean_cl_normal,
               geom = "errorbar",
               position = position_dodge(width = 0.3),
               width = 0.2) +
  theme_classic()


#what substrate are specific species of stars hanging out on?----

belt_j_species <- read_csv("data/stars_belt - Sheet1.csv (1)") %>%
  drop_na(site, date, species, substrate) %>%
  filter(!species %in% c("sunflower", "pink", "blood")) %>%
  count(site,date,species,substrate,
    name = "substrate_abundance"
  ) %>%
  group_by(site, date, species) %>%
  complete(
    substrate = substrate_levels,
    fill = list(substrate_abundance = 0)
  ) %>%
  mutate(
    species_total = sum(substrate_abundance), 
    substrate_percent = substrate_abundance /species_total * 100
  ) %>%
  ungroup()


belt_j_species <- belt_j_species %>%
  mutate(
    substrate = case_when(
      substrate %in% c("sand/boulder", "boulder/sand") ~
        "boulder/sand",
      substrate %in% c("cobble/boulder", "boulder/cobble") ~
        "boulder/cobble",
      substrate %in% c("sand/cobble", "cobble/sand") ~
        "cobble/sand",
      TRUE ~ substrate)
  )

#including species

ggplot(data = belt_j_species, aes(x = substrate, y = substrate_percent, group = species,
                                  colour = species )) + 
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.1,
    dodge.width = 0.3),
    alpha = 0.4,
    size = 2) +
  stat_summary(fun = mean,
               geom = "point",
               position = position_dodge(width = 0.3),
               size = 3) +
  stat_summary(fun.data = mean_cl_normal,
               geom = "errorbar",
               position = position_dodge(width = 0.3),
               width = 0.2) +
  facet_wrap_paginate(~species) +
  theme_classic() + theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


#ok last justification part, try the density by rocky area----

rock_data <- belt %>%
  filter(rock_keep == 1) %>%
  add_count(site, species, name = "rock_abundance")

rock <- rock_data %>%
  full_join(site, by = c("site","date"), relationship = "many-to-many") %>%
  mutate(rock_density = rock_abundance/(rock_length*2))
  
ggplot(data = rock, aes(x = rock_length, y = rock_density)) +
  geom_point() + geom_smooth(method = "lm")

#ok, what does rock_length correlate with?

#slope is unsurprising, rock walls
ggplot(data = rock, aes(x = rock_length, y = slope)) +
  geom_point() + geom_smooth(method = "lm")

#transect length is also unsurprising, walls
ggplot(data = rock, aes(x = rock_length, y = transect_length)) +
  geom_point() + geom_smooth(method = "lm")
#model time? Among----


#abundance 
#this model has deviation
among_density_model <- glmmTMB(density ~ scale(slope) + cluster + type +
                                  mussel_score + debris_score + (1|pair),
                                 family = gaussian(), data = among_site)
summary(among_density_model)
simulateResiduals(among_density_model, plot = TRUE)

#ggpredict

predict_among_model <- ggpredict(among_abundance_model, terms = "type")


ggplot(data = predict_among_model) +
  geom_pointrange(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high))

#size

among_size_model <- glmmTMB(log(average_size) ~ slope + cluster + type +
                              mussel_score + debris_score + (1|pair),
                            family = gaussian(), data = among_site)
summary(among_size_model)
simulateResiduals(among_size_model, plot = TRUE)


#richness

among_richness_model <- glmmTMB(richness ~ slope + cluster + type + transect_length + (1|pair), data = among_site)
summary(among_richness_model)
simulateResiduals(among_richness_model, plot = TRUE)





#graph time for fun

ggplot(data = within, aes(x = type, y = species_abundance, colour = species, group = species)) + 
  geom_jitter(position = position_jitterdodge(
  jitter.width = 0.1,
  dodge.width = 0.3),
  alpha = 0.4,
  size = 2) +
  stat_summary(fun = mean,
               geom = "point",
               position = position_dodge(width = 0.3),
               size = 3) +
  stat_summary(fun.data = mean_cl_normal,
               geom = "errorbar",
               position = position_dodge(width = 0.3),
               width = 0.2) +
  facet_wrap_paginate(~species) +
  theme_classic()

#graphs are actually a great way of cleaning up data

ggplot(data = within, aes(x = type, y = size, colour = species, group = species)) + 
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.1,
    dodge.width = 0.3),
    alpha = 0.4,
    size = 2) +
  stat_summary(fun = mean,
               geom = "point",
               position = position_dodge(width = 0.3),
               size = 3) +
  stat_summary(fun.data = mean_cl_normal,
               geom = "errorbar",
               position = position_dodge(width = 0.3),
               width = 0.2) +
  facet_wrap_paginate(~species) +
  theme_classic()
