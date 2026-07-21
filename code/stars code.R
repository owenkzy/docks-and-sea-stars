#code for honours project, sorry Claire, it won't be as clean as you'd like

#because of renv will need to redownload some packages


#load in packages----
library(tidyverse)
library(ggforce)
library(glmmTMB)
library(DHARMa)

#load in data----

#belt data + wrangling

belt <- read_csv("data/stars_belt - Sheet1.csv") %>%
  add_count(site, species, name = "abundance")

#get belt site data

belt_site <- read_csv("data/stars_belt - Sheet1.csv") %>%
  group_by(site, date) %>%
  mutate(average_size = mean(size, na.rm = TRUE)) %>%
  group_by(site, date, species, average_size) %>%
  summarise(species_abundance = n(), .groups = "drop") %>%
  group_by(site, date) %>%
  mutate(total_abundance = sum(species_abundance, na.rm = TRUE)) %>%
  ungroup()
                                             

#stars rove

stars_rove <- read_csv("data/stars_rove - Sheet1.csv")

#debris data


debris <- read_csv("data/debris_rove - Sheet1 (1).csv") %>%
  mutate(debris_score = as.factor(debris_score),
         mussel_score = as.factor(mussel_score))

#site data

site <- read_csv("data/stars_site - Sheet1 (1).csv") %>%
  mutate(type = factor(type, levels = c("dock", "inter-dock", "no-dock")),
         depth_change = end_depth - start_depth,
         slope = asin(depth_change/transect_length)*180/pi,
         cluster = as.factor(cluster),
         pair = as.factor(pair),
         site = as.factor(site))



#realistically, the data should exist at two levels, within site and at site


#leftjoin/fulljoin time----

#within sites df

within <- belt %>%
  full_join(site, by = c("site","date"), relationship = "many-to-one")

#among sites df

among <- site %>%
  full_join(stars_rove, by = c("site","date"), relationship = "many-to-many")

among <- among %>%
  full_join(belt_site, by = c("site", "date"), relationship = "many-to-many")

among <- among %>%
  full_join(debris, by = c("site","date"), relationship = "many-to-many")

stop

#model time? Among----

among_model <- glmmTMB(total_abundance ~ slope + cluster + type + debris_score + mussel_score + (1|pair), data = among)
summary(among_model)
simulateResiduals(among_model, plot = TRUE)








#graph time for fun

ggplot(data = within, aes(x = type, y = abundance, colour = species, group = species)) + 
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
