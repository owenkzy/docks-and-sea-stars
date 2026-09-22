#code for honours project, sorry Claire, it won't be as clean as you'd like

#because of renv will need to redownload some packages


#load in packages----
library(tidyverse)
library(ggforce)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(performance)
library(rstatix)

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

#belt_site_wide now has 0's for where species were not recorded

belt_site_wide <- belt_site %>%
  pivot_wider(names_from = species,
              values_from = species_abundance,
              values_fill = 0)

#this dataset now gives me the 0's needed for absence

belt_site_long <- belt_site_wide %>%
  pivot_longer(
    cols = c(leather, mottled, ochre, blood, pink, sunflower),
    names_to = "species",
    values_to = "species_abundance"
  ) %>%
  mutate(species = as.factor(species), )

#belt_site_long just incorporates that back into long data
#stars rove

stars_rove <- read_csv("data/stars_rove - Sheet1.csv")

#debris data


debris <- read_csv("data/debris_rove - Sheet1 (2).csv") %>%
  mutate(debris_score = as.numeric(as.character(debris_score)),
         mussel_score = as.numeric(as.character(mussel_score)))  %>%
  add_count(site, debris_type, name = "debris_count")


#site data

site <- read_csv("data/stars_site - Sheet1 (3).csv") %>%
  mutate(
    type = factor(type, levels = c("dock", "inter-dock", "no-dock")),
    depth_change = end_depth - start_depth,
    slope = asin(depth_change / transect_length) * 180 / pi,
    cluster = as.factor(cluster),
    pair = as.factor(pair),
    site = as.factor(site)
  )

#star photoquadrats

photo <- read_csv("data/stars_photoquadrat.csv") %>%
  group_by(site) %>%
  summarise(
    sand_n    = sum(sand, na.rm = TRUE),
    boulder_n = sum(boulder, na.rm = TRUE),
    cobble_n  = sum(cobble, na.rm = TRUE),
    mussel_n  = sum(mussel, na.rm = TRUE),
    kelp_n = sum(kelp, na.rm = TRUE),
    other_n = sum(other, na.rm = TRUE),
    debris_n = sum(debris, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_points = sand_n + boulder_n + cobble_n + mussel_n,
    sand_percent = 100 * sand_n / total_points,
    boulder_percent = 100 * boulder_n / total_points,
    cobble_percent = 100 * cobble_n / total_points,
    mussel_percent = 100 * mussel_n / total_points,
    kelp_percent = 100 * kelp_n / total_points,
    debris_percent = 100 * debris_n / total_points,
    other_percent = 100 * other_n / total_points
  ) %>%
  mutate(rock_percent = cobble_percent + boulder_percent)


#realistically, the data should exist at two levels, within site and at site


#leftjoin/fulljoin time----

#within sites df (UNCOMPLETED YET)

within <- belt %>%
  full_join(site, by = c("site", "date"), relationship = "many-to-many")
#quick within need to find average size of star to say certain species are larger

within <- within %>%
  group_by(species, site) %>%
  summarise(mean_size = mean(size, na.rm = TRUE),
            n = sum(!is.na(size))) %>%
  mutate(species = as.factor(species))

size_data <- belt %>%
  filter(species %in% c("leather", "mottled", "ochre")) %>%
  mutate(species = factor(species))

size <- glmmTMB(size ~ species + (1 |
                                    site),
                family = Gamma(link = "log"),
                data = size_data)
simulateResiduals(size, plot = TRUE)
summary(size)

library(emmeans)

size_emm <- emmeans(size, ~ species)

pairs(size_emm, adjust = "tukey")
#among sites df----

among <- site %>%
  full_join(stars_rove,
            by = c("site", "date"),
            relationship = "many-to-many")

among <- among %>%
  full_join(belt_site,
            by = c("site", "date"),
            relationship = "many-to-many")

among <- among %>%
  full_join(debris,
            by = c("site", "date"),
            relationship = "many-to-many")

among <- among %>%
  left_join(photo, by = "site")

#get it so only 1 row per observation
among_site <- among %>% distinct(site, .keep_all = TRUE)

among_site <- among_site %>%
  mutate(
    density = total_abundance / (transect_length * 2),
    slope_scaled = as.numeric(scale(slope)),
    distance_scaled = as.numeric(scale(dist_from_dock_new)),
    rock_scaled = as.numeric(scale(rock_percent)),
    lat_scaled = as.numeric(scale(lat)),
    area_transect = transect_length * 2,
    mussel_scaled = as.numeric(scale(mussel_percent))
  ) %>%
  distinct()

among_site$type <- relevel(factor(among_site$type), ref = "no-dock")

#trying it again so I have it with a dataset that includes 0's IS NOT ONLY 30 ROWS

trial_zero <- belt_site_long %>%
  full_join(site, by = c("site", "date"), relationship = "many-to-many")

trial_zero <- trial_zero %>%
  full_join(stars_rove,
            by = c("site", "date"),
            relationship = "many-to-many")

trial_zero <- trial_zero %>%
  full_join(debris,
            by = c("site", "date"),
            relationship = "many-to-many")

#get only one row per site

among_zero <- trial_zero %>%
  select(!average_size) %>%
  distinct()


#dock distance info
among_site %>%
  group_by(type) %>%
  summarise(
    n = n(),
    mean_distance = mean(dist_from_dock_new, na.rm = TRUE),
    sd_distance = sd(dist_from_dock_new, na.rm = TRUE),
    min_distance = min(dist_from_dock_new, na.rm = TRUE),
    max_distance = max(dist_from_dock_new, na.rm = TRUE)
  )

#site info
among_site %>%
  summarise(
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE),
    min_slope = min(slope, na.rm = TRUE),
    max_slope = max(slope, na.rm = TRUE),
    mean_rock = mean(rock_percent, na.rm = TRUE),
    sd_rock = sd(rock_percent, na.rm = TRUE),
    min_rock = min(rock_percent, na.rm = TRUE),
    max_rock = max(rock_percent, na.rm = TRUE),
    mean_sand = mean(sand_percent, na.rm = TRUE),
    sd_sand = sd(sand_percent, na.rm = TRUE),
    min_sand = min(sand_percent, na.rm = TRUE),
    max_sand = max(sand_percent, na.rm = TRUE))

#what's correlated with one another, will impact the model because of iid----

#cramers v is great for categorical stuff
#near 0 = little association
#around 0.1 = weak
#around 0.3 = moderate
#around 0.5 or higher = strong

cramer_v(table(among_site$type, among_site$debris_score))
#type and debris have strong correlation - 0.716
cramer_v(table(among_site$type, among_site$mussel_score))
#type and mussel have strong, but not as strong with debris, correlation - 0.584
cramer_v(table(among_site$mussel_score, among_site$debris_score))
#mussel score and debris have strong correlation - 0.587
cramer_v(table(among_site$cluster, among_site$type))
#type and cluster have weak

among_site <- among_site %>%
  mutate(debris_score = as.factor(debris_score),
         mussel_score = as.factor(mussel_score))

#fishers exact test, have to do it pairwise

library(effectsize)

#debris and mussel
debris_mussel <- table(among_site$debris_score, among_site$mussel_score)

set.seed(123)

fisher_debris_mussel <- fisher.test(debris_mussel, simulate.p.value = TRUE, B = 99999)

fisher_debris_mussel

#type with mussel

type_mussel <- table(among_site$type, among_site$mussel_score)

type_mussel

set.seed(123)

fisher_type_mussel <- fisher.test(type_mussel, simulate.p.value = TRUE, B = 99999)

fisher_type_mussel

cramers_v(type_mussel)

#type with debris

type_debris <- table(among_site$type, among_site$debris_score)

type_debris

set.seed(123)

fisher_type_debris <- fisher.test(type_debris, simulate.p.value = TRUE, B = 99999)

fisher_type_debris

cramers_v(type_debris)

#AIC/Differing models time DENSITY: ALL USING among_site df OFFSET----

#null model
null <- glmmTMB(
  total_abundance ~ 1 + (1 |
                           pair) + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(null)
simulateResiduals(null, plot = TRUE)


#spatial model including offset

spatial <- glmmTMB(
  total_abundance ~ cluster + offset(log(area_transect)) + (1 | pair),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(spatial)
simulateResiduals(spatial, plot = TRUE)
check_collinearity(spatial)

#latitude model

lat <- glmmTMB(
  total_abundance ~ lat_scaled + (1 | pair) +
    offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(lat)
simulateResiduals(lat, plot = TRUE)


#habitat model

habitat <- glmmTMB(
  total_abundance ~ slope_scaled + rock_scaled + cluster + (1 | pair) +
    offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(habitat)
simulateResiduals(habitat, plot = TRUE)
check_collinearity(habitat)

#dock model

dock_debris <- glmmTMB(
  total_abundance ~ debris_score + (1 |
                                      pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(dock)
simulateResiduals(dock, plot = TRUE)
check_collinearity(dock)

dock_mussel <- glmmTMB(
  total_abundance ~ mussel_score + (1 |
                                      pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
dock <- glmmTMB(
  total_abundance ~ type + (1 |
                              pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
dock_distance <- glmmTMB(
  total_abundance ~ distance_scaled + (1 |
                                         pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)

summary(distance)
simulateResiduals(distance, plot = TRUE)
#full model

full <- glmmTMB(
  total_abundance ~ type + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
full_debris <- glmmTMB(
  total_abundance ~ debris_score + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(full_debris)
simulateResiduals(full_debris, plot = TRUE)
full_mussel <- glmmTMB(
  total_abundance ~ mussel_score + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(full_mussel)

full_distance <- glmmTMB(
  total_abundance ~ distance_scaled + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(full)
simulateResiduals(full, plot = TRUE)
check_collinearity(full)


#aic time

model_set <- list(
  null = null,
  spatial = spatial,
  dock = dock,
  habitat = habitat,
  full = full,
  lat = lat,
  full_mussel = full_mussel,
  full_debris = full_debris,
  dock_mussel = dock_mussel,
  dock_debris = dock_debris,
  full_distance = full_distance,
  dock_distance = dock_distance
)

model_selection <- MuMIn::model.sel(model_set, rank = "AICc")

model_selection

df <- as.data.frame(model_selection)

#r2 stuff

r2_values <- sapply(model_set, function(x) {
  r <- suppressWarnings(performance::r2_nakagawa(x, tolerance = 1e-15, verbose = FALSE))
  
  as.numeric(unlist(r["R2_marginal"])[1])
})

round(r2_values, 3)

#models to keep

full_debris_number <- glmmTMB(
  total_abundance ~ debris_score + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(full_debris_number)
full_mussel_number <- glmmTMB(
  total_abundance ~ mussel_score + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(full_mussel_number)

#null binomial is better, so is spatial p and dockp, just gonna use nbinom then, getting rid of poisson code

#AIC COMPETETING SIZE TIME-----

#null model
null <- glmmTMB(average_size ~ 1 + (1 |pair),
                family = Gamma(link = "log"),
                data = among_site)
summary(null)
simulateResiduals(null, plot = TRUE)

#spatial model

spatial <- glmmTMB(average_size ~ cluster + (1 | pair),
                   family = Gamma(link = "log"),
                   data = among_site)
summary(spatial)
simulateResiduals(spatial, plot = TRUE)
check_collinearity(spatial)

#habitat model

habitat <- glmmTMB(
  average_size ~ slope_scaled + rock_scaled + (1 | pair),
  family = Gamma(link = "log"),
  data = among_site
)
summary(habitat)
simulateResiduals(habitat, plot = TRUE)
check_collinearity(habitat)


#latitude model

lat <- glmmTMB(average_size ~ lat_scaled + (1 | pair)
               ,
               family = Gamma(link = "log"),
               data = among_site)
summary(lat)
simulateResiduals(lat, plot = TRUE)


#dock model

dock <- glmmTMB(
  average_size ~ type + cluster + (1 |
                                     pair),
  family = Gamma(link = "log"),
  data = among_site
)
dock_mussel <- glmmTMB(
  average_size ~ mussel_score + cluster + (1 |
                                             pair),
  family = Gamma(link = "log"),
  data = among_site
)
dock_debris <- glmmTMB(
  average_size ~ debris_score + cluster + (1 |
                                             pair),
  family = Gamma(link = "log"),
  data = among_site
)
dock_distance <- glmmTMB(
  average_size ~ distance_scaled  +
    (1 |
       pair) + cluster,
  family = Gamma(link = "log"),
  data = among_site
)
summary(dock)
simulateResiduals(dock, plot = TRUE)
check_collinearity(dock)

#full model

full <- glmmTMB(
  average_size ~ cluster + type + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site
)

full_mussel <- glmmTMB(
  average_size ~ cluster + mussel_score + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site
)
full_debris <- glmmTMB(
  average_size ~ cluster + debris_score + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site
)
full_distance <- glmmTMB(
  average_size ~ distance_scaled + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster,
  family = Gamma(link = "log"),
  data = among_site
)
summary(full)
simulateResiduals(full, plot = TRUE)
check_collinearity(full)

#aic time

model_set <- list(
  null = null,
  spatial = spatial,
  dock = dock,
  habitat = habitat,
  full = full,
  lat = lat,
  full_mussel = full_mussel,
  full_debris = full_debris,
  dock_mussel = dock_mussel,
  dock_debris = dock_debris,
  full_distance = full_distance,
  dock_distance = dock_distance
  
)

model_selection <- MuMIn::model.sel(model_set, rank = "AICc")

model_selection

#r2 stuff

r2_values <- sapply(model_set, function(x) {
  r <- suppressWarnings(performance::r2_nakagawa(x, tolerance = 1e-15, verbose = FALSE))
  
  as.numeric(unlist(r["R2_marginal"])[1])
})

round(r2_values, 3)

#ggeffects time, needs the models saved so need a unique ID
dock_size <- glmmTMB(
  average_size ~ type + cluster + (1 |
                                     pair),
  family = Gamma(link = "log"),
  data = among_site
)
summary(dock_size)

full_dock_size <- glmmTMB(
  average_size ~ cluster + type + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site
)
summary(full_dock_size)
simulateResiduals(full_dock_size, plot = TRUE)


#AIC competing time RICHNESS----


#null model
null <- glmmTMB(richness ~ 1 + (1 |
                                  pair),
                family = Gamma(link = "log"),
                data = among_site)
summary(null)
simulateResiduals(null, plot = TRUE)

#spatial model

spatial <- glmmTMB(richness ~ cluster + (1 | pair),
                   family = Gamma(link = "log"),
                   data = among_site)
summary(spatial)
simulateResiduals(spatial, plot = TRUE)
check_collinearity(spatial)

#habitat model

habitat <- glmmTMB(
  richness ~ slope_scaled + rock_scaled + (1 | pair),
  family = Gamma(link = "log"),
  data = among_site
)
summary(habitat)
simulateResiduals(habitat, plot = TRUE)
check_collinearity(habitat)

#latitude model

lat <- glmmTMB(richness ~ lat_scaled + (1 | pair)
               ,
               family = Gamma(link = "log"),
               data = among_site)
summary(lat)
simulateResiduals(lat, plot = TRUE)


#dock model

dock <- glmmTMB(richness ~ type + cluster + (1 |
                                               pair),
                family = Gamma(link = "log"),
                data = among_site)
dock_mussel <- glmmTMB(
  richness ~ mussel_score + cluster + (1 |
                                         pair),
  family = Gamma(link = "log"),
  data = among_site
)
dock_debris <- glmmTMB(
  richness ~ debris_score + cluster + (1 |
                                         pair),
  family = Gamma(link = "log"),
  data = among_site
)

dock_distance <- glmmTMB(
  richness ~ distance_scaled  +
    (1 |
       pair) + cluster,
  family = Gamma(link = "log"),
  data = among_site
)
summary(dock)
simulateResiduals(dock, plot = TRUE)
check_collinearity(dock)

#full model

full <- glmmTMB(
  richness ~ cluster + type + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site
)

full_mussel <- glmmTMB(
  richness ~ cluster + mussel_score + slope_scaled + rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site
)
full_debris <- glmmTMB(
  richness ~ cluster + debris_score + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site
)
full_distance <- glmmTMB(
  richness ~ distance_scaled + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster,
  family = Gamma(link = "log"),
  data = among_site
)
summary(full)
simulateResiduals(full, plot = TRUE)
check_collinearity(full)

#aic time


model_set <- list(
  null = null,
  spatial = spatial,
  dock = dock,
  habitat = habitat,
  full = full,
  lat = lat,
  full_mussel = full_mussel,
  full_debris = full_debris,
  dock_mussel = dock_mussel,
  dock_debris = dock_debris,
  full_distance = full_distance,
  dock_distance = dock_distance
  
  
)
model_selection <- MuMIn::model.sel(model_set, rank = "AICc")

model_selection

#r2 stuff

r2_values <- sapply(model_set, function(x) {
  r <- suppressWarnings(performance::r2_nakagawa(x, tolerance = 1e-15, verbose = FALSE))
  
  as.numeric(unlist(r["R2_marginal"])[1])
})

round(r2_values, 3)

#null was the best


#something weird below, 2026-09-14, I forgot why I made this code below

null_rich <- glmmTMB(richness ~ 1 + (1 |
                                       pair),
                     family = Gamma(link = "log"),
                     data = among_site)
summary(null_rich)
simulateResiduals(null, plot = TRUE)

#NMDS time ----

#a piece of advice given is to sqrt some of the uber high densities

#getting dataframe of species

community <- belt_site_wide %>%
  select(c(site, leather, mottled, ochre, blood, sunflower, pink)) %>%
  as.data.frame()

#turn row names into site???

rownames(community) <- community$site
community$site <- NULL

#sqrt values of species abundance

community <- sqrt(community)

#vegan time?
library(vegan)

community_dist <- vegdist(community, method = "bray")

community_dist

#nmds time

set.seed(123)

nmds <- metaMDS(
  community,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE
)
#stress is 0,105 which is good! Can infer stuff
#let's make a stressplot

stressplot(
  nmds,
  pch = 19,
  p.col = "blue",
  l.col = "red",
  lwd = 2
)


print(nmds)



nmds_scores <- as.data.frame(scores(nmds, display = "sites"))

nmds_scores$site <- rownames(nmds_scores)

head(nmds_scores)

#join site type back

nmds_scores <- nmds_scores %>%
  left_join(site %>%
              select(site, type, pair), by = "site")

head(nmds_scores)
table(nmds_scores$type, useNA = "ifany")

#nmds plot???

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = type)) +
  geom_point(size = 3) +
  coord_equal() +
  theme_classic()

#trying this plot from online cause it seems cool----

species_scores <- as.data.frame(scores(nmds, display = "species"))

species_scores$species <- rownames(species_scores)

library(ggalt)

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = type)) +
  coord_equal() +
  geom_hline(yintercept = 0,
             color = 'lightgray',
             alpha = .4) +
  geom_vline(xintercept = 0,
             color = 'lightgray',
             alpha = .4) +
  geom_text(inherit.aes = FALSE,
            data = species_scores,
            aes(x = NMDS1, y = NMDS2, label = species)) +
  geom_point(aes(colour = type), size = 2) +
  geom_encircle(
    aes(colour = type, fill = type),
    s_shape = 1,
    alpha = 0.2,
    size = 1,
    expand = 0
  ) +
  scale_x_continuous(expand = c(.01, .05)) +
  theme(legend.position = 'bottom', panel.grid = element_blank())


library(ggrepel)

ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
  geom_mark_hull(
    aes(colour = type, fill = type),
    alpha = 0.10,
    expand = unit(2, "mm"),
    concavity = 4,
    show.legend = FALSE) +
  geom_hline(yintercept = 0,
             colour = "grey85",
             linewidth = 0.4) +
  geom_vline(xintercept = 0,
             colour = "grey85",
             linewidth = 0.4) +
  geom_point(aes(colour = type, shape = type),
             size = 3,
             alpha = 0.85) +
  geom_text_repel(
    data = species_scores,
    aes(NMDS1, NMDS2, label = species),
    inherit.aes = FALSE,
    size = 4,
    fontface = "italic",
    colour = "grey15",
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.colour = "grey65",
    max.overlaps = Inf) + 
  scale_colour_manual(
    values = c(
      "dock"       = "#D55E00",
      "inter-dock" = "#009E73",
      "no-dock"    = "#0072B2"
    ),
    labels = c(
      "dock"       = "Dock",
      "inter-dock" = "Inter-dock",
      "no-dock"    = "No-dock"
    )
  ) +
  scale_fill_manual(values = c(
    "dock"       = "#D55E00",
    "inter-dock" = "#009E73",
    "no-dock"    = "#0072B2"
  )) +
  scale_shape_manual(values = c(
    "dock"       = 16,
    "inter-dock" = 17,
    "no-dock"    = 15
  )) +
  coord_equal(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    colour = "Site type",
    shape = "Site type",
    subtitle = paste0("Stress = ", round(nmds$stress, 3))
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 20, 10, 20))

ggsave(
  "figures/nmds.png",
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

# Fixed limits so points remain in the same positions between slides
x_limits <- range(c(nmds_scores$NMDS1, species_scores$NMDS1), na.rm = TRUE)

y_limits <- range(c(nmds_scores$NMDS2, species_scores$NMDS2), na.rm = TRUE)

# Function for making each version
make_nmds_plot <- function(groups, show_species = FALSE) {
  plot_data <- nmds_scores %>%
    filter(type %in% groups)
  
  p <- ggplot(plot_data, aes(NMDS1, NMDS2)) +
    geom_hline(
      yintercept = 0,
      colour = "grey85",
      linewidth = 0.4) +
    geom_vline(
      xintercept = 0,
      colour = "grey85",
      linewidth = 0.4) +
    geom_encircle(
      aes(colour = type, fill = type),
      s_shape = 1,
      alpha = 0.12,
      linewidth = 0.8,
      expand = 0,
      show.legend = FALSE) +
    geom_point(aes(colour = type, shape = type),
               size = 3,
               alpha = 0.9) +
    scale_colour_manual(
      values = c(
        "dock"       = "#D55E00",
        "inter-dock" = "#009E73",
        "no-dock"    = "#0072B2"),
      labels = c(
        "dock"       = "Dock",
        "inter-dock" = "Inter-dock",
        "no-dock"    = "No-dock"),
      drop = FALSE) +
    scale_fill_manual(
      values = c(
        "dock"       = "#D55E00",
        "inter-dock" = "#009E73",
        "no-dock"    = "#0072B2"),
      drop = FALSE) +
    scale_shape_manual(
      values = c(
        "dock"       = 16,
        "inter-dock" = 17,
        "no-dock"    = 15),
      labels = c(
        "dock"       = "Dock",
        "inter-dock" = "Inter-dock",
        "no-dock"    = "No-dock"),
      drop = FALSE) +
    coord_equal(
      xlim = x_limits + c(-0.15, 0.15),
      ylim = y_limits + c(-0.15, 0.15),
      clip = "off") +
    labs(
      x = "NMDS1",
      y = "NMDS2",
      colour = "Site type",
      shape = "Site type",
      subtitle = paste0("Stress = ", round(nmds$stress, 3))) +
    ggpubr::theme_pubr(base_size = 13) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 20, 10, 15))
  
#Add species labels only when requested
  if (show_species) {
    p <- p +
      geom_text_repel(
        data = species_scores,
        aes(NMDS1, NMDS2, label = species),
        inherit.aes = FALSE,
        size = 4,
        fontface = "italic",
        colour = "grey15",
        box.padding = 0.5,
        point.padding = 0.3,
        min.segment.length = 0,
        segment.colour = "grey65",
        max.overlaps = Inf
      )
  }
  
  p
}

# Slide 1: dock sites
nmds_1 <- make_nmds_plot(groups = "dock")

# Slide 2: add inter-dock sites
nmds_2 <- make_nmds_plot(groups = c("dock", "inter-dock"))

# Slide 3: add no-dock sites and species labels
nmds_3 <- make_nmds_plot(groups = c("dock", "inter-dock", "no-dock"))

nmds_1
nmds_2
nmds_3

ggsave(
  "figures/nmds_1_dock.png",
  nmds_1,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/nmds_2_interdock.png",
  nmds_2,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/nmds_3_all_sites.png",
  nmds_3,
  width = 8,
  height = 6,
  dpi = 300
)
#permanova time?

community_metadata <- nmds_scores %>%
  select(site, type, pair)

rownames(community_metadata) <- community_metadata$site
community_metadata$site <- NULL

identical(rownames(community_metadata), rownames(community))

permanova <- adonis2(
  community ~ type,
  data = community_metadata,
  method = "bray",
  permutations = 9999
)

permanova

#so groups do differ, but not sure which one

#checking if spread inside each group is different
dispersion <- betadisper(vegdist(community, method = "bray"),
                         group = community_metadata$type,
                         sqrt.dist = TRUE)

permutest(dispersion, permutations = 9999)

#the spread within each group is not different

#pairwise time, is the community at dock different than interdock?----

#filtering
dvid <- community_metadata$type %in%
  c("dock", "inter-dock")

#raw data with the bray curtis numbers
community_dvid <- community[dvid, ]

#explanatory data paired with the bray curtis
metadata_dvid <- community_metadata[dvid, , drop = FALSE]

metadata_dvid$type <- droplevels(metadata_dvid$type)
metadata_dvid$pair <- droplevels(metadata_dvid$pair)

table(metadata_dvid$pair, metadata_dvid$type)
#run permanova time

set.seed(123)

permanova_dvid <- adonis2(
  community_dvid ~ type,
  data = metadata_dvid,
  method = "bray",
  permutations = 9999,
  strata = metadata_dvid$pair
)
permanova_dvid

#is the community at dock different than no-dock?

dvnd <- community_metadata$type %in%
  c("dock", "no-dock")

community_dvnd <- community[dvnd, , drop = FALSE]

metadata_dvnd <- community_metadata[dvnd, , drop = FALSE]

metadata_dvnd$type <- droplevels(metadata_dvnd$type)

table(metadata_dvnd$type)

#model time

set.seed(123)

permanova_dvnd <- adonis2(
  community_dvnd ~ type,
  data = metadata_dvnd,
  method = "bray",
  permutations = 9999
)

permanova_dvnd

#is the community at interdock different than no-dock?

idvnd <- community_metadata$type %in%
  c("inter-dock", "no-dock")

community_idvnd <- community[idvnd, , drop = FALSE]

metadata_idvnd <- community_metadata[idvnd, , drop = FALSE]

metadata_idvnd$type <- droplevels(metadata_idvnd$type)

table(metadata_idvnd$type)
#model time
set.seed(123)

permanova_idvnd <- adonis2(
  community_idvnd ~ type,
  data = metadata_idvnd,
  method = "bray",
  permutations = 9999
)

permanova_idvnd

#holm correction for doing it 3 times

pairwise_results <- data.frame(
  comparison = c("dock vs inter-dock", " dock vs no-dock", "inter-dock vs no-dock"),
  
  F_value = c(permanova_dvid$F[1], permanova_dvnd$F[1], permanova_idvnd$F[1]),
  R2 = c(permanova_dvid$R2[1], permanova_dvnd$R2[1], permanova_idvnd$R2[1]),
  p_unadjusted = c(
    permanova_dvid$`Pr(>F)`[1],
    permanova_dvnd$`Pr(>F)`[1],
    permanova_idvnd$`Pr(>F)`[1]
  )
)

pairwise_results$p_holm <- p.adjust(pairwise_results$p_unadjusted, method = "holm")

pairwise_results

#envfit

library(vegan)

species_fit <- envfit(nmds, belt_site_wide[, c("leather", "mottled", "ochre", "blood", "sunflower", "pink")], permutations = 999)

species_fit

#vector coordinates
vec <- as.data.frame(vegan::scores(species_fit, display = "vectors"))

#species name
vec$species <- rownames(vec)

#p values
vec$p <- species_fit$vectors$pvals

vec
vec_sig <- vec |>
  dplyr::filter(p < 0.05)

vec_sig
arrow_scale <- 0.65

library(ggrepel)

nmds_3 +
  geom_segment(
    data = vec_sig,
    aes(
      x = 0,
      y = 0,
      xend = NMDS1 * arrow_scale,
      yend = NMDS2 * arrow_scale),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.12, "cm")),
    linewidth = 0.4) +
  geom_text_repel(
    data = vec_sig,
    aes(
      x = NMDS1 * arrow_scale * 1.05,
      y = NMDS2 * arrow_scale * 1.05,
      label = species),
    inherit.aes = FALSE,
    size = 4,
    box.padding = 0.25,
    point.padding = 0.15,
    segment.color = NA) +
  coord_fixed(ratio = 1)

ggsave(
  "figures/multivariate.png",
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

vec_sig <- vec_sig |>
  mutate(
    label_x = case_when(
      species == "mottled" ~ 0.39,
      species == "ochre"   ~ 0.40,
      species == "leather" ~ 0.40
    ),
    label_y = case_when(
      species == "mottled" ~ 0.41,
      species == "ochre"   ~ 0.28,
      species == "leather" ~ -0.42
    )
  )

#Final NMDS plot
nmds_final <- nmds_3 +
geom_segment(
    data = vec_sig,
    aes(
      x = 0,
      y = 0,
      xend = NMDS1 * arrow_scale,
      yend = NMDS2 * arrow_scale
    ),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.12, "cm")),
    linewidth = 0.4
  ) +
  geom_text(
    data = vec_sig,
    aes(x = label_x, y = label_y, label = species),
    inherit.aes = FALSE,
    size = 4
  ) +
  coord_fixed(ratio = 1) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
        legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    legend.margin = margin(t = -5),
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 2,
      l = 5)) +
  guides(
    colour = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1),
    fill = guide_legend(nrow = 1)
  )

nmds_final
ggsave(
  "figures/NMDS_final.png",
  nmds_final,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)

#let's do some pycno stuff----

sun <- read_csv("data/pycno_data - Sheet1.csv") %>%
  group_by(site) %>%
  mutate(
    average_pycno_size = mean(pycno_size, na.rm = TRUE),
    average_pycno_size = replace(average_pycno_size, is.nan(average_pycno_size), NA_real_)
  ) %>%
  ungroup() %>%
  distinct(site, pycno_count, average_pycno_size)

#leftjoin it to full data

sun <- sun %>%
  left_join(among_site, by = "site")


sun <- sun %>%
  mutate(pycno_density = pycno_count / (transect_length * 6),
         transect_rove_area = transect_length * 2)

#model it AICC Pycno time using offset for count

#null model
null <- glmmTMB(
  pycno_count ~ 1 + (1 |
                       pair) + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)
summary(null)
simulateResiduals(null, plot = TRUE)


#spatial model including offset

spatial <- glmmTMB(
  pycno_count ~ cluster + offset(log(area_transect)) + (1 | pair),
  family = nbinom2(link = "log"),
  data = sun
)
summary(spatial)
simulateResiduals(spatial, plot = TRUE)
check_collinearity(spatial)

#latitude model

lat <- glmmTMB(
  pycno_count ~ lat_scaled + (1 | pair) +
    offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)
summary(lat)
simulateResiduals(lat, plot = TRUE)


#habitat model

habitat <- glmmTMB(
  pycno_count ~ slope_scaled + rock_scaled + cluster + (1 | pair) +
    offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)
summary(habitat)
simulateResiduals(habitat, plot = TRUE)
check_collinearity(habitat)

#dock model

dock_debris <- glmmTMB(
  pycno_count ~ debris_score + (1 |
                                  pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)
summary(dock)
simulateResiduals(dock, plot = TRUE)
check_collinearity(dock)

dock_mussel <- glmmTMB(
  pycno_count ~ mussel_score + (1 |
                                  pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)
dock <- glmmTMB(
  pycno_count ~ type + (1 |
                          pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)
dock_distance <- glmmTMB(
  pycno_count ~ distance_scaled + (1 |
                                     pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)

#full model

full <- glmmTMB(
  pycno_count ~ type + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)
full_debris <- glmmTMB(
  pycno_count ~ debris_score + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)
summary(full_debris)
simulateResiduals(full_debris, plot = TRUE)
full_mussel <- glmmTMB(
  pycno_count ~ mussel_score + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)
summary(full_mussel)

full_distance <- glmmTMB(
  pycno_count ~ distance_scaled + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)
summary(full)
simulateResiduals(full, plot = TRUE)
check_collinearity(full)


#aic time

model_set <- list(
  null = null,
  spatial = spatial,
  dock = dock,
  habitat = habitat,
  full = full,
  lat = lat,
  full_mussel = full_mussel,
  full_debris = full_debris,
  dock_mussel = dock_mussel,
  dock_debris = dock_debris,
  full_distance = full_distance,
  dock_distance = dock_distance
)

model_selection <- MuMIn::model.sel(model_set, rank = "AICc")

model_selection

r2_results <- lapply(model_set, r2)

r2_results


r2_results <- lapply(model_set, performance::r2, tolerance = 1e-12)

r2_results
data.frame(
  model = names(model_set),
  n = sapply(model_set, nobs),
  K = sapply(model_set, \(x) attr(logLik(x), "df")),
  convergence = sapply(model_set, \(x) x$fit$convergence),
  pdHess = sapply(model_set, \(x) x$sdr$pdHess)
)

#debris is just not happy, removing debris related models


model_set <- list(
  null = null,
  spatial = spatial,
  dock = dock,
  habitat = habitat,
  full = full,
  lat = lat,
  full_mussel = full_mussel,
  dock_mussel = dock_mussel,
  full_distance = full_distance,
  dock_distance = dock_distance
)

model_selection <- MuMIn::model.sel(model_set, rank = "AICc")

model_selection

#r2 stuff

r2_values <- sapply(model_set, function(x) {
  r <- suppressWarnings(performance::r2_nakagawa(x, tolerance = 1e-15, verbose = FALSE))
  
  as.numeric(unlist(r["R2_marginal"])[1])
})

round(r2_values, 3)

#justifying variable transect length----

#quickly do the lm for transect length

a <- lm(total_abundance ~ transect_length, data = among_site)
summary(a)

ggplot(data = among_site, aes(x = transect_length, y = total_abundance)) + geom_point()

#good to have as a back up

ggplot(data = among_site, (aes(x = transect_length, y = total_abundance))) +
  geom_point() + geom_smooth(method = "lm")

simulateResiduals(a, plot = TRUE)

#ok let's do more work to justify it, what substrate are the stars hanging out on ----

belt_j <- read_csv("data/stars_belt - Sheet1 (1).csv") %>%
  drop_na(site, species, substrate) %>%
  count(site, date, substrate, name = "substrate_abundance") %>%
  group_by(site, date) %>%
  mutate(
    total_abundance = sum(substrate_abundance),
    substrate_percent = substrate_abundance / total_abundance * 100
  ) %>%
  ungroup()

#summarize some substrate

belt_j <- belt_j %>%
  mutate(
    substrate = case_when(
      substrate %in% c("sand/boulder", "boulder/sand") ~
        "mixed",
      substrate %in% c("cobble/boulder", "boulder/cobble") ~
        "mixed",
      substrate %in% c("sand/cobble", "cobble/sand") ~
        "mixed",
      TRUE ~ substrate
    )
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

#pivot it back to normal + percent

belt_j_long <- belt_j_wider %>%
  pivot_longer(
    cols = c(boulder, sand, debris, cobble, mixed),
    names_to = "substrate",
    values_to = "substrate_abundance"
  ) %>%
  group_by(site, date) %>%
  mutate(
    total_abundance = sum(substrate_abundance, na.rm = TRUE),
    substrate_percent = if_else(
      total_abundance > 0,
      substrate_abundance / total_abundance * 100,
      NA_real_
    ),
    substrate = factor(substrate)
  ) %>%
  ungroup()


#without care for species

ggplot(data = belt_j_long, aes(x = substrate, y = substrate_percent)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.3),
    alpha = 0.4,
    size = 2
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.3),
    size = 3
  ) +
  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    position = position_dodge(width = 0.3),
    width = 0.2
  ) +
  theme_classic()


#what substrate are specific species of stars hanging out on?----

belt_j_species <- read_csv("data/stars_belt - Sheet1 (1).csv") %>%
  drop_na(site, date, species, substrate) %>%
  filter(!species %in% c("sunflower", "pink", "blood")) %>%
  count(site, date, species, substrate, name = "substrate_abundance") %>%
  group_by(site, date, species) %>%
  complete(substrate = substrate_levels,
           fill = list(substrate_abundance = 0)) %>%
  mutate(
    species_total = sum(substrate_abundance),
    substrate_percent = substrate_abundance / species_total * 100
  ) %>%
  ungroup()


belt_j_species <- belt_j_species %>%
  mutate(
    substrate = case_when(
      substrate %in% c("sand/boulder", "boulder/sand") ~
        "mixed",
      substrate %in% c("cobble/boulder", "boulder/cobble") ~
        "mixed",
      substrate %in% c("sand/cobble", "cobble/sand") ~
        "mixed",
      TRUE ~ substrate
    )
  )

#including species

ggplot(data = belt_j_species,
       aes(
         x = substrate,
         y = substrate_percent,
         group = species,
         colour = species
       )) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.3),
    alpha = 0.4,
    size = 2
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.3),
    size = 3
  ) +
  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    position = position_dodge(width = 0.3),
    width = 0.2
  ) +
  facet_wrap_paginate( ~ species) +
  theme_classic() + theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  ))


#ok last justification part, try the density by rocky area----

rock_data <- belt %>%
  filter(rock_keep == 1) %>%
  add_count(site, species, name = "rock_abundance")

rock <- rock_data %>%
  full_join(site, by = c("site", "date"), relationship = "many-to-many") %>%
  mutate(rock_density = rock_abundance / (rock_length * 2))

rock <- rock_data %>%
  left_join(photo, by = site)

ggplot(data = rock, aes(x = rock_length, y = rock_density)) +
  geom_point() + geom_smooth(method = "lm")

#ok, what does rock_length correlate with?

#slope is unsurprising, rock walls
ggplot(data = rock, aes(x = rock_length, y = slope)) +
  geom_point() + geom_smooth(method = "lm")

#transect length is also unsurprising, walls
ggplot(data = rock, aes(x = rock_length, y = transect_length)) +
  geom_point() + geom_smooth(method = "lm")

#rock percent and stuff but idrk what claire wants

ggplot(data = rock, aes(x = rock_percent))

#presentation ready graphics ----

pub_theme1 <- function() {
  theme_classic() +
    theme(
      # axis
      axis.text = element_text(size = 11, colour = "black"),
      axis.title = element_text(size = 11, colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.1),
      axis.ticks = element_line(color = "black", linewidth = 0.2),
      axis.ticks.length = unit(0.1, "lines"),
      # legend
      legend.title = element_text(size = 11, colour = "black"),
      legend.text = element_text(size = 11, colour = "black"),
      # panel border
      panel.background = element_rect(fill = "white", color  =  NA),
      panel.border = element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.5
      ),
      panel.spacing = unit(0.1, "lines")
    )
}

pub_theme1 <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = element_text(size = 11, colour = "black"),
      axis.title = element_text(size = 11, colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.1),
      axis.ticks = element_line(colour = "black", linewidth = 0.2),
      axis.ticks.length = grid::unit(0.1, "lines"),
      
      legend.title = element_text(size = 11, colour = "black"),
      legend.text = element_text(size = 11, colour = "black"),
      
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.5
      ),
      panel.spacing = grid::unit(0.1, "lines")
    )
}

library(patchwork)

#first thing first, mussel score and debris score are confounding----

score_plot <- among_site %>%
  select(site, type, debris_score, mussel_score) %>%
  pivot_longer(
    cols = c(debris_score, mussel_score),
    names_to = "variable",
    values_to = "score"
  ) %>%
  ggplot(aes(x = type, y = score)) +
  geom_jitter(
    width = 0.12,
    height = 0.05,
    size = 3,
    alpha = 0.75
  ) +
  facet_wrap(~ variable, labeller = as_labeller(c(
    debris_score = "Debris score", mussel_score = "Mussel score"
  ))) +
  scale_y_continuous(breaks = 0:4) +
  labs(x = "Site type", y = "Score") +
  pub_theme1()

ggsave(
  "figures/scores.png",
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

print(score_plot)

#dock and full_dock_size is the best model for size----
summary(dock_size)
dock_pred <- ggeffects::ggpredict(dock_size, terms = "type")

dock_size_plot <- ggplot() +
  geom_jitter(
    data = among_site,
    aes(x = type, y = average_size),
    width = 0.12,
    height = 0,
    size = 2.5,
    alpha = 0.35
  ) +
  geom_errorbar(
    data = dock_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    width = 0.12,
    linewidth = 0.9
  ) +
  geom_point(data = dock_pred, aes(x = x, y = predicted), size = 4) +
  labs(x = "Site type", y = expression("Average sea star size (cm)")) +
  pub_theme1()
print(dock_size_plot)
ggsave(
  "figures/docksize.png",
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

#now for the full model
summary(full_dock_size)
dock_pred <- ggeffects::ggpredict(full_dock_size, terms = "type")

full_size_plot <- ggplot() +
  geom_jitter(
    data = among_site,
    aes(x = type, y = average_size),
    width = 0.12,
    height = 0,
    size = 2.5,
    alpha = 0.35
  ) +
  geom_errorbar(
    data = dock_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    width = 0.12,
    linewidth = 0.9
  ) +
  geom_point(data = dock_pred, aes(x = x, y = predicted), size = 4) +
  labs(x = "Site type", y = expression("Average sea star size (cm)")) +
  pub_theme1()
print(full_size_plot)

#nice ok now to make them panel A and panel B

size_plot <- (dock_size_plot + full_size_plot + plot_annotation(tag_level = "A"))
ggsave(
  "figures/size.png",
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)


simulateResiduals(dock_size, plot = TRUE)
simulateResiduals(full_dock_size, plot = TRUE)

#COVARIATES FOR SIZE NOW----

#cluster for dock size

dock_pred <- ggeffects::ggpredict(dock_size,
                                  terms = "cluster"
                                  )

dock_size_cluster_plot <- ggplot() +
  geom_jitter(
    data = among_site,
    aes(x = cluster, y = average_size),
    width = 0.12,
    height = 0,
    size = 2.5,
    alpha = 0.35) +
  geom_errorbar(
    data = dock_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    width = 0.12,
    linewidth = 0.9) +
  geom_point(data = dock_pred, aes(x = x, y = predicted), size = 4) +
  labs(x = "Cluster", y = expression("Average sea star size (cm)")) +
  pub_theme1()
print(dock_size_cluster_plot)


#rock scaled for full dock size
dock_pred <- ggeffects::ggpredict(full_dock_size,
                                  terms = "rock_scaled"
                                  )

full_dock_size_rock_plot <- ggplot() +
  geom_point(
    data = among_site,
    aes(x = rock_scaled, y = average_size),
    size = 2.5,
    alpha = 0.35) +
  geom_ribbon(
    data = dock_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.2) +
  geom_line(
    data = dock_pred,
    aes(x = x, y = predicted),
    linewidth = 1) +
  labs(
    x = "Scaled rock cover",
    y = expression("Average sea star size (cm)")) +
  pub_theme1()
print(full_dock_size_rock_plot)

#stitch size again----

size_plot_a <-
  full_size_plot/
  full_dock_size_rock_plot +
  plot_layout(
    heights = c(1.5, 1),
    axes = "collect_y",
    axis_titles = "collect_y"
  )  + plot_annotation(tag_levels = "A")
print(size_plot_a)
ggsave(
  "figures/size_fulldock.png",
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)


size_plot_b <-
  dock_size_plot/
  dock_size_cluster_plot +
  plot_layout(
    heights = c(1.5, 1),
    axes = "collect_y",
    axis_titles = "collect_y"
  )  + plot_annotation(tag_levels = "A")
print(size_plot_b)
ggsave(
  "figures/size_dock.png",
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)
#abundance time----
#doing this in density, again...

among_site <- among_site %>%
  mutate(star_density = total_abundance / (transect_length * 2))

#graph for debris

dock_pred <- ggeffects::ggpredict(full_debris_number,
                                  terms = "debris_score",
                                  condition = c(area_transect = 1))

full_debris_number_plot <- ggplot() +
  geom_jitter(
    data = among_site,
    aes(x = debris_score, y = star_density),
    width = 0.12,
    height = 0,
    size = 2.5,
    alpha = 0.35) +
  geom_errorbar(
    data = dock_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    width = 0.12,
    linewidth = 0.9) +
  geom_point(data = dock_pred, aes(x = x, y = predicted), size = 4) +
  labs(x = "Debris score", y = expression("Sea star density (stars " * m ^
                                            {
                                              -2
                                            } * ")")) +
  pub_theme1()
print(full_debris_number_plot)

#graph for mussels


dock_pred <- ggeffects::ggpredict(full_mussel_number,
                                  terms = "mussel_score",
                                  condition = c(area_transect = 1))

full_mussel_number_plot <- ggplot() +
  geom_jitter(
    data = among_site,
    aes(x = mussel_score, y = star_density),
    width = 0.12,
    height = 0,
    size = 2.5,
    alpha = 0.35) +
  geom_errorbar(
    data = dock_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    width = 0.12,
    linewidth = 0.9) +
  geom_point(data = dock_pred, aes(x = x, y = predicted), size = 4) +
  labs(x = "Mussel score", y = expression("Sea star density (stars " * m ^
                                            {
                                              -2
                                            } * ")")) +
  pub_theme1()
print(full_mussel_number_plot)

#stitch them together

density_plot <- (full_debris_number_plot + full_mussel_number_plot + plot_annotation(tag_level = "A"))
print(density_plot)
ggsave(
  "figures/density.png",
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

summary(full_debris_number)
simulateResiduals(full_debris_number, plot = TRUE)
simulateResiduals(full_mussel_number, plot = TRUE)
summary(full_mussel_number)

#ok now we're gonna do the covariates for abundance MUSSEL----
#slope for full mussel number

dock_pred <- ggeffects::ggpredict(full_mussel_number,
                                  terms = "slope_scaled",
                                  condition = c(area_transect = 1))

full_mussel_number_slope_plot <- ggplot() +
  geom_point(
    data = among_site,
    aes(x = slope_scaled, y = star_density),
    size = 2.5,
    alpha = 0.35) +
  geom_ribbon(
    data = dock_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.2) +
  geom_line(
    data = dock_pred,
    aes(x = x, y = predicted),
    linewidth = 1) +
  labs(
    x = "Scaled slope",
    y = expression("Sea star density (stars " * m^{-2} * ")")) +
  pub_theme1()
full_mussel_number_slope_plot

#rock scaled for mussel number

dock_pred <- ggeffects::ggpredict(full_mussel_number,
                                  terms = "rock_scaled",
                                  condition = c(area_transect = 1))

full_mussel_number_rock_plot <- ggplot() +
  geom_point(
    data = among_site,
    aes(x = rock_scaled, y = star_density),
    size = 2.5,
    alpha = 0.35) +
  geom_ribbon(
    data = dock_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.2) +
  geom_line(
    data = dock_pred,
    aes(x = x, y = predicted),
    linewidth = 1) +
  labs(
    x = "Scaled rock cover",
    y = expression("Sea star density (stars " * m^{-2} * ")")) +
  pub_theme1()
full_mussel_number_rock_plot

#ok now we're gonna do the covariates for abundance debris----
#rock scaled

dock_pred <- ggeffects::ggpredict(full_debris_number,
                                  terms = "rock_scaled",
                                  condition = c(area_transect = 1))

full_debris_number_rock_plot <- ggplot() +
  geom_point(
    data = among_site,
    aes(x = rock_scaled, y = star_density),
    size = 2.5,
    alpha = 0.35) +
  geom_ribbon(
    data = dock_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.2) +
  geom_line(
    data = dock_pred,
    aes(x = x, y = predicted),
    linewidth = 1) +
  labs(
    x = "Scaled rock cover",
    y = expression("Sea star density (stars " * m^{-2} * ")")) +
  pub_theme1()
full_debris_number_rock_plot

#cluster for full debris number

dock_pred <- ggeffects::ggpredict(full_debris_number,
                                  terms = "cluster",
                                  condition = c(area_transect = 1))

full_debris_number_cluster_plot <- ggplot() +
  geom_jitter(
    data = among_site,
    aes(x = cluster, y = star_density),
    width = 0.12,
    height = 0,
    size = 2.5,
    alpha = 0.35) +
  geom_errorbar(
    data = dock_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    width = 0.12,
    linewidth = 0.9) +
  geom_point(data = dock_pred, aes(x = x, y = predicted), size = 4) +
  labs(x = "Cluster", y = expression("Sea star density (stars " * m ^
                                            {
                                              -2
                                            } * ")")) +
  pub_theme1()
print(full_debris_number_cluster_plot)

#mega stitch time for density----

density_plot_a <-
  full_debris_number_plot +
  full_debris_number_rock_plot +
  full_debris_number_cluster_plot + 
  plot_layout(
    design = "AA
    BC",
    heights = c(1.5, 1),
    axes = "collect_y",
    axis_titles = "collect_y"
  )  + plot_annotation(tag_levels = "A")
print(density_plot_a)
ggsave(
  "figures/density_fulldebris.png",
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)
density_plot_b <-
  full_mussel_number_plot +
  full_mussel_number_rock_plot +
  full_mussel_number_slope_plot +
  plot_layout(
    design = "AA
    BC",
    heights = c(1.5, 1),
    axes = "collect_y",
    axis_titles = "collect_y"
  )  + plot_annotation(tag_levels = "A")
print(density_plot_b)
ggsave(
  "figures/density_fullmussel.png",
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)


#pycno time

#incorporating more of habitat characteristics inside of the site model----

#is variation of slope explained by what type of site it is? ----
slopesites <- glmmTMB(slope_scaled ~ type + (1 |
                                               pair), data = among_site)
simulateResiduals(slopesites, plot = TRUE)
summary(slopesites)

#run the model to see if slope differs systematically due to site

slopesites_null <- glmmTMB(slope_scaled ~ 1 + (1 |
                                                 pair), data = among_site)
#compare the 2 models

anova(slopesites_null, slopesites)

#is variation of rock cover explained by what type of site it is?
rocksites <- glmmTMB(rock_scaled ~ type + (1 |
                                             pair), data = among_site)
simulateResiduals(rocksites, plot = TRUE)
summary(rocksites)

#null

rocksites_null <- glmmTMB(rock_scaled ~ 1 + (1 |
                                               pair), data = among_site)

anova(rocksites, rocksites_null)
