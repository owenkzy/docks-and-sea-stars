#Code for If you build it, will they come? Coastal residential infrastructure affects sea star communities

#load in packages----
library(tidyverse)
library(ggforce)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(performance)
library(rstatix)
library(emmeans)
library(effectsize)
library(MuMIn)
library(vegan)
library(ggalt)
library(ggrepel)
library(patchwork)
library(Hmisc)

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


#among sites df----

among <- site %>%
  full_join(stars_rove,
            by = c("site", "date"),
            relationship = "many-to-many")  %>%
  full_join(belt_site,
            by = c("site", "date"),
            relationship = "many-to-many") %>%
  full_join(debris,
            by = c("site", "date"),
            relationship = "many-to-many") %>%
  left_join(photo, by = "site")

#get it so only 1 row per observation
among_site <- among %>% distinct(site, .keep_all = TRUE)  %>%
  mutate(
    density = total_abundance / (transect_length * 2),
    slope_scaled = as.numeric(scale(slope)),
    distance_scaled = as.numeric(scale(dist_from_dock_new)),
    rock_scaled = as.numeric(scale(rock_percent)),
    lat_scaled = as.numeric(scale(lat)),
    area_transect = transect_length * 2,
    mussel_scaled = as.numeric(scale(mussel_percent)),
    debris_score = as.factor(debris_score),
             mussel_score = as.factor(mussel_score)) %>%
  distinct()

among_site$type <- relevel(factor(among_site$type), ref = "no-dock")
#general stats below----

#dock distance info
among_site %>%
  group_by(type) %>%
  summarise(
    n = n(),
    mean_distance = mean(dist_from_dock_new, na.rm = TRUE),
    sd_distance = sd(dist_from_dock_new, na.rm = TRUE),
    min_distance = min(dist_from_dock_new, na.rm = TRUE),
    max_distance = max(dist_from_dock_new, na.rm = TRUE))

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


#fishers exact test, have to do it pairwise



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
null_density <- glmmTMB(
  total_abundance ~ 1 + (1 |
                           pair) + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(null_density)
simulateResiduals(null_density, plot = TRUE)


#spatial model including offset

spatial_density <- glmmTMB(
  total_abundance ~ cluster + offset(log(area_transect)) + (1 | pair),
  family = nbinom2(link = "log"),
  data = among_site
)
summary(spatial_density)
simulateResiduals(spatial_density, plot = TRUE)
check_collinearity(spatial_density)

#latitude model

lat_density <- glmmTMB(
  total_abundance ~ lat_scaled + (1 | pair) +
    offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)

#habitat model

habitat_density <- glmmTMB(
  total_abundance ~ slope_scaled + rock_scaled + cluster + (1 | pair) +
    offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site
)


#dock model

dock_debris_density <- glmmTMB(
  total_abundance ~ debris_score + 
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site)



dock_mussel_density <- glmmTMB(
  total_abundance ~ mussel_score + 
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site)

dock_density <- glmmTMB(
  total_abundance ~ type + 
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site)

dock_distance_density <- glmmTMB(
  total_abundance ~ distance_scaled + 
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site)


#full model
full_density <- glmmTMB(
  total_abundance ~ type + slope_scaled + rock_scaled  +
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site)

full_debris_density <- glmmTMB(
  total_abundance ~ debris_score + slope_scaled + rock_scaled  +
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site)

summary(full_debris_density)
simulateResiduals(full_debris_density, plot = TRUE)

full_mussel_density <- glmmTMB(
  total_abundance ~ mussel_score + slope_scaled + rock_scaled  +
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site)

summary(full_mussel_density)

full_distance_density <- glmmTMB(
  total_abundance ~ distance_scaled + slope_scaled + rock_scaled  +
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site)


#aic time

model_set <- list(
  null = null_density,
  spatial = spatial_density,
  dock = dock_density,
  habitat = habitat_density,
  full = full_density,
  lat = lat_density,
  full_mussel = full_mussel_density,
  full_debris = full_debris_density,
  dock_mussel = dock_mussel_density,
  dock_debris = dock_debris_density,
  full_distance = full_distance_density,
  dock_distance = dock_distance_density)

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
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site)

summary(full_debris_number)
full_mussel_number <- glmmTMB(
  total_abundance ~ mussel_score + slope_scaled + rock_scaled  +
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = among_site)

summary(full_mussel_number)

#AIC COMPETING SIZE TIME-----

#this is to get the weighted by density size estimate
ochre <- read_csv("data/stars_belt - Sheet1 (1).csv")  %>%
  group_by(site, date, species) %>%
  summarise(average_size = mean(size, na.rm = TRUE),
            species_abundance = n(), .groups = "drop") %>% 
  group_by(site,date) %>%
  mutate(total_abundance = sum(species_abundance, na.rm = TRUE),
         proportion = species_abundance/total_abundance,
         weighted_size = average_size*proportion) %>% 
  ungroup()

new_size <- ochre %>%
  group_by(site, date) %>%
  summarise(site_average = sum(weighted_size, na.rm = TRUE),
            .groups = "drop")


among_site <- among_site %>% 
  full_join(new_size, by = c("site","date"))

#null model
null_size <- glmmTMB(site_average ~ 1 + (1 |pair),
                family = Gamma(link = "log"),
                data = among_site)
#spatial model

spatial_size <- glmmTMB(site_average ~ cluster + (1 | pair),
                   family = Gamma(link = "log"),
                   data = among_site)

#habitat model

habitat_size <- glmmTMB(
  site_average ~ slope_scaled + rock_scaled + (1 | pair),
  family = Gamma(link = "log"),
  data = among_site)

#latitude model

lat_size <- glmmTMB(site_average ~ lat_scaled + (1 | pair),
               family = Gamma(link = "log"),
               data = among_site)

#dock model and associated features models

dock_size <- glmmTMB(
  site_average ~ type + cluster + (1 |pair),
  family = Gamma(link = "log"),
  data = among_site)

dock_mussel_size <- glmmTMB(
  site_average ~ mussel_score + cluster + (1 |pair),
  family = Gamma(link = "log"),
  data = among_site)

dock_debris_size <- glmmTMB(
  site_average ~ debris_score + cluster + (1 |pair),
  family = Gamma(link = "log"),
  data = among_site)

dock_distance_size <- glmmTMB(
  site_average ~ distance_scaled  +
    (1 |pair) + cluster,
  family = Gamma(link = "log"),
  data = among_site)

#full model

full_size <- glmmTMB(
  site_average ~ cluster + type + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site)

full_mussel_size <- glmmTMB(
  site_average ~ cluster + mussel_score + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site)

full_debris_size <- glmmTMB(
  site_average ~ cluster + debris_score + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site)

full_distance_size <- glmmTMB(
  site_average ~ distance_scaled + slope_scaled + rock_scaled  +
    (1 |pair) + cluster,
  family = Gamma(link = "log"),
  data = among_site)

#aic time

model_set <- list(
  null = null_size,
  spatial = spatial_size,
  dock = dock_size,
  habitat = habitat_size,
  full = full_size,
  lat = lat_size,
  full_mussel = full_mussel_size,
  full_debris = full_debris_size,
  dock_mussel = dock_mussel_size,
  dock_debris = dock_debris_size,
  full_distance = full_distance_size,
  dock_distance = dock_distance_size)

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
  site_average ~ type + cluster + (1 |pair),
  family = Gamma(link = "log"),
  data = among_site)

summary(dock_size)

full_dock_size <- glmmTMB(
  site_average ~ cluster + type + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site)

summary(full_dock_size)
simulateResiduals(full_dock_size, plot = TRUE)

#adding more size stuff, so we know that the average size is larger, but why? Is it due to species or structures?----

belt_site1 <- read_csv("data/stars_belt - Sheet1 (1).csv") #%>%
#  group_by(site, date) %>%
#  mutate(average_size = mean(size, na.rm = TRUE)) %>%
#  group_by(site, date, species, average_size) %>%
#  summarise(species_abundance = n(), .groups = "drop") %>%
#  group_by(site, date) %>%
#  mutate(total_abundance = sum(species_abundance, na.rm = TRUE)) %>%
#  ungroup()

#adding pair back to it

sizing_df <- belt_site1 %>%
  full_join(site,
            by = c("site", "date"),
            relationship = "many-to-one")  %>%
  drop_na(size, type, species, pair, site) %>%
  filter(species %in% c("leather", "ochre", "mottled"))

#good enough for the model

sizing <- glmmTMB(size ~ type + species + (1|pair) + (1|site), dispformula = ~ species, data = sizing_df)
simulateResiduals(sizing, plot = TRUE)
summary(sizing)

#try it again with an interaction

sizing_int <- glmmTMB(size ~ type*species + (1|pair) + (1|site), dispformula = ~ species, data = sizing_df)
simulateResiduals(sizing_int, plot = TRUE)
summary(sizing_int)

anova(sizing, sizing_int)

#emmeans now

emm_size <- emmeans(sizing_int, ~ type | species)
pairs(emm_size, adjust = "tukey")

#AIC competing time RICHNESS----


#null model
null_rich <- glmmTMB(richness ~ 1 + (1 |pair),
                family = Gamma(link = "log"),
                data = among_site)
#spatial model

spatial_rich <- glmmTMB(richness ~ cluster + (1 | pair),
                   family = Gamma(link = "log"),
                   data = among_site)

#habitat model

habitat_rich <- glmmTMB(
  richness ~ slope_scaled + rock_scaled + (1 | pair),
  family = Gamma(link = "log"),
  data = among_site)

#latitude model

lat_rich <- glmmTMB(richness ~ lat_scaled + (1 | pair),
               family = Gamma(link = "log"),
               data = among_site)

#dock model

dock_rich <- glmmTMB(richness ~ type + cluster + (1 | pair),
                family = Gamma(link = "log"),
                data = among_site)
dock_mussel_rich <- glmmTMB(
  richness ~ mussel_score + cluster + (1 | pair),
  family = Gamma(link = "log"),
  data = among_site)

dock_debris_rich <- glmmTMB(
  richness ~ debris_score + cluster + (1 |pair),
  family = Gamma(link = "log"),
  data = among_site)

dock_distance_rich <- glmmTMB(
  richness ~ distance_scaled  +
    (1 |pair) + cluster,
  family = Gamma(link = "log"),
  data = among_site)

#full model

full_rich <- glmmTMB(
  richness ~ cluster + type + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site)

full_mussel_rich <- glmmTMB(
  richness ~ cluster + mussel_score + slope_scaled + rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site)

full_debris_rich <- glmmTMB(
  richness ~ cluster + debris_score + slope_scaled +  rock_scaled +
    (1 | pair),
  family = Gamma(link = "log"),
  data = among_site)

full_distance_rich <- glmmTMB(
  richness ~ distance_scaled + slope_scaled + rock_scaled  +
    (1 |
       pair) + cluster,
  family = Gamma(link = "log"),
  data = among_site)

#aic time


model_set <- list(
  null = null_rich,
  spatial = spatial_rich,
  dock = dock_rich,
  habitat = habitat_rich,
  full = full_rich,
  lat = lat_rich,
  full_mussel = full_mussel_rich,
  full_debris = full_debris_rich,
  dock_mussel = dock_mussel_rich,
  dock_debris = dock_debris_rich,
  full_distance = full_distance_rich,
  dock_distance = dock_distance_rich
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
community_dist <- vegdist(community, method = "bray")

community_dist

#nmds time

set.seed(123)

nmds <- metaMDS(
  community,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE)
#stress is 0,106 which is good! Can infer stuff
#let's make a stressplot

stressplot(
  nmds,
  pch = 19,
  p.col = "blue",
  l.col = "red",
  lwd = 2)

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


#plot below----

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
  permutations = 9999)

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
      yend = NMDS2 * arrow_scale),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.12, "cm")),
    linewidth = 0.4) +
  geom_text(
    data = vec_sig,
    aes(x = label_x, y = label_y, label = species),
    inherit.aes = FALSE,
    size = 4) +
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
    fill = guide_legend(nrow = 1))

nmds_final
ggsave(
  "figures/NMDS_final.png",
  nmds_final,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300)

#let's do some pycno stuff----

sun <- read_csv("data/pycno_data - Sheet1.csv") %>%
  group_by(site) %>%
  mutate(
    average_pycno_size = mean(pycno_size, na.rm = TRUE),
    average_pycno_size = replace(average_pycno_size, is.nan(average_pycno_size), NA_real_)
  ) %>%
  ungroup() %>%
  distinct(site, pycno_count, average_pycno_size) %>%
#leftjoin it to full data
  left_join(among_site, by = "site") %>%
  mutate(pycno_density = pycno_count / (transect_length * 6),
         transect_rove_area = transect_length * 2)

#model it AICC Pycno time using offset for count

#null model
null_pycno <- glmmTMB(
  pycno_count ~ 1 + 
    (1 |pair) + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun)

#spatial model including offset

spatial_pycno <- glmmTMB(
  pycno_count ~ cluster + offset(log(area_transect)) + (1 | pair),
  family = nbinom2(link = "log"),
  data = sun
)
#latitude model

lat_pycno <- glmmTMB(
  pycno_count ~ lat_scaled + (1 | pair) +
    offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun)


#habitat model

habitat_pycno <- glmmTMB(
  pycno_count ~ slope_scaled + rock_scaled + cluster + (1 | pair) +
    offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun)

#dock model

dock_debris_pycno <- glmmTMB(
  pycno_count ~ debris_score + 
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun)

dock_mussel_pycno <- glmmTMB(
  pycno_count ~ mussel_score + 
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun)

dock_pycno <- glmmTMB(
  pycno_count ~ type + 
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun)

dock_distance_pycno <- glmmTMB(
  pycno_count ~ distance_scaled + 
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun)

#full model

full_pycno <- glmmTMB(
  pycno_count ~ type + slope_scaled + rock_scaled  +
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun
)
full_debris_pycno <- glmmTMB(
  pycno_count ~ debris_score + slope_scaled + rock_scaled  +
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun)


full_mussel_pycno <- glmmTMB(
  pycno_count ~ mussel_score + slope_scaled + rock_scaled  +
    (1 |pair) + cluster + offset(log(area_transect)),
  family = nbinom2(link = "log"),
  data = sun)
summary(full_mussel)

full_distance_pycno <- glmmTMB(
  pycno_count ~ distance_scaled + slope_scaled + rock_scaled  +
    (1 | pair) + cluster + offset(log(area_transect)),
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


#dock and full_dock_size is the best model for size----
summary(dock_size)
simulateResiduals(dock_size, plot = TRUE)
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
  scale_x_discrete(labels = tools::toTitleCase) +
  labs(x = "Site type", y = expression("Abudance-weighted mean size (cm)")) +
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
simulateResiduals(full_dock_size, plot = TRUE)
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
  scale_x_discrete(labels = tools::toTitleCase) +
  labs(x = "Site type", y = expression("Abudance-weighted mean size (cm)")) +
  pub_theme1()
print(full_size_plot)

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
  scale_x_discrete(labels = toupper) +
  labs(x = "Cluster", y = expression("Abudance-weighted mean size (cm)")) +
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
    x = "Hard substrate cover (%)",
    y = expression("Abudance-weighted mean size (cm)")) +
  pub_theme1()
print(full_dock_size_rock_plot)

#stitch size again----

size_plot_a <-
  full_size_plot /
  full_dock_size_rock_plot +
  plot_layout(
    heights = c(1.5, 1),
    axes = "collect_y",
    axis_titles = "collect_y") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.location = "panel",
    plot.tag.position = c(0.03, 0.97),
    plot.tag = element_text(
      face = "bold",
      size = 14,
      hjust = 0,
      vjust = 1))
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
    axis_titles = "collect_y")  +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.location = "panel",
    plot.tag.position = c(0.03, 0.97),
    plot.tag = element_text(
      face = "bold",
      size = 14,
      hjust = 0,
      vjust = 1))
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
simulateResiduals(full_debris_number, plot = TRUE)
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
  labs(x = "Debris score", y = expression("Sea star density (individuals " * m ^
                                            {
                                              -2
                                            } * ")")) +
  pub_theme1()
print(full_debris_number_plot)

#graph for mussels
simulateResiduals(full_mussel_number, plot = TRUE)

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
  labs(x = "Mussel score", y = expression("Sea star density (individuals " * m ^
                                            {
                                              -2
                                            } * ")")) +
  pub_theme1()
print(full_mussel_number_plot)


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
    x = "Slope (°)",
    y = expression("Sea star density (individuals " * m^{-2} * ")")) +
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
    x = "Hard substrate cover (%)",
    y = expression("Sea star density (individuals " * m^{-2} * ")")) +
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
    x = "Hard substrate cover (%)",
    y = expression("Sea star density (individuals " * m^{-2} * ")")) +
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
  scale_x_discrete(labels = toupper) +
  labs(x = "Cluster", y = expression("Sea star density (individuals " * m ^
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
  )  + plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.location = "panel",
    plot.tag.position = c(0.03, 0.97),
    plot.tag = element_text(
      face = "bold",
      size = 14,
      hjust = 0,
      vjust = 1))
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
  )  + plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.location = "panel",
    plot.tag.position = c(0.03, 0.97),
    plot.tag = element_text(
      face = "bold",
      size = 14,
      hjust = 0,
      vjust = 1))
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
