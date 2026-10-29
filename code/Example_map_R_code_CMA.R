
## --------------- Example code for simple map of site locations
## --------------- Author: Claire Attridge

#adapted for docks and sea stars project

#

# Loading necessary packages

library(tidyverse)
library(MetBrewer)
library(sf)
library(ggspatial)

#

#### Initial work on the base map ----


### Step 1) Setting the map projections


# Establishing variables for my geometric and projected coordinate systems

proj <- st_crs(3005) # ESPG code for BC Albers projection
latlon <- st_crs(4326) # for baseline WGS84 geometry



### Step 2) Loading the base map shapefile and projecting into new coordinate system

## Loading the Hakai BC shapefile
# Setting the incoming crs as WGS84 (I.e., The imported map is in WGS84, so we need to establish this)
# Then, reprojecting the geometry into BC/Albers (NAD83)

## SET THE CORRECT DESTINATION FOR THIS CALL

land <- sf::st_read("data/shape_files/COAST_TEST2.shp") %>%
  st_as_sf(crs=latlon) %>%
  st_transform(proj)


## Setting bounds to a wide extent of Barkley Sound
## MANIPULATE THESE AS NEEDED TO FIT YOUR STUDY AREA WITHIN BC

ymax <- 49.47
ymin <- 49.25
xmax <- -122.82
xmin <- -123.00

## Making corners for the extent

corners <- st_multipoint(matrix(c(xmax,xmin,ymax,ymin),ncol=2)) %>% 
  st_sfc(crs=latlon) %>%
  st_sf() %>%
  st_transform(proj) 
plot(corners)


## Cropping BC map to the extent

land_cropped <- land %>%
  st_crop(st_bbox(corners))


## Saving the cropped map as a geopackage
## A geopackage polygon will essentially behave like a .shp file, but is a better way to store
## Save this to wherever you are currently working (I.e., Your R project or file folder)

write_sf(land_cropped, "data/Land_polygon_cropped.gpkg", overwrite=T)


#


#### Loading the site location/information dataframe ----

## Loading my sites and associated biodiversity info columns

data <- read_csv("data/stars_site - Sheet1 (3).csv")


## Making the dataframe into a shapefile
# Setting the incoming crs as WGS84 (because the original geometry was recorded in lat/lon)
# Then, reprojecting the geometry into BC/Albers (NAD83)

## In this example, I keep "sitetype" and "year" as column categories that I may be 
# interested to plot.
## I keep "site_code" and "site_name" in case I want to label my site points with text.

data_shp <- st_as_sf(data, coords = c("lon", "lat"),
                     crs=latlon) %>%
  st_transform(proj) %>%
  mutate(
         cluster = as.factor(cluster),
         site = as.factor(site),
         type = as.factor(type))

#

#### Plotting map of site locations ----

## Loading the cropped land map

## SET THE CORRECT DESTINATION FOR THIS CALL
## I.e., Call from wherever you are currently working (I.e., Your R project or file folder)

land_plot <- sf::st_read("data/Land_polygon_cropped.gpkg")


## Creating ggplot of my sites (n = 35)
## Saving as a .tiff file

## Lots of different elements going on in the theme() call
## Try adding/subtracting them or looking them up online to see what they do!

## When saving a file this way, the tiff() call 'opens' an image window
## then, whatever ggplot you run will be saved into that open call of tiff().
## The final dev.off() call is necessary to say that 'I'm done', which saves the
## plot to the location specified in the tiff() call

## SET THE CORRECT DESTINATION FOR THIS SAVE CALL

tiff(file="./figures/Map_Sitetype.tiff", height = 7, width = 6.5, units = "in", res=400)


annotation_scale(
  location = "br",
  style = "ticks",
  width_hint = 0.2,
  unit_category = "metric",
  plot_unit = "m"
)

land_plot <- sf::st_transform(land_plot, 3005)
data_shp  <- sf::st_transform(data_shp, 3005)


st_bbox(data_shp)
bb <- st_bbox(data_shp)

coord_sf(
  xlim = c(bb["xmin"] - 500,bb["xmax"] + 500),
  ylim = c(bb["ymin"] - 300, bb["ymax"] + 300),
  expand = FALSE
)
ggplot(land_plot) +
  geom_sf(fill = "grey65", color = NA) +
  geom_sf(
    data = data_shp,
    shape = 21,
    size = 3.0,
    aes(fill = type)
  ) +
  scale_fill_manual(values = met.brewer("Hokusai1", 3)) +
  
  coord_sf(
    xlim = c(bb["xmin"] - 700, bb["xmax"] + 700),
    ylim = c(bb["ymin"] - 300, bb["ymax"] + 300),
    expand = FALSE
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_line(colour = "transparent"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(color = "black", size = 10),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 10),
    plot.margin = unit(c(0,0.4,0.4,0.4), "cm"),
    panel.border = element_rect(colour = "black", fill = NA)
  ) +
  annotation_scale(
    location = "br",
    style = "ticks",
    width_hint = 0.2
  ) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    height = unit(0.4, "in"),
    width = unit(0.4, "in"),
    pad_x = unit(0.1, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  )

#here's the map
main_map <- ggplot(land_plot) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(
    data = data_shp,
    aes(fill = type, shape = type),
    size = 2.5,
    color = "black",
    stroke = 0.4) +
  scale_fill_manual(
    values = c(
      "dock" = "#F28E2B",
      "inter-dock" = "#59A14F",
      "no-dock" = "#4E79A7"),
    labels = c(
      "dock" = "Dock",
      "inter-dock" = "Inter-dock",
      "no-dock" = "No-dock")) +
  scale_shape_manual(
    values = c(
      "dock" = 21,
      "inter-dock" = 24,
      "no-dock" = 22),
    labels = c(
      "dock" = "Dock",
      "inter-dock" = "Inter-dock",
      "no-dock" = "No-dock")) +
  coord_sf(
    xlim = c(bb["xmin"] - 1000, bb["xmax"] + 1000),
    ylim = c(bb["ymin"] - 500, bb["ymax"] + 500),
    expand = FALSE) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "#A3D1D4", color = NA),
    panel.grid.major = element_line(colour = "transparent"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(color = "black", size = 8),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 8),
    plot.margin = unit(c(0,0.4,0.4,0.4), "cm"),
    panel.border = element_rect(colour = "black", fill = NA)) +
  annotation_scale(location = "br", style = "ticks", width_hint = 0.2) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    height = unit(0.28, "in"),
    width = unit(0.28, "in"),
    pad_x = unit(0.1, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering)

main_map

ggsave("figures/map.png",  width = 6.5, height = 7, units = "in", dpi = 300)


main_extent <- st_bbox(
  c(
    xmin = bb[["xmin"]] - 1000,
    xmax = bb[["xmax"]] + 1000,
    ymin = bb[["ymin"]] - 500,
    ymax = bb[["ymax"]] + 500),
  crs = st_crs(data_shp)) %>%
  st_as_sfc()

#this is the inset!!!

#---------------------------

inset_bb <- st_bbox(land_plot)

inset_map <- ggplot(land_plot) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(
    data = main_extent,
    fill = NA,
    color = "red",
    linewidth = 0.8) +
  coord_sf(
    xlim = c(inset_bb["xmin"], inset_bb["xmax"]),
    ylim = c(inset_bb["ymin"], inset_bb["ymax"]),
    expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(
      fill = "#A6D4D8",
      colour = NA),
    plot.background = element_rect(
      fill = "#A6D4D8",
      colour = NA),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8),
    plot.margin = margin(0, 0, 0, 0)) +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8),
    plot.margin = margin(0, 0, 0, 0))

inset_map

final_map <- main_map +
  inset_element(
    inset_map,
    left = 0.64,
    bottom = 0.27,
    right = 0.99,
    top = 0.66)

final_map

ggsave("figures/map.png",  width = 6.5, height = 7, units = "in", dpi = 300)
