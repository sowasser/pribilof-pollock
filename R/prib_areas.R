# Looking through and setting up the shapefiles for Pribilof areas

library(here)
library(dplyr)
library(sf)
library(akgfmaps)

load(here("shapefiles", "Island_complex_polygons.RData"))

df <- STG_polygon

# Clean up labels
df$associated_circle_radius_meters <- df$associated_circle_radius_meters / 1000
df$associated_circle_radius_meters[1] <- "all"


# Create expansion grid for each area
grid_by_area <- function(area) {
  polygon <- df$geometry[area]
  grid <- make_2d_grid(obj = polygon,
                       resolution = c(3704, 3704),  # default resolution - 2x2nm
                       output_type = "point",
                       include_tile_center = TRUE) %>%
    st_transform(crs = "EPSG:32602") 
  
  grid[, c("LON_UTM", "LAT_UTM")] <- st_coordinates(grid)
  
  grid <- data.frame(grid) %>%
    select(LON_UTM, LAT_UTM, AREA) %>%
    mutate(X = LON_UTM / 1000,
           Y = LAT_UTM / 1000,
           area_km2 = as.numeric(AREA)/1e6) %>%
    select(X, Y, area_km2)
  grid <- as.data.frame(as.matrix(grid)) # drop attributes
  grid$stratum <- df$associated_circle_radius_meters[area]
  return(grid)
}
grid_list <- lapply(1:nrow(df), grid_by_area)
grids <- do.call(rbind, grid_list)
grids$stratum  <- factor(grids$stratum,
                         levels = c("25", "50", "75", "100", "125", "150", "175", "200", "225", "250", "all"),
                         labels = c("25km", "50km", "75km", "100km", "125km", "150km", "175km", "200km", "225km", "250km", "all")
)

ggplot(grids, aes(X, Y, colour = area_km2)) +
  geom_tile(width = 2, height = 2, fill = NA) +
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.5) +
  coord_fixed() +
  xlab("") + ylab("") +
  facet_wrap(~stratum)

ggsave(file = here(results_wd, "pred_grids.png"),
       height = 6, width = 7.5, units = c("in"))
