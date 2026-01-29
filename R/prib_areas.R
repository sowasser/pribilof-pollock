#' Polygons for the Pribilof Islands have been separated into regions, then into
#' radii out from a central point. There are 11 polygons within each region and
#' 7 regions. This script reads in the regions and creates a grid for each 
#' radius (sub-area) within each region. The code produces a plot of the sub-
#' area and combines all grids into a single dataframe for later use.

library(here)
library(dplyr)
library(sf)

# install_github("afsc-gap-products/akgfmaps", build_vignettes = TRUE)
library(akgfmaps)

# ggplot theme
# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
theme_set(theme_sleek())

# Run grid creation at the sub-area level for each polygon
region_grids <- function(polygon, label, save_plot = TRUE) {
  df <- polygon
  
  # Clean up labels
  df$associated_circle_radius_meters <- df$associated_circle_radius_meters / 1000
  df$associated_circle_radius_meters[1] <- "all"
  
  # Create expansion grid for each sub-area
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
  
  # Better stratum layers
  grids$stratum  <- factor(grids$stratum,
                           levels = c("25", "50", "75", "100", "125", "150", "175", "200", "225", "250", "all"),
                           labels = c("25km", "50km", "75km", "100km", "125km", "150km", "175km", "200km", "225km", "250km", "all")
  )
  grids$region <- label
  
  plot <- ggplot(grids, aes(X, Y, colour = area_km2)) +
    geom_tile(width = 2, height = 2, fill = NA) +
    scale_colour_viridis_c(direction = -1) +
    geom_point(size = 0.5) +
    coord_fixed() +
    xlab("") + ylab("") +
    facet_wrap(~stratum) +
    ggtitle(label)
  
  if(save_plot == TRUE) {
    ggsave(plot, file = here("shapefiles", "processed", paste0("pred_grids_", label, ".png")),
           height = 6, width = 7.5, units = c("in"))
  }
  
  return(grids)
}

# Load in (and get vector of names) for each region
regions <- load(here("shapefiles", "Island_complex_polygons.RData"))
regions <- sub("_polygon$", "", regions)  # clean up names

# Run function and combine all grids into a data frame
all_grids <- bind_rows(
  region_grids(STG_North_polygon, regions[1]),
  region_grids(STG_polygon, regions[2]),
  region_grids(STG_South_polygon, regions[3]),
  region_grids(STP_East_polygon, regions[4]),
  region_grids(STP_EB_polygon, regions[5]),
  region_grids(STP_polygon, regions[6]),
  region_grids(STP_Reef_Point_polygon, regions[7])
)

# Save as .csv for use in pribs_pollock.R script
write.csv(all_grids, 
          file = here("shapefiles", "processed", "island_complex_grids.csv"),
          row.names = FALSE)
