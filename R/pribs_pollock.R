#' Code in development for pollock abundance within 100km of the Pribilof 
#' Islands for research on fur seal foraging success. Code managed and updated 
#' by Sophia Wassermann, based on sdmTMB code adapted from bs_indices.Rmd 
#' (developed by Lewis Barnett)

library(sdmTMB)
library(dplyr)
library(ggplot2)
library(here)
library(devtools)
library(sf)

# install_github("afsc-gap-products/akgfmaps", build_vignettes = TRUE)
library(akgfmaps)

# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
theme_set(theme_sleek())

# Get pollock CPUE data -------------------------------------------------------
this_year <- as.numeric(format(Sys.Date(), "%Y"))

kts <- 250  # Number of knots for the index model mesh

# Whether abundance will be in "numbers" or "biomass"
data_type <- "numbers"

# Make a new directory for the model output
results_wd <- here("results", data_type, paste0(kts, "kts"))
dir.create(here(results_wd), recursive = TRUE, showWarnings = FALSE)

# Read in data
if(data_type == "numbers") {
  dat <- read.csv(here("data", "pollock_num.csv"))
}

if(data_type == "biomass") {
  dat <- read.csv(here("data", "VAST_ddc_all_2025.csv"))
  colnames(dat)[3:5] <- c("lat", "lon", "cpue")
}

# detect if any years have occurrences at every haul and fix params as needed
if(data_type == "numbers") {
  mins <- dat %>% group_by(year) %>% summarize(min = min(cpue))
  if (sum(mins$min, na.rm = TRUE) == 0) {
    control = sdmTMBcontrol()
  } else {
    no_zero_yr <- as.integer(mins %>% filter(min > 0) %>% select(year))
    # set up map and fix value of p(occurrence) to slightly less than 1:
    yrs <- sort(unique(factor(dat$year)))
    .map <- seq_along(yrs)
    .map[yrs %in% no_zero_yr] <- NA
    .map <- factor(.map)
    .start <- rep(0, length(yrs))
    .start[yrs %in% no_zero_yr] <- 20
    
    control =  sdmTMBcontrol(
      map = list(b_j = .map),
      start = list(b_j = .start)
    )
  }
}

# Set up cold pool covariate
# devtools::install_github("afsc-gap-products/coldpool")
env_df <- coldpool::cold_pool_index
if(max(!this_year %in% env_df$YEAR < this_year)) {
  print("You may need to remove and reinstall the coldpool package to get the newest data!")
}

env <- cbind(env_df, env = scale(coldpool::cold_pool_index$AREA_LTE2_KM2)) %>%
  mutate(year = as.integer(YEAR)) %>%
  select(year, env)

dat <- left_join(dat, env, by = "year") 

# Final data manipulation steps
dat$year_f <- as.factor(dat$year)
dat <- add_utm_columns(dat, ll_names = c("lon", "lat"), utm_crs = 32602, units = "km")

# Fit model (if needed) -------------------------------------------------------
f1 <- here(results_wd, "fit.RDS")
if (file.exists(f1)) {
  fit <- readRDS(f1)
  } else {
    mesh <-  make_mesh(dat, xy_cols = c("X", "Y"), 
                       mesh = fmesher::fm_as_fm(readRDS(file = here("shapefiles", 
                                                                    paste0("ebs_vast_mesh_", kts, "_knots.RDS")))))
    fit <- sdmTMB( 
      cpue ~ 0 + year_f,
      spatial_varying = ~ env,
      data = dat, 
      mesh = mesh,
      family = delta_gamma(type = "poisson-link"), 
      time = "year", 
      spatial = "on",
      spatiotemporal = "ar1",
      extra_time = 2020L, 
      silent = FALSE,
      anisotropy = TRUE,
      control = control
    )
  }

# Check fit
sanity(fit)
summary(fit)
saveRDS(fit, file = here(results_wd, "fit.RDS"))

# Make predictions and index --------------------------------------------------
# Read in fit if object is not already in environment
if(!exists("fit")) {
  fit <- readRDS(here(results_wd, "fit.RDS"))
}

# Read in polygons
island_grids <- read.csv(here("shapefiles", "processed", "island_complex_grids.csv"))

# Calculate index for each area
index_by_area <- function(strat = "all", reg) {
  df <- island_grids %>% filter(stratum == strat, region == reg)
  # Create a folder for each index area
  dir_name <- paste0("radius", reg, strat)
  dir.create(here(results_wd, dir_name), recursive = TRUE, showWarnings = FALSE)
  
  # replicate prediction grid for each year in data
  pred_grid <- replicate_df(data.frame(df), "year_f", unique(dat$year_f))
  pred_grid$year <- as.integer(as.character(factor(pred_grid$year_f)))
  
  # join with environmental covariate (cold pool)
  pred_grid <- left_join(pred_grid, env, by = "year")
  
  # get prediction
  p <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE)
  save(p, file = here(results_wd, dir_name, "pred.Rdata"))
  
  # get index
  ind <- get_index(p, bias_correct = TRUE, area = pred_grid$area_km2)
  ind$stratum <- strat
  ind$region <- reg
  write.csv(ind, file = here(results_wd, dir_name, "index.csv"), row.names = FALSE)
  print(paste0("Completed index for ", reg, " ", strat))
  
  return(ind)
}

indices <- bind_rows(index_by_area(reg = "STG"),
                     index_by_area(reg = "STG_North"),
                     index_by_area(reg = "STG_South"),
                     index_by_area(reg = "STP"),
                     index_by_area(reg = "STP_East"),
                     index_by_area(reg = "STP_EB"),
                     index_by_area(reg = "STP_Reef_Point"))

# Plot index, scaled from kg to Mt
ggplot(indices, aes(x = year, y = (est / 1e9))) +
  geom_line() +
  ylim(0, NA) +
  geom_ribbon(aes(ymin = (lwr / 1e9), ymax = (upr / 1e9)), alpha = 0.4) +
  xlab("") + 
  facet_wrap(~region, scales = "free") +
  {
    if (data_type == "numbers") {
      ylab("Abundance (billions)")
    } else if (data_type == "biomass") {
      ylab("Biomass (Mt)")
    }
  }
ggsave(file = here(results_wd, "index.png"), 
       height = 6, width = 10, units = "in")

# Plot predicted density maps and fit diagnostics -----------------------------
# q-q plot
pdf(file = here(results_wd, "qq.pdf"),
    width = 5, height = 5)
sims <- simulate(fit, nsim = 500, type = "mle-mvn") 
sims |> dharma_residuals(fit, test_uniformity = FALSE) 
# previous is not working: 
# Error in `predict()`:
# ! Prediction offset vector does not equal number of rows in prediction dataset.
dev.off()

#residuals on map plot, by year
resids <- sims |>
  dharma_residuals(fit, test_uniformity = FALSE, return_DHARMa = TRUE)
fit$data$resids <- resids$scaledResiduals

ggplot(subset(fit$data, !is.na(resids) & is.finite(resids)), aes(X, Y, col = resids)) +
  scale_colour_gradient2(name = "residuals", midpoint = 0.5) +
  geom_point(size = 0.7) +
  scale_x_continuous(breaks = c(250, 750)) +
  scale_y_continuous(breaks = c(6000, 6500, 7000)) +
  facet_wrap(~year) +
  coord_fixed() 
ggsave(file = here(results_wd, "residuals_map.pdf"),
       height = 9, width = 6.5, units = "in")

# predictions on map plot, by year
for(i in 1:nrow(final_combined_hr_polygons_projected_sf)) {
  dir_name <- paste0("radius", final_combined_hr_polygons_projected_sf$associated_circle_radius_meters[i])
  load(here(results_wd, dir_name, "pred.Rdata"))
  p <- p$data
  p$radius <- final_combined_hr_polygons_projected_sf$associated_circle_radius_meters[i]
  pred_map <- ggplot(p, aes(X, Y, fill = exp(est1 + est2))) +
    geom_tile(width = 10, height = 10) +
    scale_fill_viridis_c(trans = "sqrt", name = "") +
    scale_x_continuous(breaks = c(250, 750)) +
    scale_y_continuous(breaks = c(6000, 6500, 7000)) +
    facet_wrap(~year) +
    coord_fixed() +
    {
      if (data_type == "numbers") {
        ggtitle("Predicted densitites (numbers / square km)")
      } else if (data_type == "biomass") {
        ggtitle("Predicted densitites (kg / square km)")
      }
    }
  ggsave(pred_map, file = here(results_wd, dir_name, "predictions_map.pdf"),
         height = 7, width = 7, units = "in")
}
