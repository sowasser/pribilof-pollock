#' Code in development for pollock abundance within 100km of the Pribilof 
#' Islands for research on fur seal foraging success. Code managed and updated 
#' by Sophia Wassermann, based on sdmTMB code adapted from bs_indices.Rmd 
#' (developed by Lewis Barnett)

library(sdmTMB)
library(dplyr)
library(ggplot2)
library(here)
library(sf)

# pak::pkg_install("afsc-gap-products/coldpool")

# pak::pkg_install("afsc-gap-products/akgfmaps", build_vignettes = TRUE)
library(akgfmaps)

# pak::pkg_install("seananderson/ggsidekick")
library(ggsidekick)
theme_set(theme_sleek())

# Get pollock CPUE data -------------------------------------------------------
this_year <- as.numeric(format(Sys.Date(), "%Y"))
# this_year <- 2025

kts <- 250  # Number of knots for the index model mesh

# Whether abundance will be in "numbers" or "biomass"
data_type <- "numbers"

# Make a new directory for the model output
results_wd <- here("results", data_type, paste0(kts, "kts"))
dir.create(here(results_wd), recursive = TRUE, showWarnings = FALSE)

# Read in data
if(data_type == "numbers") {
  dat <- read.csv(here("data", "ddc_num.csv"))
}

if(data_type == "biomass") {
  file <- here("data", paste0("VAST_ddc_all_", this_year, ".csv"))
  if(!file.exists(file)) {
    dir.create(here("data"), recursive = TRUE, showWarnings = FALSE)
    library(googledrive)
    
    # Update stale Google Drive authorization if needed
    gdrive_email <- rstudioapi::showPrompt(title = "Email",
                                           message = "Email for Google Drive",
                                           default = "")
    drive_auth(token = gargle::credentials_user_oauth2(
      scopes = "https://www.googleapis.com/auth/drive",
      email = gdrive_email))
    drive_user()  # check user account
    
    # Download from google drive
    file_id <- "1Kn_wftvGyylCTMTD-daHpFqoxanTMBXz"
    
    drive_download(file = as_id(file_id),
                   path = here("data", paste0("VAST_ddc_all_", this_year, ".csv")),
                   overwrite = TRUE)
    dat <- read.csv(file)
  } else {
    dat <- read.csv(file)
  }
  colnames(dat)[3:5] <- c("lat", "lon", "cpue")
}

# Check for NAs
unique(is.na(dat))

# detect if any years have occurrences at every haul and fix params as needed
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

# Set up cold pool covariate
# pak::pak("afsc-gap-products/coldpool")
env_df <- coldpool::cold_pool_index
# if(max(!this_year %in% env_df$YEAR < this_year)) {
#   print("You may need to remove and reinstall the coldpool package to get the newest data!")
# }

env <- cbind(env_df, env = scale(coldpool::cold_pool_index$AREA_LTE2_KM2)) %>%
  mutate(year = as.integer(YEAR)) %>%
  select(year, env)

dat <- left_join(dat, env, by = "year") 

# Final data manipulation steps
dat$year_f <- as.factor(dat$year)
dat <- add_utm_columns(dat, ll_names = c("lon", "lat"), utm_crs = 32602, units = "km")

# Fit model (if needed) -------------------------------------------------------
start <- Sys.time()
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
end <- Sys.time()
fit_time <- difftime(end, start, units = "hours")
fit_time

# Check fit
sanity(fit)
summary(fit)

# Safer file saving (hopefully)
temp_f1 <- tempfile(fileext = ".RDS")
saveRDS(fit, file = here(results_wd, temp_f1), compress = TRUE)

# Extra optimization if needed
# fit_opt <- run_extra_optimization(fit)
# sanity(fit_opt)

# Make predictions and index --------------------------------------------------
# Read in fit if object is not already in environment
if(!exists("fit")) {
  fit <- sdmTMB:::reload_model(here(results_wd, "fit.RDS"))
}

# Calculate index for each area
index_by_area <- function(reg) {
  start <- Sys.time()  # start timer
  # Load in prediction grid created in prib_areas.R
  load(here("shapefiles", "processed", reg, "grid_5nm.Rdata"))  # this is grid_list

  # Subset grid_list to include the "all" polygon and 50-125km.
  grid_list <- grid_list[c(1:6)]
  
  # Create a folder for each index area
  if(!dir.exists(here(results_wd, reg))) {
    dir.create(here(results_wd, reg), recursive = TRUE, showWarnings = FALSE)
  }
  
  ind_list <- vector("list", length = length(grid_list))
  for(i in 1:length(grid_list)) {
    df <- grid_list[[i]]
    stratum <- unique(df$stratum)  # for later labelling

    # Create prediction grid
    pred_grid <- replicate_df(data.frame(df), "year_f", unique(dat$year_f))
    pred_grid$year <- as.integer(as.character(factor(pred_grid$year_f)))
      
    # join with environmental covariate (cold pool)
    pred_grid <- left_join(pred_grid, env, by = "year")
    
    # Get predictions (if hasn't been done already) or load saved file
    pred_file <- here(results_wd, reg, paste0("pred_", stratum,".Rdata"))
    if (!file.exists(pred_file)) {
      # get predictions
      print(paste0("Predicting for ", reg, " ", stratum))
      p <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE,
                  offset = rep(0, nrow(pred_grid)))
      save(p, file = here(results_wd, reg, paste0("pred_", stratum,".Rdata")))

    } else {
      print(paste0("Predictions already exist for ", reg, " ", stratum, "; loading from file"))
      load(pred_file)  # this is p
    }

    # get index
    ind_file <- here(results_wd, reg, paste0("index_", stratum, ".csv"))
    ind <- get_index(p, bias_correct = TRUE, area = p$data$area_km2)
    ind$stratum <- stratum
    ind$region <- reg

    write.csv(ind, 
        file = here(results_wd, reg, paste0("index_", stratum, ".csv")), 
        row.names = FALSE
      )
    
        # Map of predicted density
    pdata <- p$data
    # pdata$stratum <- stratum
    pred_map <- ggplot(pdata, aes(X, Y, fill = est1 + est2)) +
      geom_tile(width = 10, height = 10) +
      scale_fill_viridis_c(name = "") +
      scale_x_continuous(breaks = c(250, 750)) +
      scale_y_continuous(breaks = c(6000, 6500, 7000)) +
      facet_wrap(~year) +
      xlab("") + ylab("") +
      coord_fixed() +
      {
        if (data_type == "numbers") {
          ggtitle(expression("Predicted log density (numbers / km"^2*")"))
        } else if (data_type == "biomass") {
          ggtitle(expression("Predicted log density (kg / km"^2*")"))
        }
      }
    
    ggsave(pred_map, file = here(results_wd, reg, paste0("log_pred_map_", stratum, ".pdf")),
           height = 7, width = 7, units = "in")
    
      ind_list[[i]] <- ind
  }

  ind_df <- bind_rows(ind_list)

  end <- Sys.time()
  print(paste0("Completed index for ", reg, " in ", round(difftime(end, start, units = "hours"), 2), " hours"))
  
  return(ind_df)
}

stg <- index_by_area("STG")
stg_n <- index_by_area("STG_North")
stg_s <- index_by_area("STG_South")
stp <- index_by_area("STP")
stp_e <- index_by_area("STP_East")
stp_eb <- index_by_area("STP_EB")
stp_rp <- index_by_area("STP_Reef_Point")

# Alternatively, read in existing index files ---------------------------------
# read_index <- function(reg) {
#   files <- list.files(
#     path = here(results_wd, reg),
#     pattern = "^index_.*\\.csv$", 
#     full.names = TRUE
#   )
  
#   # Check if any files were found
#   if (length(files) == 0) {
#     warning("No index CSV files found in: ", here(results_wd, reg))
#     return(NULL)
#   }
  
#   message("Reading ", length(files), " index CSVs for region: ", here(results_wd, reg))
  
#   # Read and combine all CSVs into a single data frame
#   combined_df <- files %>%
#   lapply(function(file) {
#     df <- read.csv(file)
#     # Force stratum to character type so bind_rows can combine numeric & text labels
#     df$stratum <- as.character(df$stratum)
#     return(df)
#   }) %>%
#   bind_rows()
  
#   return(combined_df)
# }

# stg <- read_index("STG")
# stg_n <- read_index("STG_North")
# stg_s <- read_index("STG_South")
# stp <- read_index("STP")
# stp_e <- read_index("STP_East")
# stp_eb <- read_index("STP_EB")
# stp_rp <- read_index("STP_Reef_Point")

# Combine all indices and plot ------------------------------------------------
indices <- bind_rows(stg, stg_n, stg_s, stp, stp_e, stp_eb, stp_rp)
indices$stratum <- factor(indices$stratum, 
                          levels = c("all", "250", "225", "200", "175", "150", 
                                     "125", "100", "75", "50", "25"))

write.csv(indices, file = here(results_wd, paste0("index_", data_type, ".csv")), row.names = FALSE)

ggplot(indices, aes(x = year, y = (est / 1e9))) +
  geom_line() +
  geom_ribbon(aes(ymin = (lwr / 1e9), ymax = (upr / 1e9)), alpha = 0.4) +
  ylim(0, NA) +
  xlab("") + 
  facet_grid(stratum ~ region, scales = "free") +
  {
    if (data_type == "numbers") {
      ylab("Abundance (billions)")
    } else if (data_type == "biomass") {
      ylab("Biomass (Mt)")
    }
  }
ggsave(file = here(results_wd, "index.png"), 
       height = 10, width = 12, units = "in")

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
  xlab("") + ylab("") +
  coord_fixed() 
ggsave(file = here(results_wd, "residuals_map.pdf"),
       height = 9, width = 6.5, units = "in")


# # Extra code for just plotting predicted density maps 
# den_maps <- function(reg) {
#   load(here("shapefiles", "processed", reg, "grid.Rdata"))  # this is grid_list
  
#   for(i in 1:length(grid_list)) {
#     df <- grid_list[[i]]
#     stratum <- unique(df$stratum)  # for later labelling

#     load(file = here(results_wd, reg, paste0("pred_", stratum, ".Rdata")))  # this is p

#     # Map of predicted density
#     pdata <- p$data
#     pred_map <- ggplot(pdata, aes(X, Y, fill = est1 + est2)) +
#       geom_tile(width = 10, height = 10) +
#       scale_fill_viridis_c(name = "") +
#       scale_x_continuous(breaks = c(250, 750)) +
#       scale_y_continuous(breaks = c(6000, 6500, 7000)) +
#       facet_wrap(~ year) +
#       xlab("") + ylab("") +
#       coord_fixed() +
#       {
#         if (data_type == "numbers") {
#           ggtitle("Predicted log density (numbers / square km)")
#         } else if (data_type == "biomass") {
#           ggtitle("Predicted log density (kg / square km)")
#         }
#       }
    
#     ggsave(pred_map, file = here(results_wd, reg, paste0("pred_map_", stratum, ".pdf")),
#            height = 7, width = 7, units = "in")
#   }
# }

# for(reg in c("STG", "STG_North", "STG_South", "STP", "STP_East", "STP_EB", "STP_Reef_Point")) {
#   den_maps(reg)
# }
