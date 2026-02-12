#' Script for pulling CPUE by numbers for Pribilof pollock indices project. 
#' Based on the Prepare_Bering_data.R script.

# devtools::install_github(repo = "afsc-gap-products/gapindex@v3.0.2", 
#                          dependencies = TRUE)
library(gapindex)
library(here)
library(dplyr)
library(ggplot2)
library(ggsidekick)
theme_set(theme_sleek())

# Connect to Oracle
if (file.exists("Z:/Projects/ConnectToOracle.R")) {
  source("Z:/Projects/ConnectToOracle.R")
} else {
  # For those without a ConnectToOracle file
  channel <- gapindex::get_connected(check_access = F)
}

## checks to see if connection has been established
odbcGetInfo(channel)

# Pull catch and effort data --------------------------------------------------
species_code <- c(21740, 21741)

# First, pull data from the standard EBS stations
ebs_standard_data <- get_data(year_set = 1982:as.integer(format(Sys.Date(), "%Y")),
                              survey_set = "EBS",
                              spp_codes = species_code,
                              pull_lengths = FALSE, 
                              haul_type = 3, 
                              abundance_haul = "Y",
                              channel = channel,
                              remove_na_strata = TRUE)

#' Next, pull data from hauls that are not included in the design-based index
#' production (abundance_haul == "N") but are included in VAST. By default, the 
#' gapindex::get_data() function will filter out hauls with negative performance 
#' codes (i.e., poor-performing hauls).
ebs_other_data <- get_data(year_set = c(1994, 2001, 2005, 2006),
                           survey_set = "EBS",
                           spp_codes = species_code,
                           pull_lengths = FALSE, 
                           haul_type = 3, 
                           abundance_haul = "N",
                           channel = channel, 
                           remove_na_strata = TRUE)

# Combine the EBS standard and EBS other data into one list. 
ebs_data <- list(survey = ebs_standard_data$survey,
                 survey_design = ebs_standard_data$survey_design,
                 #' Some cruises are shared between the standard and other EBS cruises, so the 
                 #' unique() wrapper is there to remove duplicate cruise records. 
                 cruise = unique(rbind(ebs_standard_data$cruise, ebs_other_data$cruise)),
                 haul = rbind(ebs_standard_data$haul, ebs_other_data$haul),
                 catch = rbind(ebs_standard_data$catch, ebs_other_data$catch),
                 species = ebs_standard_data$species,
                 strata = ebs_standard_data$strata)

# Calculate CPUE and export ---------------------------------------------------
ebs_cpue <- calc_cpue(gapdata = ebs_data) %>%
  select("YEAR", "LATITUDE_DD_START",
         "LONGITUDE_DD_START", "CPUE_NOKM2") %>%
  transmute(cpue = CPUE_NOKM2, 
            year = as.integer(YEAR),
            lat = LATITUDE_DD_START,
            lon = LONGITUDE_DD_START)

write.csv(ebs_cpue, 
          here("data", "pollock_num.csv"), 
          row.names = FALSE)

# Explore numerical CPUE data - check for NAs
nas <- ebs_cpue %>%
  filter(is.na(cpue))

# Compare to DDC data ---------------------------------------------------------
# Read in DDC in numbers (created for ALK)
ddc_alk <- read.csv(here("data", "VAST_ddc_alk_2025.csv")) %>%
  group_by(Year, Lat, Lon) %>%
  summarize(CPUE_num = sum(CPUE_num)) %>%
  # filter to lat/lon in gapindex pull
  filter(Lat %in% ebs_cpue$lat & Lon %in% ebs_cpue$lon) 

# Get annual values for comparison
annual_ddc <- ddc_alk %>%
  group_by(Year) %>%
  summarize(CPUE = sum(CPUE_num)) %>%
  mutate(data = "DDC")
colnames(annual_ddc)[1:2] <- c("year", "cpue")

annual_num <- ebs_cpue %>%
  filter(!is.na(cpue)) %>%  # for now, removing NAs
  group_by(year) %>%
  summarize(cpue = sum(cpue)) %>%
  mutate(data = "raw") %>%
  mutate(cpue = cpue / 100)  # convert to square kilometers

annual_combined <- bind_rows(annual_ddc, annual_num)

ggplot(annual_combined, aes(x = year, y = cpue, color = data, fill = data)) +
  geom_point() +
  geom_line()

# Compare biomass cpue
ebs_biom <- calc_cpue(gapdata = ebs_data) %>%
  select("YEAR", "LATITUDE_DD_START",
         "LONGITUDE_DD_START", "CPUE_KGKM2") %>%
  transmute(cpue = CPUE_KGKM2, 
            year = as.integer(YEAR),
            lat = LATITUDE_DD_START,
            lon = LONGITUDE_DD_START)

ddc_biom <- read.csv(here("data", "VAST_ddc_EBSonly_2025.csv")) %>%
  rename(lat = start_latitude, lon = start_longitude) %>%
  group_by(year, lat, lon) %>%
  summarize(CPUE_kg = sum(ddc_cpue_kg_ha)) %>%
  filter(lat %in% ebs_biom$lat & lon %in% ebs_biom$lon)

annual_ddc_biom <- ddc_biom %>%
  group_by(year) %>%
  summarize(cpue = sum(CPUE_kg)) %>%
  mutate(data = "DDC")

annual_biom <- ebs_biom %>%
  filter(!is.na(cpue)) %>%  
  group_by(year) %>%
  summarize(cpue = sum(cpue)) %>%
  mutate(data = "raw") %>%
  mutate(cpue = cpue / 100)

annual_combined_biom <- bind_rows(annual_ddc_biom, annual_biom)

ggplot(annual_combined_biom, aes(x = year, y = cpue, color = data, fill = data)) +
  geom_point() +
  geom_line()

# Calculate difference between the two data sources and plot together
diff <- bind_rows(
  bind_cols(year = annual_ddc$year, 
            cpue_diff = (annual_ddc$cpue - annual_num$cpue) / annual_num$cpue,
            type = "numbers"),
  bind_cols(year = annual_ddc_biom$year, 
            cpue_diff = (annual_ddc_biom$cpue - annual_biom$cpue) / annual_biom$cpue,
            type = "biomass"),
) 

ggplot(diff, aes(x = year, y = cpue_diff, color = type)) +
  geom_point() +
  geom_line()
