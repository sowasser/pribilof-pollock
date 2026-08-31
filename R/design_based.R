#' Script for pulling design-based (area-weighted) estimates for pollock from 
#' the standard gapindex package.

# pak::pkg_install("afsc-gap-products/gapindex")
library(RODBC)
library(gapindex)
library(dplyr)
library(ggplot2)
library(here)

library(ggsidekick)
theme_set(theme_sleek())

# channel <- gapindex::get_connected(check_access = F)

if (file.exists("Z:/Projects/ConnectToOracle.R")) {
  source("Z:/Projects/ConnectToOracle.R")
} else {
  # For those without a ConnectToOracle file
  channel <- odbcConnect(dsn = "AFSC",
                         uid = rstudioapi::showPrompt(title = "Username",
                                                      message = "Oracle Username",
                                                      default = ""),
                         pwd = rstudioapi::askForPassword("Enter Password"),
                         believeNRows = FALSE)
}

# Pull data
gapindex_data <- gapindex::get_data(
  year_set = c(1982:2025),
  survey_set = "EBS",
  spp_codes = 21740,   
  haul_type = 3,
  abundance_haul = "Y",
  pull_lengths = F,
  channel = channel)

# Fill in zeros and calculate CPUE
cpue <- gapindex::calc_cpue(gapdata = gapindex_data)

# Calculate stratum-level biomass, abundance, and variance
biomass_stratum <- gapindex::calc_biomass_stratum(
  gapdata = gapindex_data,
  cpue = cpue) %>%
  select(STRATUM, YEAR, BIOMASS_MT, BIOMASS_VAR, POPULATION_COUNT, POPULATION_VAR)
write.csv(biomass_stratum, here("biomass_stratum.csv"), row.names = FALSE)

prib_biomass <- biomass_stratum %>%
  filter(STRATUM %in% c(42, 32, 50))

numbers <- ggplot(prib_biomass, aes(x = YEAR, y = (POPULATION_COUNT / 1e9))) +
  geom_line() +
  ylim(0, NA) +
  xlab("") + ylab("Population (Billions)") +
  facet_grid(~ STRATUM, scales = "free")
numbers
ggsave(numbers, filename = here("design_based_numbers.png"), width = 9, height = 3, units = "in")

biomass <- ggplot(prib_biomass, aes(x = YEAR, y = (BIOMASS_MT / 1e6))) +
  geom_line() +
  ylim(0, NA) +
  xlab("") + ylab("Biomass (Mt)") +
  facet_grid(~ STRATUM, scales = "free_y")
biomass

# Comparison to indices from the model
biomass_model <- read.csv(here("results", "biomass", "250kts", "index_biomass.csv")) %>%
  filter(stratum == "all")
