# ======================================================================
# 01 — Data preparation
#
# Project: Hydrological vs geomorphological controls on the frequency
# of high-flow events (fhqe) in UK river catchments.
# Dataset: CAMELS-GB v2 (Coxon et al., 2025), 671 catchments.
#
# This script loads the v2 static attribute tables, the catchment
# boundary shapefile, and the daily hydrometeorological time series
# for every catchment. It produces the analysis-ready tables used by
# 02_analysis.R and 03_visualisation.R.
#
# Outputs in the global environment:
#   camel_table  — full attribute table (all attribute CSVs joined on ID)
#   main_table   — analysis-ready subset of variables
#   main_table_1 — 3-sigma trimmed table for multivariate modelling
#                  (with Box-Cox / Yeo-Johnson / log transforms applied)
#   lmd_clean    — change-in-discharge frame (1990 -> 2022, trimmed)
#   dyu_clean    — long-format discharge/urban frame for the interaction model
# ======================================================================

# --- Packages ---------------------------------------------------------
# Run install_packages.R once before sourcing this script.
suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(lubridate)
  library(sf)
  library(tools)
  library(car)        # powerTransform, bcPower, yjPower
})

# --- Paths ------------------------------------------------------------
# `here()` resolves paths relative to the project root, so this works on
# any machine without editing.
data_dir  <- here::here("data")
daily_dir <- here::here("data", "Daily")

# --- Load static attribute tables -------------------------------------
climate          <- read.csv(file.path(data_dir, "CAMELS_GB_v2_climatic_attributes.csv"))
hydrogeology     <- read.csv(file.path(data_dir, "CAMELS_GB_v2_hydrogeology_attributes.csv"))
hydrologic       <- read.csv(file.path(data_dir, "CAMELS_GB_v2_hydrologic_attributes.csv"))
hydrometry       <- read.csv(file.path(data_dir, "CAMELS_GB_v2_hydrometry_attributes.csv"))
landcover        <- read.csv(file.path(data_dir, "CAMELS_GB_v2_landcover_attributes.csv"))
soil             <- read.csv(file.path(data_dir, "CAMELS_GB_v2_soil_attributes.csv"))
topographic      <- read.csv(file.path(data_dir, "CAMELS_GB_v2_topographic_attributes.csv"))
human_influence  <- read.csv(file.path(data_dir, "CAMELS_GB_v2_humaninfluences_attributes.csv"))
groundwater_well <- read.csv(file.path(data_dir, "CAMELS_GB_v2_groundwaterwell_attributes.csv"))

catchment_boundaries <- read_sf(file.path(data_dir, "CAMELS_GB_catchment_boundaries.shp"))
drain_dens           <- read.csv(file.path(data_dir, "camels_OS_drainagedens.csv"))

# --- Merge all attribute tables into a single sf table ----------------
# `all = TRUE` keeps every catchment even if individual attribute files
# disagree on which gauges are present. The boundary geometry comes last
# so the result is an sf object keyed on ID.
camel_table <- merge(climate, human_influence, all = TRUE)
camel_table <- merge(camel_table, hydrogeology,  all = TRUE)
camel_table <- merge(camel_table, hydrologic,    all = TRUE)
camel_table <- merge(camel_table, hydrometry,    all = TRUE)
camel_table <- merge(camel_table, landcover,     all = TRUE)
camel_table <- merge(camel_table, soil,          all = TRUE)
camel_table <- merge(camel_table, topographic,   all = TRUE)
camel_table <- merge(catchment_boundaries, camel_table,
                     by.x = "ID", by.y = "gauge_id")
camel_table$drain_dens <- drain_dens$drainage_dens

# --- Build the focused analysis table (main_table) --------------------
# A trimmed table holding only the variables used downstream. Kept as an
# sf object so that catchment geometries propagate to plots and maps.
main_table <- data.frame(ID = camel_table$ID, dpsbar = camel_table$dpsbar)
main_table <- merge(catchment_boundaries, main_table)

main_table$freq_highq    <- camel_table$high_q_freq
main_table$high_q_dur    <- camel_table$high_q_dur
main_table$cond_hypres   <- soil$conductivity_hypres
main_table$tawc          <- soil$tawc
main_table$p_mean        <- climate$p_mean
main_table$temp          <- camel_table$temperature      # may be NA in some rows
main_table$drain_dens    <- camel_table$drain_dens
main_table$freq_q_high_events <- main_table$freq_highq / main_table$high_q_dur
main_table$ID_STRING     <- NULL

# Replace NAs with zero for downstream consistency. This is a coarse
# choice (it biases means and inflates the n of "zero" catchments);
# revisit if you extend the analysis.
main_table[is.na(main_table)] <- 0

# --- Process daily hydromet time series -------------------------------
# CAMELS-GB ships one daily file per catchment under data/Daily/.
# We need two summaries:
#   (1) per-catchment long-run means (for joining onto main_table)
#   (2) per-catchment means for 1990 and 2022 specifically (for the
#       urbanisation change analysis)

daily_files <- list.files(daily_dir, pattern = "\\.csv$", full.names = TRUE)

# (1) Long-run means across the full period of record per catchment.
mean_table <- daily_files %>%
  map_df(function(file) {
    data <- read_csv(file, show_col_types = FALSE)
    means <- data %>%
      summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
    mutate(means, Site = file_path_sans_ext(basename(file)))
  }) %>%
  dplyr::select(Site, everything()) %>%
  mutate(
    Site_ID = as.numeric(str_extract(Site, "(?<=_)\\d{3,6}(?=_)"))
  ) %>%
  arrange(Site_ID)
mean_table$Site <- NULL

# Pull the variables we care about onto main_table.
main_table$mean_temp           <- mean_table$temperature_haduk
main_table$mean_discharge_spec <- mean_table$discharge_spec
main_table$arid                <- camel_table$aridity
main_table$urban_perc          <- camel_table$urban_perc_2022
main_table[is.na(main_table)]  <- 0

# Coordinates as "lat, lon" for convenience.
main_table <- main_table %>%
  mutate(coordinates = paste(camel_table$gauge_lat, camel_table$gauge_lon, sep = ", "))

# (2) Per-catchment, per-year means for 1990 and 2022.
process_catchment <- function(file) {
  df <- read.csv(file)
  df %>%
    mutate(date = as.Date(date),
           year = year(date)) %>%
    filter(year %in% c(1990, 2022)) %>%
    group_by(year) %>%
    summarise(mean_discharge = mean(discharge_spec, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(catchment = tools::file_path_sans_ext(basename(file)))
}

mean_discharge <- map_df(daily_files, process_catchment) %>%
  mutate(Site_ID = as.numeric(str_extract(catchment, "(?<=_)\\d{3,6}(?=_)")),
         ID      = Site_ID)

# --- Long format: year-by-year discharge with matching urban % --------
discharge_joined <- mean_discharge %>%
  left_join(camel_table, by = "ID") %>%
  mutate(urban = case_when(
    year == 1990 ~ urban_perc_1990,
    year == 2022 ~ urban_perc_2022,
    TRUE         ~ NA_real_
  ))

discharge_year_urban <- data.frame(
  ID        = discharge_joined$ID,
  year      = discharge_joined$year,
  urban     = discharge_joined$urban,
  discharge = discharge_joined$mean_discharge
)
discharge_year_urban[is.na(discharge_year_urban)] <- 0

# IQR-trim per year, then restrict to catchments with >=10% urban cover
# so the interaction effect isn't dominated by near-natural catchments.
dyu_clean <- discharge_year_urban %>%
  group_by(year) %>%
  mutate(
    Q1          = quantile(discharge, 0.25, na.rm = TRUE),
    Q3          = quantile(discharge, 0.75, na.rm = TRUE),
    IQR         = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(discharge >= lower_bound & discharge <= upper_bound) %>%
  ungroup() %>%
  filter(urban >= 10)

# --- Wide format: 1990 vs 2022 deltas per catchment ------------------
long_mean_discharge <- mean_discharge %>%
  pivot_wider(names_from = year, values_from = mean_discharge,
              names_prefix = "mean_") %>%
  mutate(Site_ID = as.numeric(str_extract(catchment, "(?<=_)\\d{3,6}(?=_)"))) %>%
  arrange(Site_ID)
long_mean_discharge$catchment            <- NULL
long_mean_discharge$freq_q_high_events   <- main_table$freq_q_high_events
long_mean_discharge$urban_perc_1990      <- landcover$urban_perc_1990
long_mean_discharge$urban_perc_2022      <- landcover$urban_perc_2022
long_mean_discharge[is.na(long_mean_discharge)] <- 0

# Trim upper-tail discharge outliers in each year independently.
lmd_clean <- long_mean_discharge
upper_1990 <- quantile(lmd_clean$mean_1990, 0.75, na.rm = TRUE) +
              1.5 * IQR(lmd_clean$mean_1990, na.rm = TRUE)
upper_2022 <- quantile(lmd_clean$mean_2022, 0.75, na.rm = TRUE) +
              1.5 * IQR(lmd_clean$mean_2022, na.rm = TRUE)
lmd_clean <- lmd_clean %>%
  dplyr::filter(mean_1990 <= upper_1990, mean_2022 <= upper_2022) %>%
  mutate(
    delta_urban_perc     = urban_perc_2022 - urban_perc_1990,
    delta_mean_discharge = mean_2022 - mean_1990
  )

# Attach catchment area and drop the smallest 25% to focus the
# urban-change signal on larger, less variable catchments.
camel_table$area <- camel_table$area.x
area_df          <- data.frame(area = camel_table$area, ID = camel_table$ID)
lmd_clean        <- merge(lmd_clean, area_df) %>%
  filter(area >= quantile(area, 0.25, na.rm = TRUE))

# Subset to catchments with >=4 percentage-point urban growth, used for
# focused correlation tests in 02_analysis.R.
lmd_clean_4 <- subset(lmd_clean, delta_urban_perc >= 4)

# --- Trim main_table for univariate analyses --------------------------
# Use 3-sigma bounds for each predictor of interest. The bounds are
# computed once here so 02_analysis.R can use them throughout.

bounds_3sigma <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  c(low = m - 3 * s, high = m + 3 * s)
}

b_freq      <- bounds_3sigma(main_table$freq_highq)
b_pmean     <- bounds_3sigma(main_table$p_mean)
b_drain     <- bounds_3sigma(main_table$drain_dens)
b_dpsbar    <- bounds_3sigma(main_table$dpsbar)
b_urban     <- bounds_3sigma(main_table$urban_perc)

trimmed_main_table <- subset(main_table,
  freq_highq >= b_freq["low"]  & freq_highq <= b_freq["high"])

trimmed_main_table_1 <- subset(trimmed_main_table,
  p_mean >= b_pmean["low"] & p_mean <= b_pmean["high"])
trimmed_main_table_1$freq_q_high_events <-
  trimmed_main_table_1$freq_highq / trimmed_main_table_1$high_q_dur

trimmed_main_table_2 <- subset(trimmed_main_table,
  drain_dens >= b_drain["low"] & drain_dens <= b_drain["high"])

trimmed_main_table_3 <- subset(trimmed_main_table,
  dpsbar >= b_dpsbar["low"] & dpsbar <= b_dpsbar["high"])

trimmed_main_table_4 <- subset(trimmed_main_table,
  urban_perc >= b_urban["low"] & urban_perc <= b_urban["high"])

# --- Build the multivariate-model table -------------------------------
# Apply all four 3-sigma trims simultaneously, then add transformed
# versions of the right-skewed predictors.
main_table_1 <- main_table %>%
  filter(
    p_mean     >= b_pmean["low"]  & p_mean     <= b_pmean["high"],
    drain_dens >= b_drain["low"]  & drain_dens <= b_drain["high"],
    dpsbar     >= b_dpsbar["low"] & dpsbar     <= b_dpsbar["high"],
    urban_perc >= b_urban["low"]  & urban_perc <= b_urban["high"]
  )

# Box-Cox for p_mean (strictly positive).
boxcox_pmean <- powerTransform(lm(p_mean ~ 1, data = main_table_1))
main_table_1$p_mean_boxcox <- bcPower(main_table_1$p_mean, lambda = boxcox_pmean$lambda)

# Yeo-Johnson for urban_perc (handles zeros).
yeo_urban <- powerTransform(lm(urban_perc ~ 1, data = main_table_1), family = "yjPower")
main_table_1$urban_perc_yeojohnson <- yjPower(main_table_1$urban_perc, lambda = yeo_urban$lambda)
main_table_1$urban_perc_log        <- log(main_table_1$urban_perc + 1)

# Log dpsbar (right-skewed).
main_table_1$dpsbar_log <- log(main_table_1$dpsbar)

# Bring in two additional candidate predictors from camel_table.
main_table_1 <- main_table_1 %>%
  dplyr::left_join(
    camel_table %>%
      sf::st_drop_geometry() %>%
      dplyr::select(ID, slope_fdc, baseflow_index),
    by = "ID"
  )

main_table_1[is.na(main_table_1)] <- 0

message("01_data_preparation.R complete.")
message("Catchments in main_table:          ", nrow(main_table))
message("Catchments in trimmed main_table_1:", nrow(main_table_1))
message("Catchments in lmd_clean:           ", nrow(lmd_clean))
message("Rows in dyu_clean:                 ", nrow(dyu_clean))
