# ======================================================================
# 03 — Visualisation
#
# Histograms of (transformed) predictors, predicted-vs-observed plots
# for each model iteration, univariate scatter plots, and the
# choropleth map of model residuals across UK catchments.
#
# All figures are written to outputs/ as PNGs.
#
# Assumes 01_data_preparation.R and 02_analysis.R have been sourced.
# ======================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggcorrplot)
  library(sf)
  library(tmap)
  library(viridis)
  library(here)
})

out_dir <- here::here("outputs")
dir.create(out_dir, showWarnings = FALSE)

save_plot <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  ggsave(file.path(out_dir, filename), plot = plot,
         width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: outputs/", filename)
}

# --- 1. Distributions of (transformed) predictors --------------------
H1 <- ggplot(main_table_1, aes(x = drain_dens)) +
  geom_histogram(fill = "steelblue", colour = "white", bins = 30) +
  labs(title = "Drainage density", x = "Drainage density",
       y = "Frequency") + theme_minimal()

H2 <- ggplot(main_table_1, aes(x = dpsbar_log)) +
  geom_histogram(fill = "steelblue", colour = "white", bins = 30) +
  labs(title = "Catchment slope (log dpsbar)", x = "log(dpsbar)",
       y = "Frequency") + theme_minimal()

H3 <- ggplot(main_table_1, aes(x = p_mean_boxcox)) +
  geom_histogram(fill = "steelblue", colour = "white", bins = 30) +
  labs(title = "Mean precipitation (Box-Cox)",
       x = "Box-Cox mean precipitation", y = "Frequency") + theme_minimal()

H4 <- ggplot(main_table_1, aes(x = urban_perc_yeojohnson)) +
  geom_histogram(fill = "steelblue", colour = "white", bins = 30) +
  labs(title = "Urban land cover (Yeo-Johnson)",
       x = "Yeo-Johnson urban %", y = "Frequency") + theme_minimal()

distributions <- (H1 | H2) / (H3 | H4)
save_plot(distributions, "01_predictor_distributions.png", width = 10, height = 7)

# --- 2. Predicted vs observed for each model iteration ---------------
pred_obs_plot <- function(data, x, y, title) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(color = "darkred", alpha = 0.7, size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
    labs(title = title,
         x = "Observed high-flow frequency",
         y = "Predicted") +
    theme_minimal()
}

save_plot(
  pred_obs_plot(main_table_1, "freq_q_high_events", "pred_multi",
                "4-predictor model (drainage density, slope, urban %, precipitation)"),
  "02_predicted_vs_observed_initial.png"
)
save_plot(
  pred_obs_plot(main_table_1, "freq_q_high_events", "pred_2ndfinal",
                "6-predictor model (adds slope_fdc, baseflow index)"),
  "03_predicted_vs_observed_improved.png"
)
save_plot(
  pred_obs_plot(main_table_1_trimmed_cooks, "freq_q_high_events", "pred_final",
                "Final model (IQR-trimmed + Cook's-distance filter)"),
  "04_predicted_vs_observed_final.png"
)

# --- 3. Univariate scatter plots -------------------------------------
p_pmean <- ggplot(trimmed_main_table_1,
                  aes(x = p_mean, y = freq_q_high_events)) +
  geom_point(shape = 21, fill = "steelblue", colour = "black",
             size = 1, stroke = 0.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, colour = "darkred", linewidth = 1) +
  labs(title = "Mean precipitation vs high-flow event frequency",
       x = "Mean precipitation (mm)",
       y = "Frequency of high flow events (events/yr)") +
  theme_minimal()
save_plot(p_pmean, "05_scatter_precipitation.png")

p_drain <- ggplot(trimmed_main_table_2,
                  aes(x = drain_dens, y = freq_q_high_events)) +
  geom_point(shape = 21, fill = "purple", colour = "black",
             size = 1, stroke = 0.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, colour = "darkred", linewidth = 1) +
  labs(title = "Drainage density vs high-flow event frequency",
       x = expression("Drainage density (m"^-1*")"),
       y = "Frequency of high flow events (events/yr)") +
  theme_minimal()
save_plot(p_drain, "06_scatter_drainage_density.png")

p_dpsbar <- ggplot(trimmed_main_table_3,
                   aes(x = dpsbar, y = freq_q_high_events)) +
  geom_point(shape = 21, fill = "lightgreen", colour = "black",
             size = 1, stroke = 0.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, colour = "darkred", linewidth = 1) +
  labs(title = "Catchment steepness (dpsbar) vs high-flow event frequency",
       x = "dpsbar (m/km)",
       y = "Frequency of high flow events (events/yr)") +
  theme_minimal()
save_plot(p_dpsbar, "07_scatter_dpsbar.png")

# --- 4. Urbanisation change plots ------------------------------------
p_urban_delta <- ggplot(lmd_clean_4,
                        aes(x = delta_urban_perc, y = delta_mean_discharge)) +
  geom_point(shape = 21, fill = "steelblue", colour = "black",
             size = 1, stroke = 0.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, colour = "darkred", linewidth = 1) +
  labs(title = "Change in urbanisation vs change in mean discharge (1990 -> 2022)",
       x = expression(Delta * " urban land cover (>=4 percentage points)"),
       y = expression(Delta * " mean specific discharge (mm/day)")) +
  theme_minimal()
save_plot(p_urban_delta, "08_urbanisation_delta_discharge.png")

p_urban_area <- ggplot(lmd_clean,
                       aes(x = delta_urban_perc, y = delta_mean_discharge,
                           colour = area)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(method = "lm", colour = "black") +
  scale_colour_viridis_c(option = "plasma", trans = "log10") +
  labs(title = "Urbanisation vs discharge change, coloured by catchment area",
       x = expression(Delta * " urban land cover (%)"),
       y = expression(Delta * " mean discharge"),
       colour = expression(log[10] * " area (km"^2*")")) +
  theme_minimal()
save_plot(p_urban_area, "09_urbanisation_by_area.png")

p_dyu <- ggplot(dyu_clean,
                aes(x = urban, y = discharge, colour = factor(year))) +
  geom_point(alpha = 0.8, size = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_colour_manual(values = c("1990" = "steelblue", "2022" = "darkred"),
                      name = "Year") +
  labs(title = "Discharge vs urban land cover (catchments >10% urban)",
       x = "Urban land cover (%)",
       y = expression("Discharge (m"^3 * "/s)")) +
  theme_minimal()
save_plot(p_dyu, "10_discharge_vs_urban_by_year.png")

# --- 5. Predictor correlation matrix ---------------------------------
corr_plot <- ggcorrplot(cor_matrix, hc.order = TRUE, type = "lower",
                        lab = TRUE, lab_size = 4,
                        colors = c("#AEC6CF", "white", "#FFB347"),
                        title = "Correlation matrix of model predictors",
                        ggtheme = theme_minimal())
save_plot(corr_plot, "11_predictor_correlations.png", width = 7, height = 6)

# --- 6. Choropleth: residuals of the precipitation model -------------
# Spatial structure in residuals indicates predictors the simple model
# is missing — a useful diagnostic for the multivariate extension.
catchments_merged <- catchment_boundaries %>%
  dplyr::left_join(
    main_table_1 %>% sf::st_drop_geometry() %>%
      dplyr::select(ID, resid_m4),
    by = "ID"
  )
st_crs(catchments_merged) <- 27700
catchments_merged <- na.omit(catchments_merged)
uk_boundary <- read_sf(here::here("data", "UK_boundary_NG.shp"))

resid_map <-
  tm_shape(uk_boundary) +
    tm_polygons(fill = "darkolivegreen3") +
  tm_shape(catchments_merged) +
    tm_polygons(
      fill = "resid_m4",
      fill.scale = tm_scale(
        breaks = c(-10, -5, -2, 2, 5, 10),
        values = c("#08306B", "#2171B5", "#6BAED6", "white",
                   "#FC9272", "#FB6A4A", "#CB181D"),
        continuous = TRUE),
      fill.legend = tm_legend(title = "Residuals (p_mean ~ high flow)")
    ) +
  tm_layout(
    inner.margins         = c(0.02, 0.02, 0.02, 0.3),
    frame                 = TRUE,
    legend.outside        = TRUE,
    legend.outside.position = "bottom",
    legend.text.size      = 0.7,
    legend.title.size     = 0.9
  ) +
  tm_scalebar(text.size = 0.5, position = c("RIGHT", "BOTTOM")) +
  tm_compass(type = "8star", position = c("RIGHT", "top"), size = 3)

tmap_save(resid_map, filename = file.path(out_dir, "12_residual_map_uk.png"),
          width = 7, height = 9, dpi = 300)
message("Saved: outputs/12_residual_map_uk.png")

message("\n03_visualisation.R complete. Twelve figures written to outputs/.")
