# ======================================================================
# 02 — Statistical analysis
#
# All correlation tests, univariate and multivariate linear models,
# model diagnostics, and the model-iteration loop (extra predictors
# from the Spearman residual scan, then Cook's-distance trimming).
#
# Headline result of the analysis: baseflow_index dominates the final
# model with ~78% relative importance — contradicting the hypothesis
# that geomorphological factors (Dd, dpsbar) would be the strongest
# controls on high-flow event frequency.
#
# Assumes 01_data_preparation.R has already been sourced.
# ======================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lmtest)
  library(car)        # Anova(), vif()
  library(relaimpo)   # calc.relimp()
})

# --- 1. Urbanisation & discharge change (1990 -> 2022) ---------------
# Focus on catchments with >=4 pp increase in urban cover.
cat("\n--- Urbanisation change correlations -----------------------\n")
print(cor.test(lmd_clean_4$delta_urban_perc, lmd_clean_4$delta_mean_discharge))
print(cor.test(lmd_clean_4$delta_urban_perc, lmd_clean_4$freq_q_high_events))

# Full sample with year as a factor: interaction between urban% and year.
discharge_model <- lm(discharge ~ urban * year, data = dyu_clean)
cat("\n--- discharge ~ urban * year (Anova, type III) ---------------\n")
print(Anova(discharge_model, type = "III"))
cat("\n--- discharge ~ urban * year (summary) -----------------------\n")
print(summary(discharge_model))

# --- 2. Univariate correlations with high-flow-event frequency -------
cat("\n--- Univariate correlations ---------------------------------\n")

cat("\nMean precipitation -- frequency of high flow events:\n")
print(cor.test(trimmed_main_table_1$p_mean, trimmed_main_table_1$freq_q_high_events))

cat("\nDrainage density -- frequency of high flow events:\n")
print(cor.test(trimmed_main_table_2$drain_dens, trimmed_main_table_2$freq_q_high_events))

cat("\nCatchment slope (dpsbar) -- frequency of high flow events:\n")
print(cor.test(trimmed_main_table_3$dpsbar, trimmed_main_table_3$freq_q_high_events))

cat("\nUrban cover -- frequency of high flow events:\n")
print(cor.test(trimmed_main_table_4$urban_perc, trimmed_main_table_4$freq_q_high_events))

# Same as simple linear models (kept for reporting slope estimates).
cat("\nSlope estimates (single-predictor lm fits):\n")
print(lm(freq_q_high_events ~ p_mean,     data = trimmed_main_table_1))
print(lm(freq_q_high_events ~ dpsbar,     data = trimmed_main_table_3))

# --- 3. Initial multivariate model (4 predictors) --------------------
m1 <- lm(freq_q_high_events ~ drain_dens,            data = main_table_1)
m2 <- lm(freq_q_high_events ~ dpsbar_log,            data = main_table_1)
m3 <- lm(freq_q_high_events ~ urban_perc_yeojohnson, data = main_table_1)
m4 <- lm(freq_q_high_events ~ p_mean_boxcox,         data = main_table_1)

m_full_scaled <- lm(
  freq_q_high_events ~
    scale(drain_dens) + scale(dpsbar_log) +
    scale(urban_perc_yeojohnson) + scale(p_mean_boxcox),
  data = main_table_1
)

cat("\n--- 4-predictor multivariate model (scaled) -----------------\n")
print(summary(m_full_scaled))

cat("\nBreusch-Pagan heteroscedasticity test:\n")
print(bptest(m_full_scaled))

cat("\nRelative importance (lmg method):\n")
print(calc.relimp(m_full_scaled, type = "lmg", rela = TRUE))

# Store predictions for plotting (used in 03_visualisation.R).
main_table_1$pred_multi <- predict(m_full_scaled)
main_table_1$resid_m4   <- residuals(m4)

# --- 4. Spearman scan of residuals against all numeric attributes ----
# Drives the choice of additional predictors below.
cat("\n--- Spearman scan: residuals of p_mean model vs all attrs ---\n")
camel_table$freq_q_high_events <- main_table$freq_q_high_events
boxcox_pmean_full <- powerTransform(lm(p_mean.x ~ 1, data = camel_table))
camel_table$p_mean_boxcox <- bcPower(camel_table$p_mean.x,
                                     lambda = boxcox_pmean_full$lambda)
camel_table[is.na(camel_table)] <- 0
lm_camel <- lm(freq_q_high_events ~ p_mean_boxcox,
               data = camel_table, na.action = na.fail)
camel_table$resid_p_mean <- residuals(lm_camel)

predictors <- camel_table %>%
  sf::st_drop_geometry() %>%
  select_if(is.numeric) %>%
  mutate(across(everything(), as.numeric))

spearman_results <- lapply(names(predictors), function(v) {
  test <- cor.test(camel_table$resid_p_mean, predictors[[v]], method = "spearman")
  data.frame(variable = v, rho = test$estimate, p_value = test$p.value)
}) %>% bind_rows() %>%
  arrange(desc(abs(rho)))

cat("\nTop 15 variables by |rho|:\n")
print(head(spearman_results, 15))

# --- 5. Improved 6-predictor model -----------------------------------
# Add slope_fdc and baseflow_index identified above.
m_full_scaled_2 <- lm(
  freq_q_high_events ~
    scale(drain_dens) + scale(dpsbar_log) +
    scale(urban_perc_yeojohnson) + scale(p_mean_boxcox) +
    scale(slope_fdc) + scale(baseflow_index),
  data = main_table_1
)

cat("\n--- 6-predictor multivariate model --------------------------\n")
print(summary(m_full_scaled_2))
cat("\nCondition number (kappa):", kappa(m_full_scaled_2), "\n")
cat("VIFs:\n"); print(vif(m_full_scaled_2))
cat("\nBreusch-Pagan:\n"); print(bptest(m_full_scaled_2))

main_table_1$pred_multi_2 <- predict(m_full_scaled_2)

# --- 6. Trim to interquartile range of the response ------------------
# Removes both heavy upper-tail leverage and very-low-frequency
# catchments to address residual heteroscedasticity.
freq_Q1 <- quantile(main_table_1$freq_q_high_events, 0.25, na.rm = TRUE)
freq_Q3 <- quantile(main_table_1$freq_q_high_events, 0.75, na.rm = TRUE)
main_table_1_trimmed <- subset(main_table_1,
  freq_q_high_events >= freq_Q1 & freq_q_high_events <= freq_Q3)

m_full_scaled_2_trimmed <- lm(
  freq_q_high_events ~
    scale(drain_dens) + scale(dpsbar_log) +
    scale(urban_perc_yeojohnson) + scale(p_mean_boxcox) +
    scale(slope_fdc) + scale(baseflow_index),
  data = main_table_1_trimmed
)

main_table_1_trimmed$pred_multi_2 <- predict(m_full_scaled_2_trimmed)
cat("\n--- IQR-trimmed model: Breusch-Pagan ------------------------\n")
print(bptest(m_full_scaled_2_trimmed))

# --- 7. Cook's-distance trim then final model ------------------------
main_table_1_trimmed$cooks <- cooks.distance(m_full_scaled_2_trimmed)
cooks_limit <- 3 * mean(main_table_1_trimmed$cooks)
main_table_1_trimmed_cooks <- subset(main_table_1_trimmed, cooks <= cooks_limit)

m_final <- lm(
  freq_q_high_events ~
    scale(drain_dens) + scale(dpsbar_log) +
    scale(urban_perc_yeojohnson) + scale(p_mean_boxcox) +
    scale(slope_fdc) + scale(baseflow_index),
  data = main_table_1_trimmed_cooks
)

main_table_1_trimmed_cooks$pred_final <- predict(m_final)
main_table_1$pred_2ndfinal <- predict(m_full_scaled_2)

cat("\n--- FINAL MODEL ---------------------------------------------\n")
print(summary(m_final))
cat("\nVIFs:\n"); print(vif(m_final))
cat("\nBreusch-Pagan:\n"); print(bptest(m_final))
cat("\nRelative importance (lmg):\n")
m_final_results <- calc.relimp(m_final, type = "lmg", rela = TRUE)
print(m_final_results)

# --- 8. Correlation matrix of predictors -----------------------------
vars <- c("drain_dens", "dpsbar_log", "urban_perc_yeojohnson",
          "p_mean_boxcox", "slope_fdc", "baseflow_index")
numeric_vars <- main_table_1[, vars]
numeric_vars$geometry <- NULL
cor_matrix <- cor(numeric_vars, use = "complete.obs")
cat("\n--- Predictor correlation matrix ----------------------------\n")
print(round(cor_matrix, 2))

message("\n02_analysis.R complete.")
