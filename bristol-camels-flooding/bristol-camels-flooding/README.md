# Hydrological vs geomorphological controls on high-flow frequency in UK river catchments

A reproducible quantitative-geography analysis using **CAMELS-GB v2** (Coxon
et al., 2025; 671 catchments) to test which catchment attributes most
strongly determine the frequency of high-flow events. Built as an assessed
report at the University of Bristol.

---

## TL;DR

- **Question.** What is the relative importance of hydrological and
  geomorphological factors in determining the frequency of high-flow
  events in UK river catchments?
- **Hypothesis (disproven).** Geomorphological factors (drainage density,
  catchment slope) would dominate climatic and hydrological factors,
  because they directly govern catchment response time.
- **Headline finding.** The hypothesis was incorrect. In the only
  homoscedastic model produced (Breusch-Pagan p = 0.050), **baseflow
  index alone accounted for ~78% of the explained variance** in
  high-flow event frequency — by far the strongest single predictor.
  Drainage density's apparent dominance in the 4-predictor model
  collapsed from a relative importance of 64.5% to 2.8% (a ~2,300%
  reduction) once the two additional hydrological predictors
  (baseflow index, slope of the flow-duration curve) were included.
- **Urbanisation.** No robust effect on mean discharge between 1990 and
  2022 after appropriate filtering (p = 0.643 after restricting the
  sample as in Over et al., 2025) — contradicting US-based studies that
  report a 1.2–3.3% increase in peak flow per 1% increase in
  impervious cover.

> *Caveat the author flags in the discussion: the linear model is
> homoscedastic only after IQR and Cook's-distance trimming, so it
> does not generalise to catchments with infrequent but extended
> floods. A polynomial regression is suggested as the natural next
> step.*

---

## Background

High-flow events drive fluvial flood risk, channel geomorphology, and
aquatic ecology. From 1970 to 2019, flooding caused 31% of all economic
losses globally (Caretta et al., 2022), and Bates et al. (2023) estimate
52.3% of UK flood costs come specifically from fluvial events. The
relative importance of climatic forcing (precipitation), geomorphological
controls (drainage density `Dd`, mean drainage-path slope `dpsbar`), and
land-use (`%U`) is a long-standing open question — particularly in the UK,
where the Flood Estimation Handbook approach uses a similar descriptor
selection.

The hypothesis tested here, drawing on Pallard et al. (2009) and
De Niel & Willems (2019), was that geomorphological descriptors would
dominate, because they directly govern hydrological lag time and
catchment response. The result reported below contradicts that
hypothesis.

## Data

**CAMELS-GB v2** (Coxon et al., 2025): 671 catchments in Great Britain
with static catchment attributes across nine themes (climate, hydrology,
hydrogeology, hydrometry, land cover, soil, topography, human influences,
groundwater wells), plus daily hydrometeorological time series and
catchment boundary polygons in British National Grid (EPSG:27700).

Auxiliary data:
- **Drainage density (`Dd`)** derived per catchment from
  Ordnance Survey Open Rivers (Ordnance Survey, 2025) using QGIS.
- **UK national outline** for the choropleth map.

CAMELS-GB v2 is not redistributed in this repository — see
[`data/README.md`](data/README.md) for download instructions.

## Methods

### Pipeline

The analysis runs in three logical stages, one R script each:

| Stage | Script | What it does |
|------:|:-------|:-------------|
| 1 | `R/01_data_preparation.R` | Loads and merges the nine attribute tables and catchment boundaries; processes daily files into both long-run and year-specific (1990, 2022) means; applies 3-sigma and IQR trimming; computes Box-Cox and Yeo-Johnson transformations for skewed predictors. |
| 2 | `R/02_analysis.R` | Univariate correlations and linear fits; the 4-predictor multivariate model; a Spearman residual scan to identify additional predictors; the improved 6-predictor model; IQR + Cook's-distance trimming to produce the final model; full diagnostic suite. |
| 3 | `R/03_visualisation.R` | Predictor histograms, predicted-vs-observed plots for each model, univariate scatter plots, urbanisation-change plots, the predictor correlation matrix, and the UK residual choropleth. All figures saved to `outputs/`. |

### Response variable

`fhqe = high_q_freq / high_q_dur` — the frequency of high-flow events per
year, normalised by mean event duration. Normalising by duration reduces
sensitivity to short-term flash flooding (which is largely pluvial and
unrelated to antecedent catchment conditions, per Ran et al., 2022) and
better represents catchment-driven extended high-flow regimes.

### Predictors

| Symbol | Variable | Transformation | Justification |
|:-------|:---------|:---------------|:--------------|
| `Dd`   | Drainage density (m⁻¹) | none | approximately normal |
| `dpsbar_log` | log(dpsbar) | natural log | right-skewed |
| `pbc`  | Mean precipitation | Box-Cox | strictly positive, right-skewed |
| `%Uyj` | Urban land cover, 2022 | Yeo-Johnson | contains zeros |
| `slope_fdc` | Slope of flow-duration curve | none | added after Spearman scan |
| `baseflow_index` | Baseflow index | none | added after Spearman scan |

### Modelling sequence

1. **Univariate.** Pearson correlations and simple linear fits for each
   of the four primary descriptors, with 3-sigma trimming on each
   predictor in turn.
2. **Initial 4-predictor model** on the 3-sigma-trimmed sample using
   scaled predictors. Diagnostics revealed heteroscedasticity and
   precipitation was statistically insignificant.
3. **Spearman residual scan.** Residuals of the univariate
   `fhqe ~ pbc` model were tested against every numeric attribute in
   CAMELS-GB v2. The two strongest absolute correlations
   (`baseflow_index` ρ = −0.82, `slope_fdc` ρ = 0.56) were added to
   the model as proxies for the missing variance.
4. **6-predictor model**, followed by IQR trimming on the response and
   Cook's-distance trimming (`cooks > 3 × mean(cooks)`) to produce the
   final, homoscedastic model.

### Diagnostics

- Breusch–Pagan test for heteroscedasticity (final model
  p = 0.050, borderline homoscedastic).
- Variance Inflation Factors (max VIF in final model = 2.16, well below
  O'Brien's (2007) threshold of 4).
- Predictor correlation matrix (greatest pairwise correlation
  ρ = −0.68 between `baseflow_index` and `slope_fdc` — approaching but
  below Kim's (2019) collinearity threshold of 0.8–0.9).
- `relaimpo::calc.relimp` with the `lmg` decomposition to apportion R²
  among correlated predictors.

### Urbanisation analysis

For the 1990 vs 2022 urbanisation question, daily discharge was aggregated
to per-catchment annual means and joined to the corresponding year's urban
land cover. Two formulations were tested:

1. A **delta model** (`lmd_clean`): ΔU vs ΔQ on the wide-format frame,
   restricted to catchments larger than the 25th-percentile area. A
   second test restricted to catchments with ≥4 pp urban change,
   following Over et al. (2025).
2. A **long-format interaction model** (`dyu_clean`):
   `discharge ~ urban * year`, with year as a factor and the sample
   restricted to catchments with ≥10% urban cover.

## Key results

### Final model (6 predictors, IQR + Cook's-distance trimmed)

| Variable | Coefficient | Relative importance (%) | p-value |
|:---------|------------:|------------------------:|--------:|
| `baseflow_index`        | −1.549   | **77.95** | < 2 × 10⁻¹⁶ |
| `slope_fdc`             | −0.1357  | 11.24     | 0.069       |
| `dpsbar_log`            | −0.4798  | 5.32      | 5.2 × 10⁻⁸  |
| `drain_dens`            |  3.8231  | 2.78      | < 2 × 10⁻¹⁶ |
| `urban_perc_yeojohnson` |  0.1262  | 2.01      | 0.115       |
| `p_mean_boxcox`         | −0.04853 | 0.70      | 0.495       |

Breusch–Pagan p = 0.050; max VIF = 2.16.

### Spearman residual scan (top variables vs residuals of `fhqe ~ pbc`)

| Variable | ρ | p-value |
|:---------|--:|--------:|
| `baseflow_index` | −0.822 | 8.4 × 10⁻¹⁶⁶ |
| `low_q_freq`     |  0.655 | 1.5 × 10⁻⁸³  |
| `slope_fdc`      |  0.557 | 6.7 × 10⁻⁵⁶  |
| `bulkdens_50`    | −0.471 | 2.2 × 10⁻³⁸  |
| `tawc`           |  0.456 | 1.1 × 10⁻³⁵  |

This step is what produced the central finding of the analysis: the
strongest predictor of variance unexplained by precipitation is *not*
geomorphology, but a hydrological signature (baseflow index) reflecting
how much of a catchment's discharge is groundwater-supported.

### Spatial pattern

The choropleth of univariate `fhqe ~ pbc` residuals
(`outputs/12_residual_map_uk.png`) shows weak spatial structure overall
but small clusters of underprediction in northern England and the south
east, consistent with the conclusion that a single climatic predictor is
inadequate.

## Limitations (from the discussion)

- The final model is homoscedastic only after IQR and Cook's-distance
  trimming — it does not generalise to catchments with infrequent but
  extended floods.
- `p_mean_boxcox` and `urban_perc_yeojohnson` remain statistically
  insignificant in the final model and likely drive heteroscedasticity in
  the outer quartiles. They could legitimately be dropped to avoid
  inflating R².
- The hypothesis test relied on a linear specification; the residual
  distributions suggest a polynomial regression would be more appropriate.
- Several relationships contradict known literature (notably the
  near-zero precipitation effect and the absent urban–discharge
  relationship). This raises questions about cross-region generalisability
  of US-based studies (Over et al., 2016, 2025; Blum et al., 2020) to UK
  catchments, and possibly about specific attributes in CAMELS-GB v2.
- An instantaneous-peak-flow measure (e.g. percentile-based threshold,
  per Bartens et al., 2024) would likely better represent extended
  high-flow regimes than the daily-mean-derived `fhqe`.

## Reproducing the analysis

### Requirements

- R ≥ 4.2.
- ~2–3 GB free disk space for CAMELS-GB v2.
- A C/C++ toolchain for compiled packages (Rtools on Windows,
  Xcode CLT on macOS, build-essential on Linux).

### Steps

```r
# 1. Clone the repo
# git clone https://github.com/<your-username>/bristol-camels-flooding.git
# cd bristol-camels-flooding

# 2. Open the project in RStudio (or set working directory to the repo root)

# 3. Install required packages (one-off)
source("install_packages.R")

# 4. Download CAMELS-GB v2 into data/ as described in data/README.md

# 5. Run the full pipeline
source("main.R")
```

The full pipeline takes a few minutes. Figures land in `outputs/`; model
summaries, diagnostic tests, and relative-importance scores are printed to
the console.

## Repository layout

```
bristol-camels-flooding/
├── README.md                    <- you are here
├── LICENSE                      <- MIT
├── .gitignore                   <- excludes data and generated outputs
├── main.R                       <- runs the full pipeline
├── install_packages.R           <- one-off package install helper
├── R/
│   ├── 01_data_preparation.R    <- load, merge, transform
│   ├── 02_analysis.R            <- correlations, models, diagnostics
│   └── 03_visualisation.R       <- ggplot and tmap outputs
├── data/
│   └── README.md                <- where to get CAMELS-GB v2
└── outputs/
    └── (12 PNG figures, generated)
```

## Skills demonstrated

- **R / tidyverse / sf** — multi-table joins, spatial-data handling,
  functional-style processing of 671 per-catchment CSVs via
  `purrr::map_df`.
- **Quantitative hydrology** — large-sample-hydrology methodology
  applied to a recently released national dataset, separation of
  climatic from catchment-intrinsic controls, baseflow-index
  interpretation.
- **Regression methodology** — appropriate transformations (Box-Cox,
  Yeo-Johnson, log), heteroscedasticity diagnostics, multicollinearity
  checks (VIF, condition number), relative-importance scoring via the
  `lmg` decomposition, Spearman residual scanning, Cook's-distance
  outlier handling.
- **Reproducible workflow** — scripted pipeline, `here()` for portable
  paths, explicit data-acquisition instructions, no hard-coded local
  paths.
- **Scientific cartography** — catchment-scale residual choropleth in
  `tmap` with diverging palette, custom breaks, north arrow and
  scale bar.

## Citation

If you use this code, please cite the CAMELS-GB v2 dataset:

> Coxon, G., Zheng, Y., Barbedo, R., Fileni, F., Fowler, H., Fry, M.,
> Green, A., Harfoot, H., Lewis, E., Qiu, X., Salwey, S., and Wendt, D.
> (2025). CAMELS-GB v2: hydrometeorological time series and landscape
> attributes for 671 catchments in Great Britain. *EGU General
> Assembly 2025*, EGU25-4371.
> <https://doi.org/10.5194/egusphere-egu25-4371>

A reading list of references cited above (Pallard, De Niel & Willems,
Over, Blum, Bates, Ran, Bartens, O'Brien, Kim, and others) is included
in the original report.

## Author

**Shaun Turner** — 2nd-year Geography undergraduate, University of
Bristol. Contact: shaunturner680@gmail.com

## License

[MIT](LICENSE).
