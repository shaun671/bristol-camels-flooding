# Data

The analysis uses the **CAMELS-GB v2** dataset — a large-sample hydrological
dataset for Great Britain published by UKCEH. The dataset is freely available
but **not redistributed in this repository** due to its size and licensing.

## Downloading CAMELS-GB v2

1. Go to the CAMELS-GB v2 record on the EIDC catalogue:
   <https://catalogue.ceh.ac.uk/documents/db37d4b5-6cb1-4f25-bdc8-e7f81ffc2cbc>
   *(Search "CAMELS-GB v2" on the EIDC if the link has moved.)*
2. Register for a free EIDC account if you don't already have one.
3. Download the full archive and extract it.
4. Place the contents into the `data/` directory of this repository so the
   layout matches the structure below.

## Expected file layout

```
data/
├── CAMELS_GB_v2_climatic_attributes.csv
├── CAMELS_GB_v2_hydrogeology_attributes.csv
├── CAMELS_GB_v2_hydrologic_attributes.csv
├── CAMELS_GB_v2_hydrometry_attributes.csv
├── CAMELS_GB_v2_landcover_attributes.csv
├── CAMELS_GB_v2_soil_attributes.csv
├── CAMELS_GB_v2_topographic_attributes.csv
├── CAMELS_GB_v2_humaninfluences_attributes.csv
├── CAMELS_GB_v2_groundwaterwell_attributes.csv
├── CAMELS_GB_catchment_boundaries.shp     (plus .dbf, .shx, .prj)
├── camels_OS_drainagedens.csv             (derived from OS data — see note)
├── UK_boundary_NG.shp                     (plus .dbf, .shx, .prj)
└── Daily/
    ├── CAMELS_GB_hydromet_timeseries_<ID>_<period>.csv
    ├── ...
```

## Notes on auxiliary data

- **`camels_OS_drainagedens.csv`** is a per-catchment drainage density layer
  derived from Ordnance Survey hydrographic data. If you don't have it, the
  drainage-density analysis (`drain_dens`) can be commented out, or replaced
  with an alternative source such as the EU-Hydro stream network.
- **`UK_boundary_NG.shp`** is a national outline used only for the
  choropleth map in `R/03_visualisation.R`. Any UK boundary shapefile in
  British National Grid (EPSG:27700) will work — e.g. the ONS Countries
  (December 2022) BGC dataset from the Open Geography Portal.

## Citation

If you use the dataset, please cite the v2 release:

> Coxon, G., Zheng, Y., Barbedo, R., Fileni, F., Fowler, H., Fry, M.,
> Green, A., Harfoot, H., Lewis, E., Qiu, X., Salwey, S., and Wendt, D.
> (2025). CAMELS-GB v2: hydrometeorological time series and landscape
> attributes for 671 catchments in Great Britain. *EGU General
> Assembly 2025*, EGU25-4371.
> <https://doi.org/10.5194/egusphere-egu25-4371>

The v2 release builds on the original CAMELS-GB paper:

> Coxon, G., Addor, N., Bloomfield, J. P., et al. (2020). CAMELS-GB:
> hydrometeorological time series and landscape attributes for 671
> catchments in Great Britain. *Earth System Science Data*, 12(4),
> 2459–2483. <https://doi.org/10.5194/essd-12-2459-2020>
