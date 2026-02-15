# SPCOX — Spatial Penalized Cox Models on Triangulated Meshes

R code to fit **penalized spatial Cox proportional hazards** models over 2D triangulated meshes with cross-validated smoothing. The repository includes two empirical case studies (Campi Flegrei and London) and a full simulation pipeline, plus reproducible plotting for meshes, spatial fields, and risk maps.

---

## Contents

- **empirical_campi_flegrei.R** — End-to-end analysis for Campi Flegrei:
  - Builds a mesh from the Campania administrative boundary shapefile.
  - Joins `active.csv` (population at risk) with `trigger_mod.csv` (event times).
  - Fits a penalized spatial PH model (triangulated FEM basis with a sum-to-zero constraint).
  - Produces maps of mesh, event times, estimated spatial field, and risk.
- **empirical_london.R** — End-to-end analysis for London:
  - Uses London ward/boundary shapefiles and the `spatsurv::fstimes` dataset.
  - Fits the same spatial PH model with diurnal terms (Fourier sines/cosines).
  - Generates mesh/time maps, spatial field maps, and risk maps at multiple times of day.
- **empirical_san_francisco.R** — End-to-end analysis for San Francisco:
  - Uses San Franscio data taken from: https://data.sfgov.org/Public-Safety/Fire-Department-and-Emergency-Medical-Services-Dis/nuek-vuh3/about_data [accessed 10/01/2026]
  - Fits the same spatial PH model with diurnal terms (Fourier sines/cosines).
  - Generates mesh/time maps, spatial field maps, and risk maps at multiple times of day.
- **simulation_script.R** — Full simulation pipeline:
  - Simulates data on a “horseshoe” domain (from **fdaPDE**) with varying sample sizes and censoring.
  - Selects the smoothing parameter via K-fold CV of the partial likelihood (parallel option available).
  - Saves a summary table of bias/MSE/coverage and L2 errors.
- **generate_italy_colored.R**, **generate_uk_colored.R** — Country/region basemap helpers.
- **generate_FEM_picture.py** — Helper to render an FEM mesh picture.
- **Shapefiles** — Campania (`SITRC_COMUNI_CAMPANIA.*`) and London (`London_*.*`) boundaries used in the empirical scripts.
- **Data** — `active.csv`, `trigger_mod.csv` for the Campi Flegrei analysis.

> See the repo file list for the exact artifacts and data that ship with the repository. :contentReference[oaicite:0]{index=0}

---

## Installation & Requirements

- **R ≥ 4.x** recommended.
- The scripts auto-install needed CRAN packages on first run (e.g., `ggplot2`, `survival`, `fdaPDE`, `sf`, `dplyr`, `rmapshaper`, `doParallel`, `doRNG`, `Matrix`, `pracma`, `spatsurv`, `viridis`, `ggspatial`, `fields`, `raster`, `lubridate`, `geosphere`, etc.).
- For the optional helper, **Python 3** (only for `generate_FEM_picture.py`).

---

## Quick Start

Clone the repository:
```bash
git clone https://github.com/tedescolor/SPCOX
cd SPCOX
