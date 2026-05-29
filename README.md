# Shoreline Model Comparison (3-Year Segments)

This repository contains the MATLAB code and data used to run the Yates09 and ShoreFor equilibrium shoreline models in 3-year moving windows at Hasaki Beach, Japan.

## Files

| File | Description |
|------|-------------|
| `Models_3y.m` | Main script (self-contained; all functions included) |
| `Hasaki_Shore_ref1.40.csv` | Pier-based shoreline observations |
| `FT_mean_shoreline.csv` | Alongshore-averaged (FT mean) shoreline position |
| `Hasaki_Wave_JRA55.csv` | Offshore wave data (JRA-55 reanalysis) |

## Model Calibration

- **Yates09**: Five free parameters (*a*, *b*, *c_acr*, *c_ero*, *Y_ini*) are optimized by minimizing RMSE between modeled and observed shoreline positions. The optimization uses `fmincon` with the Sequential Quadratic Programming (SQP) algorithm, wrapped in `MultiStart` with 5 random initial points to avoid local minima.
- **ShoreFor**: The beach memory parameter *phi* is selected by grid search over a predefined set of candidate values. Accretion and erosion coefficients are then obtained via linear regression with a two-step trend correction.

## Data Sources

- **Single-transect shoreline (pier)**: Hasaki Oceanographical Research Station (HORS) pier survey data, available at <https://pari.mpat.go.jp/bdhome/Hasaki/>.
- **Multi-transect shoreline (FT mean)**: Alongshore-averaged shoreline position derived from multi-transect survey data, available at <https://doi.org/10.5281/zenodo.20439203>.
- **Wave data**: JRA-55-based WAVEWATCH III wave reanalysis, available at <https://doi.org/10.5281/zenodo.15622122>.

## Requirements

- MATLAB R2016b or later
- Optimization Toolbox (`fmincon`, `MultiStart`)
- Statistics and Machine Learning Toolbox (`signrank`; optional — the script falls back gracefully if unavailable)

## Usage

1. Place all files in the same directory.
2. Open MATLAB and set the current folder to that directory.
3. Run `Models_3y.m`.

## Outputs

- `Model_<model>_pier_<period>.csv` — Daily modeled shoreline position for each 3-year window
- `model_performance_3yr.csv` — Summary table of R² and RMSE for all periods and models
- Time series comparison figures (one per model)

