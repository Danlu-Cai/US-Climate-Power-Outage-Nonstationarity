Data sources and instructions for accessing research datasets. Raw datasets are not included.
# Data Documentation

This directory documents the datasets used in the climate–power outage analysis.

## Primary Data Sources

**EAGLE-I:** County-level electricity outage observations for identifying high-impact outages and characterizing outage patterns.

**CONUS404:** Meteorological data used to calculate temperature indicators, heating and cooling degree days, and extreme temperature events.

## Processed Analysis Datasets

The machine-learning and figure-generation workflows use merged county-level datasets at daily, monthly, and annual resolutions.

Expected filenames include:

- `Digicare_Merged_Daily_2015_2022_CONUS_spatial.csv`
- `Digicare_Merged_Monthly_2015_2022_CONUS_spatial.csv`
- `Digicare_Merged_Yearly_2015_2022_CONUS_spatial.csv`
- `Digicare_Merged_Monthly_2015_2022_CONUS.csv`
- `Digicare_Merged_Yearly_2015_2022_CONUS.csv`

These are filenames referenced by the analysis scripts, not a confirmation that the files are distributed in this repository.

## Data Access

Users should obtain underlying datasets from their authorized sources and follow the corresponding data-use and redistribution requirements.

Processed files should be generated using the documented workflow or obtained through an approved data-sharing mechanism.

## Reproducibility

The data-merging step connecting the original outage and climate datasets to the model-ready files must be documented for end-to-end reproducibility.
