
# Data Processing and Analysis Scripts

**Project:** Climate–Power Outage Relationships Across the United States  
**Author:** Danlu Cai

This directory contains R scripts for processing county-level power outage observations, climate data, and machine-learning analyses used to investigate high-impact power outages and climate–power outage relationships across the contiguous United States.

## Scripts

### 01_process_eaglei_high_impact_outages.R
Processes EAGLE-I 15-minute county-level power outage observations (2015–2024), identifies high-impact outage events, and generates hourly, daily, monthly, and annual summaries.

High-impact county-days are defined by daily customer outage hours (COH) exceeding 5,000 and overlap with a continuous county-level outage event lasting more than 8 hours.

### 02_aggregate_conus404_degree_days.R
Processes CONUS404 county-level daily mean air temperature data (October 1979–September 2022) to calculate cooling degree days (CDD) and heating degree days (HDD).

Generates daily, monthly, calendar-year, and water-year summaries using a base temperature of 18.3°C.

### 03_identify_conus404_temperature_extremes.R
Identifies extreme heat and cold events using county-specific historical temperature percentiles from CONUS404 (October 1979–September 2022).

Uses a two-tier extreme-event identification approach and generates daily, monthly, and annual summaries.

### 04_national_ml_attribution_grouped_cv.R
Evaluates associations between environmental predictors and county-level power outage metrics using Random Forest and XGBoost at daily, monthly, and annual resolutions (2015–2022).

The analysis includes:
- Five-fold county-grouped cross-validation.
- Three predictor configurations: unweighted degree days, weighted degree days, and their combination.
- Customer outage hours and outage-hour outcomes.
- Model performance metrics (R² and MAE).
- Feature importance and annual diagnostic analyses.
- CSV exports and diagnostic figures.

County-grouped cross-validation evaluates spatial generalization to held-out counties. It does not directly evaluate temporal transferability to unseen years.

## Dependencies

Install the required R packages from the repository root:

```r
source("requirements.R")
```

## Data and Reproducibility

Input and output directories can be configured using environment variables.

The scripts do not require institution-specific paths when the corresponding environment variables are configured.

The merged datasets used by Script 04 must be prepared separately. They are not generated directly by Scripts 01–03.

Some underlying datasets may require separate access or preprocessing.
