# Data Processing Scripts

**Project:** Climate–Power Outage Relationships Across the United States  
**Author:** Danlu Cai

This directory contains R scripts for processing county-level power outage observations and climate data used to investigate high-impact power outages and the temporal stability of climate–outage relationships across the contiguous United States.

| Script | Description |
|---|---|
| `01_process_eaglei_high_impact_outages.R` | Processes EAGLE-I 15-minute county-level outage observations (2015–2024), identifies high-impact outage events, and generates hourly, daily, monthly, and annual summaries. |
| `02_aggregate_conus404_degree_days.R` | Processes CONUS404 county-level daily mean air temperature (October 1979–September 2022), calculates cooling and heating degree days, and generates daily, monthly, calendar-year, and water-year summaries. |

Both scripts use configurable input and output directories to support reproducible execution without exposing institution-specific file paths.
