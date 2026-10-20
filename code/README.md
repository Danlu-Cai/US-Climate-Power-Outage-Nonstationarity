# Research Analysis Scripts

This directory contains the R scripts used to analyze climate–power outage relationships across the contiguous United States.

## 01. EAGLE-I Outage Data Processing

**Script:** `01_process_eaglei_high_impact_outages.R`

**Purpose:** Process EAGLE-I outage observations from 2015–2024 to construct county-level outage metrics and identify high-impact outage days.

**High-impact outage criteria:**

- A continuous county-level outage event lasting more than 8 consecutive clock hours, including events crossing midnight.
- Daily customer-outage-hours (COH) exceeding 5,000.
- For events spanning multiple days, each calendar day must independently satisfy the COH threshold.

**Outputs:** Hourly, daily, monthly and annual outage summaries.

**Reproducibility:** The script requires EAGLE-I input data and appropriate R package dependencies. Timestamp handling, missing observations and event continuity should be validated against the source data before using the outputs for publication.
