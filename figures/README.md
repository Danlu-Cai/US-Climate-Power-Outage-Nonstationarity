Scripts and outputs for manuscript figures

# Figure Generation Scripts

**Project:** Climate–Power Outage Relationships Across the United States  
**Author:** Danlu Cai

## Figure 1: Spatiotemporal Patterns of High-Impact Power Outages

**Script:** `Figure_1_outage_trends.R`

Generates Figure 1 panels using EAGLE-I county-level daily outage summaries (2015–2024).

- **Figure 1a:** Spatial trends in county-level outage hours.
- **Figure 1b:** Regional heterogeneity in outage-hour trends.
- **Figure 1c:** National annual outage duration and customer outage hours.
- **Figure 1d(i):** Monthly changes in high-impact outage county-days.
- **Figure 1d(ii):** Seasonal contributions to total outage hours.

Outputs include individual high-resolution PNG panels and supporting CSV tables.

County trend estimates use linear regression with monthly controls and heteroskedasticity/autocorrelation-consistent (HAC) standard errors.
