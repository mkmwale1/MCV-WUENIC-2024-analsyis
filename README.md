README.md for MCV-WUENIC-2024-Analysis
Overview
This repository contains the R script for analyzing measles-containing vaccine (MCV) coverage using the WHO/UNICEF Estimates of National Immunization Coverage (WUENIC) 2024 revision. The analysis focuses on trends, inequities, and projections in low- and middle-income countries (LMICs), as described in the associated manuscript.
The script processes publicly available data from WUENIC and World Bank sources to generate descriptive statistics, statistical tests, Gini coefficients, k-means clustering, and projections to 2030.
Requirements

R Version: 4.4.1 or higher
Required Packages:

rio (1.2.3)
tidyverse (2.0.0)
readxl (1.4.3)
sf (1.0-16)
rnaturalearth (1.0.1)
countrycode (1.6.0)
ineq (0.2-13)
broom (1.0.6)
ggalluvial (0.12.5)
cluster (2.1.6)
factoextra (1.0.7)
ggplot2 (3.5.1)
corrplot (0.92)



Install packages using:
rinstall.packages(c("rio", "tidyverse", "readxl", "sf", "rnaturalearth", "countrycode", "ineq", "broom", "ggalluvial", "cluster", "factoextra", "ggplot2", "corrplot"))
Data Sources

WUENIC Data: Download from WHO Immunization Coverage (2024 revision).
World Bank Data: Fragility status and income classifications from World Bank FCS List and Income Classifications.
UN Population Data: Surviving infants from UN World Population Prospects.

Place downloaded files in a data/ directory or update paths in the script accordingly.
Usage

Clone the repository:
bashgit clone https://github.com/mkmwale1/MCV-WUENIC-2024-analysis.git

Navigate to the directory:
bashcd MCV-WUENIC-2024-analysis

Run the script in R:
rsource("script.R")  # Replace with your script filename


The script performs data cleaning, imputation, statistical analyses, clustering, and generates figures/tables as described in the manuscript.
Output

Figures: Saved as PNG/PDF in an output/figures/ directory (e.g., alluvial diagrams, choropleth maps).
Tables: Exported as CSV or included in console output (e.g., fragility gaps, regional burdens).
Projections: 2030 MCV1 estimates with standard errors.

Replication Notes

Ensure data files are up-to-date from sources.
Sensitivity analyses (e.g., imputation methods) are included; adjust parameters as needed.
For questions, refer to the manuscript or open an issue.

License
This project is licensed under the MIT License - see the LICENSE file for details.
Contact
Author: Moses Mwale
Affiliation: World Health Organization, Lusaka, Zambia
Email: mwalem@who.int
Last Updated: October 15, 2025
