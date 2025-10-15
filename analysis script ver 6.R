# ============================================================
# Title: Persistent Gaps in Measles Immunisation Coverage in LMICs (WUENIC 2024)
# Author: Moses Mwale
# Created: 2025-08-04
# Last edited: 2025-10-13
# Purpose: Reproducible analysis of MCV1/MCV2 coverage, inequality, fragility gaps,
#          clustering, and IA2030 tracking for low- and middle‑income countries (LMICs).
# Notes:
#   - Uses 2019–2024 WUENIC measles data.
#   - 2024 surviving‑infants values are substituted with 2023 where missing, else regional/global medians.
#   - MCV2 set to 0% for countries that have not introduced the second dose (verified list below).
#   - Linear projections to 2030 are per‑country OLS fits.
#   - Outputs (tables/figures) are written to the "outputs/" directory.
#   - UK English spelling is used throughout.
# Licence: MIT (add a LICENCE file in the repository for clarity).
# ============================================================

# ------------------------------------------------------------
# 0) Setup and packages
# ------------------------------------------------------------

# Recommended: use renv for reproducible environments in the repo.
# install.packages("pacman") # if needed
suppressPackageStartupMessages({
  pacman::p_load(
    here, rio, readxl, tidyverse, janitor,
    sf, rnaturalearth, rnaturalearthdata, countrycode,
    ineq, broom, ggplot2, ggalluvial, cluster, factoextra,
    corrplot, patchwork, gt, webshot, scales
  )
})

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20251013)

dir.create(here::here("outputs"), showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 11))

# ------------------------------------------------------------
# 1) Configuration: file paths (edit for your repo layout)
# ------------------------------------------------------------

PATH_WUENIC   <- here::here("data", "MCV DATA.xlsx")
PATH_WB_INC   <- here::here("data", "World bank CLASS_2025_07_02.xlsx")
PATH_POP      <- here::here("data", "live births.xlsx")

# ------------------------------------------------------------
# 2) Helper functions
# ------------------------------------------------------------

pivot_coverage <- function(data, year_cols, value_name) {
  data %>%
    select(unicef_region, iso3, country, all_of(year_cols)) %>%
    pivot_longer(cols = all_of(year_cols), names_to = "Year", values_to = value_name) %>%
    mutate(Year = as.numeric(Year))
}

dropout_flag <- function(x) {
  dplyr::case_when(
    is.na(x)   ~ "Missing",
    x < 0      ~ "<0% (check data)",
    x > 60     ~ ">60% (check data)",
    TRUE       ~ "OK"
  )
}

.add_sizes <- function(df, tidy_row, metric_label) {
  sz <- df %>% count(Fragility_Status) %>%
    pivot_wider(names_from = Fragility_Status, values_from = n, values_fill = 0)
  names(sz) <- gsub("\\s|-", "_", names(sz))
  tidy_row %>% mutate(metric = metric_label) %>% bind_cols(sz)
}

check_iso_mismatches <- function(df_iso3, world_sf) {
  anti_join(df_iso3 %>% distinct(iso3),
            world_sf %>% st_drop_geometry() %>% select(iso_a3) %>% distinct(),
            by = c("iso3" = "iso_a3"))
}

ensure_phantomjs <- function() {
  if (!webshot::is_phantomjs_installed()) {
    message("Installing PhantomJS (once-off) for gt -> PNG export…")
    webshot::install_phantomjs()
  }
}

# ------------------------------------------------------------
# 3) Load data
# ------------------------------------------------------------

wuenic_mcv1 <- readxl::read_excel(PATH_WUENIC, sheet = "MCV1") %>%
  mutate(iso3 = toupper(trimws(iso3)))

wuenic_mcv2 <- readxl::read_excel(PATH_WUENIC, sheet = "MCV2") %>%
  mutate(iso3 = toupper(trimws(iso3)))

wuenic_global <- readxl::read_excel(PATH_WUENIC, sheet = "regional_global_mcv") %>%
  rename(Year = year) %>% mutate(Year = as.numeric(Year))

wuenic_mcv1_long <- pivot_coverage(wuenic_mcv1, as.character(2019:2024), "Coverage_MCV1")
wuenic_mcv2_long <- pivot_coverage(wuenic_mcv2, as.character(2019:2024), "Coverage_MCV2")

wuenic <- left_join(
  wuenic_mcv1_long, wuenic_mcv2_long,
  by = c("unicef_region", "iso3", "country", "Year")
)

# World Bank income groups and regions
income_raw <- readxl::read_excel(PATH_WB_INC, sheet = "composition") %>%
  mutate(WB_Country_Code = toupper(trimws(WB_Country_Code)))

income <- income_raw %>%
  filter(WB_Group_Name %in% c("Low income", "Lower middle income", "Upper middle income",
                              "Middle income", "High income")) %>%
  mutate(
    Income_Level = WB_Group_Name,
    Middle_or_High = case_when(
      Income_Level %in% c("Lower middle income", "Upper middle income", "Middle income") ~ "Middle income",
      Income_Level == "High income" ~ "High income",
      TRUE ~ Income_Level
    )
  ) %>%
  arrange(WB_Country_Code, Income_Level) %>%
  group_by(WB_Country_Code) %>% slice(1) %>% ungroup() %>%
  select(ISO3 = WB_Country_Code, Country = WB_Country_Name, Income_Level, Middle_or_High)

wb_regions <- income_raw %>%
  filter(str_detect(WB_Group_Name, "excluding high income")) %>%
  select(ISO3 = WB_Country_Code, WB_Region = WB_Group_Name) %>% distinct()

# Fragility list (World Bank FCS 2024)
fragile_list <- c(
  "Afghanistan","Burkina Faso","Cameroon","Central African Republic",
  "Democratic Republic of the Congo","Ethiopia","Haiti","Iraq","Lebanon","Mali",
  "Mozambique","Myanmar","Niger","Nigeria","Somalia","South Sudan","Sudan",
  "Syria","Ukraine","Yemen","Burundi","Chad","Comoros","Republic of the Congo",
  "Eritrea","Guinea-Bissau","Libya"
)
fragile_iso <- countrycode::countrycode(fragile_list, "country.name", "iso3c", warn = FALSE)

# Population: Surviving infants
pop_data <- readxl::read_excel(PATH_POP) %>%
  filter(Type == "Country/Area", Year %in% 2019:2023) %>%
  select(ISO3 = `ISO3 Alpha-code`, Year,
         Surviving_Infants = `Live Births Surviving to Age 1 (thousands)`) %>%
  mutate(ISO3 = toupper(trimws(ISO3)), Surviving_Infants = Surviving_Infants * 1000)

pop_data_2023 <- pop_data %>% filter(Year == 2023) %>%
  select(ISO3, Surviving_Infants_2023 = Surviving_Infants)

# ------------------------------------------------------------
# 4) Harmonise, enrich, and impute
# ------------------------------------------------------------

wuenic <- wuenic %>%
  left_join(income,    by = c("iso3" = "ISO3")) %>%
  left_join(wb_regions,by = c("iso3" = "ISO3")) %>%
  mutate(LMIC_Flag = if_else(Income_Level %in% c("Low income","Lower middle income","Upper middle income"), "LMIC", "HIC")) %>%
  mutate(Fragility_Status = if_else(iso3 %in% fragile_iso, "Fragile", "Non-Fragile")) %>%
  left_join(pop_data, by = c("iso3" = "ISO3", "Year")) %>%
  left_join(pop_data_2023, by = c("iso3" = "ISO3")) %>%
  mutate(
    used_2023_for_2024 = (Year == 2024 & is.na(Surviving_Infants) & !is.na(Surviving_Infants_2023)),
    Surviving_Infants   = if_else(Year == 2024 & is.na(Surviving_Infants), Surviving_Infants_2023, Surviving_Infants)
  )

# Regional median imputation for any remaining missing surviving infants
wuenic <- wuenic %>%
  group_by(unicef_region, Year) %>%
  mutate(Surviving_Infants = if_else(is.na(Surviving_Infants),
                                     median(Surviving_Infants, na.rm = TRUE),
                                     Surviving_Infants)) %>%
  ungroup()

# Countries that have not introduced MCV2 (set to 0% if NA)
non_intro_countries <- c("BEN", "CAF", "GAB", "SSD")

wuenic <- wuenic %>%
  mutate(
    Coverage_MCV2 = if_else(iso3 %in% non_intro_countries & is.na(Coverage_MCV2), 0, Coverage_MCV2),
    Zero_Dose_Prevalence = if_else(!is.na(Coverage_MCV1), 1 - Coverage_MCV1/100, NA_real_),
    Dropout_Rate = if_else(!is.na(Coverage_MCV1) & !is.na(Coverage_MCV2) & Coverage_MCV1 > 0,
                           ((Coverage_MCV1 - Coverage_MCV2)/Coverage_MCV1) * 100, NA_real_),
    Unvaccinated_MCV1 = Zero_Dose_Prevalence * Surviving_Infants,
    Unvaccinated_MCV2 = if_else(!is.na(Coverage_MCV2), (1 - Coverage_MCV2/100)*Surviving_Infants, NA_real_)
  ) %>%
  select(-Surviving_Infants_2023)

imputed_count <- wuenic %>%
  filter(Year == 2024) %>%
  summarise(
    used_2023_substitution = sum(used_2023_for_2024, na.rm = TRUE),
    remaining_imputed_via_region_median = sum(is.na(Zero_Dose_Prevalence) & !is.na(Surviving_Infants), na.rm = TRUE)
  )
readr::write_csv(imputed_count, here::here("outputs", "Supp_Table_S5b_imputed_count_2024.csv"))

# 2024 snapshot with QA flags
data_2024 <- wuenic %>%
  filter(Year == 2024) %>%
  mutate(
    Income_Level = factor(Income_Level, levels = c("Low income","Lower middle income","Upper middle income","Middle income","High income")),
    MCV2_Missing = if_else(is.na(Coverage_MCV2), "MCV2: No Data", "MCV2: Reported"),
    Dropout_Outlier = dropout_flag(Dropout_Rate)
  )

# ------------------------------------------------------------
# 5) Descriptive summaries & IA2030 tracking
# ------------------------------------------------------------

summary_lmics <- data_2024 %>%
  filter(LMIC_Flag == "LMIC") %>%
  group_by(Income_Level, unicef_region, Fragility_Status) %>%
  summarise(
    Mean_MCV1 = mean(Coverage_MCV1, na.rm = TRUE),
    Median_MCV1 = median(Coverage_MCV1, na.rm = TRUE),
    Mean_MCV2 = mean(Coverage_MCV2, na.rm = TRUE),
    Total_Unvaccinated = sum(Unvaccinated_MCV1, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(summary_lmics, here::here("outputs","Supp_Table_S1_summary_lmics_2024.csv"))

zd_region <- data_2024 %>%
  filter(LMIC_Flag == "LMIC") %>%
  group_by(unicef_region) %>%
  summarise(
    Mean_ZD = mean(Zero_Dose_Prevalence, na.rm = TRUE) * 100,
    Total_Unvaccinated = sum(Unvaccinated_MCV1, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(desc(Total_Unvaccinated))
readr::write_csv(zd_region, here::here("outputs","Supp_Table_S2_zero_dose_by_region_2024.csv"))

ia2030 <- wuenic %>%
  filter(LMIC_Flag == "LMIC", Year %in% c(2019, 2024)) %>%
  mutate(Meets_95 = if_else(!is.na(Coverage_MCV1) & !is.na(Coverage_MCV2) &
                              Coverage_MCV1 >= 95 & Coverage_MCV2 >= 95, 1L, 0L)) %>%
  group_by(Year) %>% summarise(Share_Meeting_95 = mean(Meets_95, na.rm = TRUE) * 100, .groups = "drop")
readr::write_csv(ia2030, here::here("outputs","Supp_Table_S3_ia2030_95.csv"))

# Publication-friendly bar
p_ia2030_bar <- ggplot(ia2030, aes(x = factor(Year), y = Share_Meeting_95, fill = factor(Year))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Share_Meeting_95)), vjust = -0.5, size = 4.2) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "LMICs Meeting ≥95% Coverage for MCV1 and MCV2", x = NULL, y = "Countries (%)") +
  ylim(0, max(15, max(ia2030$Share_Meeting_95, na.rm = TRUE) + 5))

ggsave(here::here("outputs","fig_ia2030_meeting95_bar.png"), p_ia2030_bar, width = 6, height = 4, dpi = 300)

# ------------------------------------------------------------
# 6) Maps: MCV1/MCV2/Dropout/Zero-dose (2024)
# ------------------------------------------------------------

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Diagnostics for ISO3 mismatches (optional)
iso_miss <- check_iso_mismatches(data_2024 %>% select(iso3), world)
if (nrow(iso_miss) > 0) readr::write_csv(iso_miss, here::here("outputs","_diagnostic_iso_mismatches.csv"))

map_df <- function(var, title, palette = "viridis", lims = c(0,100), legend = title) {
  left_join(world, data_2024 %>% select(iso3, !!sym(var)) %>% rename(iso_a3 = iso3), by = "iso_a3") %>%
    ggplot() +
    geom_sf(aes(fill = .data[[var]]), colour = "white", size = 0.1) +
    scale_fill_viridis_c(option = palette, na.value = "grey90", limits = lims, name = legend) +
    labs(title = title, subtitle = "Grey indicates no data") +
    theme(legend.box = "vertical")
}

p_map_mcv1   <- map_df("Coverage_MCV1",  "Global MCV1 Coverage, 2024", palette = "viridis")
p_map_mcv2   <- map_df("Coverage_MCV2",  "MCV2 Coverage, 2024",     palette = "viridis")
p_map_dropout<- map_df("Dropout_Rate",    "Dropout: MCV1 → MCV2, 2024", palette = "plasma")
map_zero <- left_join(world, data_2024 %>% transmute(iso3, ZeroDosePct = Zero_Dose_Prevalence*100) %>% rename(iso_a3 = iso3), by = "iso_a3")
p_map_zero <- ggplot(map_zero) +
  geom_sf(aes(fill = ZeroDosePct), colour = "white", size = 0.1) +
  scale_fill_viridis_c(option = "cividis", na.value = "grey90", limits = c(0, 100), name = "Zero-dose (%)") +
  labs(title = "Zero-dose (MCV1), 2024", subtitle = "Grey indicates no data") +
  theme(legend.box = "vertical")

# Save maps
walk2(list(p_map_mcv1, p_map_mcv2, p_map_dropout, p_map_zero),
      c("map_mcv1_2024.png","map_mcv2_2024.png","map_dropout_2024.png","map_zero_dose_2024.png"),
      ~ ggsave(here::here("outputs", .y), .x, width = 10, height = 6, dpi = 300))

# ------------------------------------------------------------
# 7) Trends and linear projections to 2030
# ------------------------------------------------------------

trends <- wuenic %>%
  filter(Income_Level %in% c("Low income","Lower middle income","Upper middle income","Middle income","High income")) %>%
  group_by(Year, Income_Level) %>%
  summarise(
    Mean_MCV1 = mean(Coverage_MCV1, na.rm = TRUE),
    Mean_Dropout = mean(Dropout_Rate, na.rm = TRUE),
    Total_Unvaccinated = sum(Unvaccinated_MCV1, na.rm = TRUE), .groups = "drop"
  )

p_trends <- ggplot(trends, aes(Year, Mean_MCV1, colour = Income_Level)) +
  geom_line(size = 1.1) +
  labs(title = "MCV1 Trend by Income Group (2019–2024)", y = "Coverage (%)", x = "Year") +
  scale_color_brewer(palette = "Set2")

ggsave(here::here("outputs","trend_mcv1_income.png"), p_trends, width = 8, height = 5, dpi = 300)

# Per‑country OLS projection to 2030
proj_coefs <- wuenic %>%
  filter(Year %in% 2019:2024, !is.na(Coverage_MCV1)) %>%
  group_by(iso3) %>%
  nest() %>%
  mutate(
    fit  = map(data, ~ lm(Coverage_MCV1 ~ Year, data = .x)),
    pred = map(fit,  ~ predict(.x, newdata = tibble(Year = 2030), se.fit = TRUE))
  ) %>%
  transmute(
    iso3,
    proj_2030 = map_dbl(pred, ~ .x$fit[1]) %>% pmin(100) %>% pmax(0),
    proj_2030_se = map_dbl(pred, ~ .x$se.fit[1])
  )

data_2024 <- data_2024 %>% left_join(proj_coefs, by = "iso3")

map_proj <- left_join(world, data_2024 %>% select(iso3, proj_2030) %>% rename(iso_a3 = iso3), by = "iso_a3")
p_map_proj <- ggplot(map_proj) +
  geom_sf(aes(fill = proj_2030), colour = "white", size = 0.1) +
  scale_fill_viridis_c(option = "magma", na.value = "grey90", limits = c(0, 100), name = "Projected MCV1 (%)") +
  labs(title = "Projected MCV1 in 2030 (Linear Model)") +
  theme(legend.box = "vertical")

ggsave(here::here("outputs","map_projected_mcv1_2030.png"), p_map_proj, width = 10, height = 6, dpi = 300)

# ------------------------------------------------------------
# 8) Group comparisons (Fragile vs Non‑Fragile) + Income strata
# ------------------------------------------------------------

# Welch t‑tests (and Wilcoxon for robustness) in LMICs
lmics_2024 <- data_2024 %>% filter(LMIC_Flag == "LMIC")

# Overall
tt_mcv1  <- t.test(Coverage_MCV1 ~ Fragility_Status, data = lmics_2024, var.equal = FALSE)
tt_mcv2  <- t.test(Coverage_MCV2 ~ Fragility_Status, data = lmics_2024 %>% filter(!is.na(Coverage_MCV2)), var.equal = FALSE)
wx_mcv1  <- wilcox.test(Coverage_MCV1 ~ Fragility_Status, data = lmics_2024, exact = FALSE)
wx_mcv2  <- wilcox.test(Coverage_MCV2 ~ Fragility_Status, data = lmics_2024 %>% filter(!is.na(Coverage_MCV2)), exact = FALSE)

overall_gap <- tibble(
  Stratification = "Overall (LMICs)",
  Metric = c("MCV1","MCV2"),
  `Mean Difference (Fragile − Non‑fragile)` = c(broom::tidy(tt_mcv1)$estimate, broom::tidy(tt_mcv2)$estimate),
  `95% CI Low`  = c(broom::tidy(tt_mcv1)$conf.low,  broom::tidy(tt_mcv2)$conf.low),
  `95% CI High` = c(broom::tidy(tt_mcv1)$conf.high, broom::tidy(tt_mcv2)$conf.high),
  `p (Welch t)` = c(broom::tidy(tt_mcv1)$p.value,   broom::tidy(tt_mcv2)$p.value),
  `p (Wilcoxon)`= c(broom::tidy(wx_mcv1)$p.value,   broom::tidy(wx_mcv2)$p.value)
)
readr::write_csv(overall_gap, here::here("outputs","Table3_fragile_nonfragile_overall.csv"))

# Income‑level stratified tests (skip safely where only one level exists)

tt_by_income <- data_2024 %>%
  filter(!is.na(Income_Level)) %>%
  group_by(Income_Level) %>%
  group_modify(~{
    df  <- .x %>% filter(!is.na(Fragility_Status))
    rows <- list()
    if (n_distinct(df$Fragility_Status) == 2) {
      t1 <- broom::tidy(t.test(Coverage_MCV1 ~ Fragility_Status, data = df, var.equal = FALSE))
      rows[[1]] <- .add_sizes(df, t1, "MCV1")
      df2 <- df %>% filter(!is.na(Coverage_MCV2))
      if (nrow(df2) > 1 && n_distinct(df2$Fragility_Status) == 2) {
        t2 <- broom::tidy(t.test(Coverage_MCV2 ~ Fragility_Status, data = df2, var.equal = FALSE))
        rows[[2]] <- .add_sizes(df2, t2, "MCV2")
      }
    }
    if (length(rows) == 0) tibble(estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_, metric = "MCV1") else bind_rows(rows)
  }) %>%
  ungroup()

# Build Table 3 (CSV + gt image)

table3_fragility <- bind_rows(
  overall_gap %>% transmute(
    Stratification, Metric,
    `Mean Difference (Fragile - Non-Fragile)` = `Mean Difference (Fragile − Non‑fragile)`,
    `95% CI Low`, `95% CI High`, `p-value` = `p (Welch t)`
  ),
  tt_by_income %>% transmute(
    Stratification = as.character(Income_Level),
    Metric = metric,
    `Mean Difference (Fragile - Non-Fragile)` = estimate,
    `95% CI Low` = conf.low,
    `95% CI High` = conf.high,
    `p-value` = p.value
  )
) %>%
  filter(Metric %in% c("MCV1","MCV2")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)), `p-value` = scales::pvalue(`p-value`))

readr::write_csv(table3_fragility, here::here("outputs", "Table3_fragility_gaps_2024.csv"))

ensure_phantomjs()

tbl3 <- table3_fragility %>%
  gt() %>%
  tab_header(
    title = md("**Fragility Gaps in Measles Coverage (LMICs, 2024)**"),
    subtitle = md("Mean differences (fragile − non‑fragile) with 95% CIs; Welch tests")
  ) %>%
  cols_label(
    Stratification ~ "Stratification",
    Metric ~ "Metric",
    `Mean Difference (Fragile - Non-Fragile)` ~ "Mean Diff (pp)",
    `95% CI Low` ~ "95% CI Low",
    `95% CI High` ~ "95% CI High",
    `p-value` ~ "p-value"
  ) %>%
  tab_options(table.font.size = px(10)) %>%
  opt_table_lines()

gtsave(tbl3, here::here("outputs", "Table3_fragility_gaps_2024.html"))
webshot(here::here("outputs", "Table3_fragility_gaps_2024.html"),
        here::here("outputs", "Table3_fragility_gaps_2024.png"), vwidth = 2000, vheight = 1000, zoom = 2)

# ------------------------------------------------------------
# 9) Alluvial: MCV1 coverage group transitions (2019 → 2024)
# ------------------------------------------------------------

coverage_groups <- wuenic %>%
  filter(Year %in% c(2019, 2024)) %>%
  mutate(Coverage_Group = case_when(
    Coverage_MCV1 >= 90 ~ "≥90%",
    Coverage_MCV1 >= 70 ~ "70–89%",
    TRUE ~ "<70%"
  )) %>%
  select(iso3, Year, Coverage_Group) %>%
  filter(!is.na(Coverage_Group)) %>%
  group_by(iso3, Year) %>% slice(1) %>% ungroup() %>%
  group_by(iso3) %>% filter(n() == 2) %>% ungroup()

alluvial_data <- coverage_groups %>%
  pivot_wider(names_from = Year, values_from = Coverage_Group, names_prefix = "Y") %>%
  filter(!is.na(Y2019), !is.na(Y2024)) %>%
  count(Y2019, Y2024)

alluvial_data_lab <- alluvial_data %>%
  mutate(
    Y2019_lab = factor(paste0("2019: ", Y2019), levels = c("2019: <70%", "2019: 70–89%", "2019: ≥90%")),
    Y2024_lab = factor(paste0("2024: ", Y2024), levels = c("2024: <70%", "2024: 70–89%", "2024: ≥90%")),
    legend_2019 = factor(Y2019, levels = c("<70%", "70–89%", "≥90%"))
  )

p_alluvial <- ggplot(alluvial_data_lab, aes(axis1 = Y2019_lab, axis2 = Y2024_lab, y = n)) +
  ggalluvial::geom_alluvium(aes(fill = legend_2019), width = 1/12, alpha = 0.85, colour = NA) +
  ggalluvial::geom_stratum(width = 1/12, fill = "grey90", colour = "black") +
  scale_x_discrete(limits = c("2019", "2024"), labels = c("2019", "2024"), expand = c(.02, .02)) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(name = "2019 Coverage Group",
                    breaks = c("<70%","70–89%","≥90%"),
                    values = c("<70%"  = "#E69F00", "70–89%"= "#56B4E9", "≥90%"= "#0072B2")) +
  labs(title = "Shift in MCV1 Coverage Groups (2019 → 2024)", y = "Number of Countries", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.ticks = element_blank())

ggsave(here::here("outputs","alluvial_group_shift.png"), p_alluvial, width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------
# 10) Inequality and correlations
# ------------------------------------------------------------

gini_mcv1 <- ineq::Gini(data_2024$Coverage_MCV1, na.rm = TRUE)
readr::write_lines(paste0("Gini (MCV1, 2024): ", round(gini_mcv1, 3)), here::here("outputs","_gini_mcv1_2024.txt"))

num_vars <- data_2024 %>% select(Coverage_MCV1, Unvaccinated_MCV1, Surviving_Infants) %>% drop_na()
cor_matrix <- cor(num_vars)
readr::write_csv(as.data.frame(cor_matrix) %>% rownames_to_column("var"), here::here("outputs","Supp_Table_S4_cor_matrix.csv"))

# ------------------------------------------------------------
# 11) Clustering (LMICs): features, k choice, labels, and plots
# ------------------------------------------------------------

features <- wuenic %>%
  filter(LMIC_Flag == "LMIC") %>%
  group_by(iso3) %>%
  summarise(
    MCV1_2024 = mean(Coverage_MCV1[Year == 2024], na.rm = TRUE),
    Mean_MCV1 = mean(Coverage_MCV1, na.rm = TRUE),
    Mean_Dropout = mean(Dropout_Rate, na.rm = TRUE),
    Change_MCV1 = mean(Coverage_MCV1[Year == 2024], na.rm = TRUE) - mean(Coverage_MCV1[Year == 2019], na.rm = TRUE),
    Mean_Unvaccinated = mean(Unvaccinated_MCV1, na.rm = TRUE), .groups = "drop"
  ) %>%
  left_join(wuenic %>% filter(Year == 2024) %>% select(iso3, Fragility_Status, country), by = "iso3") %>%
  mutate(Pct_Fragile = if_else(Fragility_Status == "Fragile", 100, 0))

# Choose k via elbow + silhouette
set.seed(123)
wss <- map_dbl(1:10, ~ kmeans(scale(features %>% select(-iso3, -Fragility_Status, -country)), centers = .x, nstart = 25)$tot.withinss)
elbow_data <- tibble(k = 1:10, wss = wss)

ggsave(here::here("outputs","cluster_elbow.png"),
       ggplot(elbow_data, aes(k, wss)) + geom_line() + geom_point() + labs(title = "Elbow Plot", x = "k", y = "WSS"),
       width = 6, height = 4, dpi = 300)

sil_scores <- map_dbl(2:8, function(kk) {
  km <- kmeans(scale(features %>% select(-iso3, -Fragility_Status, -country)), centers = kk, nstart = 25)
  ss <- cluster::silhouette(km$cluster, dist(scale(features %>% select(-iso3, -Fragility_Status, -country))))
  mean(ss[, 3], na.rm = TRUE)
})

sil_tbl <- tibble(k = 2:8, mean_silhouette = sil_scores)

ggsave(here::here("outputs","cluster_silhouette.png"),
       ggplot(sil_tbl, aes(k, mean_silhouette)) + geom_line() + geom_point() + labs(title = "Average Silhouette by k", x = "k", y = "Mean silhouette"),
       width = 6, height = 4, dpi = 300)

# Pick k (inspect elbow/silhouette; default 4)
k_opt <- 4
km <- kmeans(scale(features %>% select(-iso3, -Fragility_Status, -country)), centers = k_opt, nstart = 25)
features$Cluster <- factor(km$cluster)

# Label clusters (rule‑based)
tmp_clust <- features %>% group_by(Cluster) %>% summarise(
  MCV1_2024 = mean(MCV1_2024, na.rm = TRUE),
  Mean_Dropout = mean(Mean_Dropout, na.rm = TRUE),
  Pct_Fragile = mean(Pct_Fragile, na.rm = TRUE), .groups = "drop")

label_map <- tmp_clust %>% mutate(Label = case_when(
  MCV1_2024 >= 85 & Mean_Dropout <= 15 ~ "High coverage, low dropout",
  MCV1_2024 < 70 & Mean_Dropout > 30 & Pct_Fragile >= 50 ~ "Very low coverage, high dropout, high fragility",
  MCV1_2024 < 80 & Mean_Dropout > 20 ~ "Low–moderate coverage, elevated dropout",
  TRUE ~ "Moderate coverage, moderate dropout"
)) %>% select(Cluster, Label)

features <- features %>% left_join(label_map, by = "Cluster")

cluster_summary <- features %>%
  left_join(wuenic %>% filter(Year == 2024) %>% select(iso3, Income_Level, Unvaccinated_MCV1), by = "iso3") %>%
  group_by(Cluster, Label) %>%
  summarise(
    Countries = n(),
    Mean_MCV1_2024 = mean(MCV1_2024, na.rm = TRUE),
    Mean_Dropout = mean(Mean_Dropout, na.rm = TRUE),
    Mean_Change_MCV1 = mean(Change_MCV1, na.rm = TRUE),
    Total_Unvaccinated = sum(Unvaccinated_MCV1, na.rm = TRUE),
    Pct_Fragile_Avg = mean(Pct_Fragile, na.rm = TRUE),
    Income_Modal = names(sort(table(Income_Level), decreasing = TRUE))[1], .groups = "drop"
  ) %>% arrange(Cluster)

readr::write_csv(cluster_summary, here::here("outputs","cluster_summary_2024.csv"))

# Recommendations by cluster label
recs <- features %>%
  left_join(wuenic %>% filter(Year == 2024) %>% select(iso3, Unvaccinated_MCV1), by = "iso3") %>%
  group_by(Label) %>%
  summarise(
    Mean_MCV1 = mean(MCV1_2024, na.rm = TRUE),
    Mean_Dropout = mean(Mean_Dropout, na.rm = TRUE),
    Total_Unvaccinated = sum(Unvaccinated_MCV1, na.rm = TRUE), .groups = "drop"
  ) %>%
  mutate(Recommendation = case_when(
    str_detect(Label, "Very low coverage") ~ "High‑intensity SIAs; security‑adapted delivery; supply restoration; social mobilisation.",
    str_detect(Label, "Low–moderate coverage") ~ "Catch‑up + microplanning; address access gaps; demand generation; 2nd‑dose tracking.",
    str_detect(Label, "High coverage, low dropout") ~ "Sustain gains; monitor heterogeneity; address pockets; outbreak readiness.",
    TRUE ~ "Mixed strategies: subnational SIAs, DQAs, outreach to zero‑dose communities."
  ))
readr::write_csv(recs, here::here("outputs","policy_recommendations_by_cluster.csv"))

# Table 2 (cluster profiles + recs)

table2_clusters <- cluster_summary %>%
  select(Cluster, Label, Countries, `Mean MCV1 (%)` = Mean_MCV1_2024,
         `Mean Dropout (%)` = Mean_Dropout, Total_Unvaccinated, `Pct Fragile (%)` = Pct_Fragile_Avg) %>%
  left_join(recs %>% select(Label, Recommendation), by = "Label") %>%
  mutate(`Mean MCV1 (%)` = round(`Mean MCV1 (%)`, 1), `Mean Dropout (%)` = round(`Mean Dropout (%)`, 1),
         `Unvaccinated (M)` = round(Total_Unvaccinated/1e6, 1), `Pct Fragile (%)` = round(`Pct Fragile (%)`, 1)) %>%
  arrange(as.numeric(as.character(Cluster)))

readr::write_csv(table2_clusters, here::here("outputs", "Table2_cluster_profiles_recommendations_2024.csv"))

# Scatter and map by cluster label
p_cluster_scatter <- ggplot(features, aes(MCV1_2024, Mean_Dropout, colour = Label, shape = Pct_Fragile > 50)) +
  geom_point(size = 2.8, alpha = 0.7, position = position_jitter(width = 0.5, height = 0.5), na.rm = TRUE) +
  labs(title = "LMIC Clusters: Coverage vs Dropout (2024)", x = "MCV1 (%)", y = "Dropout (%)", shape = "High fragility (>50%)") +
  scale_color_brewer(palette = "Set2")

ggsave(here::here("outputs","cluster_scatter.png"), p_cluster_scatter, width = 8, height = 5, dpi = 300)

map_clusters <- left_join(world, features %>% rename(iso_a3 = iso3), by = "iso_a3")
p_map_clusters <- ggplot(map_clusters) +
  geom_sf(aes(fill = Label), colour = "white", size = 0.1) +
  scale_fill_brewer(palette = "Set2", na.value = "grey90") +
  labs(title = "LMIC Cluster Profiles (2024)", fill = "Cluster profile")

ggsave(here::here("outputs","map_clusters_2024.png"), p_map_clusters, width = 10, height = 6, dpi = 300)

# ------------------------------------------------------------
# 12) Top‑N diagnostics and gap tables (supporting material)
# ------------------------------------------------------------

top_dropout <- data_2024 %>% filter(LMIC_Flag == "LMIC", !is.na(Dropout_Rate)) %>%
  arrange(desc(Dropout_Rate)) %>% transmute(country, iso3, Coverage_MCV1, Coverage_MCV2, Dropout_Rate, Dropout_Outlier) %>% slice_head(n = 10)
readr::write_csv(top_dropout, here::here("outputs","top10_dropout_2024.csv"))

negative_dropout <- data_2024 %>% filter(Dropout_Outlier == "<0% (check data)") %>% select(country, iso3, Coverage_MCV1, Coverage_MCV2, Dropout_Rate)
readr::write_csv(negative_dropout, here::here("outputs","Supp_Table_S6_negative_dropout_2024.csv"))

gap_table <- data_2024 %>% filter(LMIC_Flag == "LMIC", Coverage_MCV1 < 95) %>%
  transmute(country, iso3, Coverage_MCV1, Gap_to_95 = 95 - Coverage_MCV1) %>% arrange(desc(Gap_to_95))
readr::write_csv(gap_table, here::here("outputs","gap_to_herd_immunity_mcv1.csv"))

top_unvax <- data_2024 %>% filter(LMIC_Flag == "LMIC") %>% arrange(desc(Unvaccinated_MCV1)) %>% select(country, iso3, Unvaccinated_MCV1) %>% slice_head(n = 10)
readr::write_csv(top_unvax, here::here("outputs","top10_unvaccinated_mcv1_2024.csv"))

# ------------------------------------------------------------
# 13) Combined figure (IA2030 bar + Alluvial)
# ------------------------------------------------------------

base_title_size <- 11
p_ia2030_bar2 <- p_ia2030_bar + labs(title = "A) LMICs meeting ≥95% for MCV1 and MCV2", x = NULL, y = "Countries (%)") +
  theme(plot.title = element_text(face = "plain", size = base_title_size), legend.position = "bottom")

p_alluvial2 <- p_alluvial + labs(title = "B) Transitions between MCV1 coverage groups, 2019–2024", x = NULL, y = "Countries") +
  theme(plot.title = element_text(face = "plain", size = base_title_size), legend.position = "bottom")

fig3 <- (p_ia2030_bar2 | p_alluvial2) + plot_layout(widths = c(1, 1.2)) +
  plot_annotation(title = "Figure 3. IA2030 status and MCV1 coverage transitions in LMICs",
                  theme = theme(plot.title = element_text(face = "plain", size = base_title_size + 1)))

ggsave(here::here("outputs","Fig3_bar_then_alluvial.pdf"), fig3, width = 12, height = 6.5, dpi = 300)

# ------------------------------------------------------------
# 14) Session info (for reproducibility)
# ------------------------------------------------------------

capture.output(sessionInfo(), file = here::here("outputs","_sessionInfo.txt"))

message("Analysis completed. Outputs written to ./outputs/")
