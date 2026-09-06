# ============================================================================
# Pine marten camera-trap protocol comparison
# Standalone final analysis: explicit 30-minute independent-event threshold
# ============================================================================

# Expected project layout:
# PineMarten_ESE/
#   dataraw2023/deployment-tableSCPM.csv
#   dataraw2023/observation-tableSCPM.csv
#   dataraw2024/deployments_Collection1_ks.csv
#   dataraw2024/deployments_Collection2_ks.csv
#   dataraw2024/deployments_Collection3_ks.csv
#   dataraw2024/observations.csv
#   scripts/pine_marten_analysis_final.R
#
# Run from the PineMarten_ESE project root with:
# source("scripts/pine_marten_analysis_final.R")

required_packages <- c(
  "tidyverse", "lubridate", "MASS", "survival", "survminer",
  "survRM2", "broom", "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running the analysis: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(MASS)
  library(survival)
  library(survminer)
  library(survRM2)
  library(broom)
  library(patchwork)
})

# ---- 1. Settings and paths -------------------------------------------------

independence_minutes <- 30
analysis_start <- ymd_hms("2023-07-01 00:00:00", tz = "UTC")
ground_analysis_end <- ymd_hms("2023-12-01 00:00:00", tz = "UTC")
tree_analysis_end <- ymd_hms("2024-12-01 00:00:00", tz = "UTC")

paths <- list(
  ground_deployments = "dataraw2023/deployment-tableSCPM.csv",
  ground_observations = "dataraw2023/observation-tableSCPM.csv",
  
  tree_r1 = "dataraw2024/deployments_Collection1_ks.csv",
  tree_r2 = "dataraw2024/deployments_Collection2_ks.csv",
  tree_r3 = "dataraw2024/deployments_Collection3_ks.csv",
  tree_observations = "dataraw2024/observations.csv"
)

missing_files <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing_files) > 0) {
  stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
}

out_dir <- "outputs"
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

protocol_levels <- c(
  "Ground-based lure",
  "Tree-mounted bait"
)

protocol_colours <- c(
  "Ground-based lure" = "#FDAE61",
  "Tree-mounted bait" = "#66BD63"
)

# ---- 2. Helpers ------------------------------------------------------------

parse_datetime <- function(x) {
  ymd_hms(as.character(x), quiet = TRUE, tz = "UTC")
}

assert_columns <- function(data, required, object_name) {
  absent <- setdiff(required, names(data))
  if (length(absent) > 0) {
    stop(object_name, " is missing columns: ", paste(absent, collapse = ", "))
  }
}

prepare_deployments <- function(data, protocol, analysis_end) {
  assert_columns(
    data,
    c("deploymentID", "locationID", "deploymentStart", "deploymentEnd"),
    protocol
  )
  
  data %>%
    transmute(
      deployment_id = str_to_lower(as.character(deploymentID)),
      site = str_to_upper(as.character(locationID)),
      deployment_start = parse_datetime(deploymentStart),
      deployment_end = parse_datetime(deploymentEnd),
      protocol = protocol
    ) %>%
    mutate(
      crop_start = pmax(deployment_start, analysis_start),
      crop_end = pmin(deployment_end, analysis_end),
      effort_days = as.numeric(difftime(crop_end, crop_start, units = "days"))
    ) %>%
    filter(!is.na(crop_start), !is.na(crop_end), effort_days > 0)
}

prepare_media_observations <- function(data, deployments) {
  assert_columns(
    data,
    c(
      "deploymentID", "mediaID", "eventStart", "eventEnd",
      "observationLevel", "scientificName"
    ),
    "observation table"
  )
  
  data %>%
    transmute(
      deployment_id = str_to_lower(as.character(deploymentID)),
      media_id = as.character(mediaID),
      observation_level = str_to_lower(as.character(observationLevel)),
      scientific_name = str_squish(as.character(scientificName)),
      record_start = parse_datetime(eventStart),
      record_end = parse_datetime(eventEnd)
    ) %>%
    filter(
      observation_level == "media",
      scientific_name == "Martes martes",
      !is.na(record_start)
    ) %>%
    mutate(record_end = coalesce(record_end, record_start)) %>%
    inner_join(
      deployments %>%
        dplyr::select(deployment_id, site, crop_start, crop_end, protocol),
      by = "deployment_id"
    ) %>%
    filter(record_start >= crop_start, record_start < crop_end) %>%
    # A media file can occasionally appear more than once in an export.
    group_by(protocol, site, deployment_id, media_id) %>%
    summarise(
      record_start = min(record_start),
      record_end = max(record_end),
      .groups = "drop"
    )
}

make_independent_events <- function(media_records, threshold_minutes = 30) {
  media_records %>%
    arrange(protocol, site, record_start, record_end) %>%
    group_by(protocol, site) %>%
    mutate(
      previous_record_end = lag(record_end),
      gap_minutes = as.numeric(
        difftime(record_start, previous_record_end, units = "mins")
      ),
      starts_new_event = is.na(gap_minutes) | gap_minutes >= threshold_minutes,
      independent_event_id = cumsum(starts_new_event)
    ) %>%
    group_by(protocol, site, independent_event_id) %>%
    summarise(
      event_start = min(record_start),
      event_end = max(record_end),
      media_records = n(),
      .groups = "drop"
    )
}

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p, format = "f", digits = 4))
}

# ---- 3. Import and prepare deployments ------------------------------------

ground_dep_raw <- read_csv(paths$ground_deployments, show_col_types = FALSE)
tree_r1_raw <- read_csv(paths$tree_r1, show_col_types = FALSE)
tree_r2_raw <- read_csv(paths$tree_r2, show_col_types = FALSE)
tree_r3_raw <- read_csv(paths$tree_r3, show_col_types = FALSE)

ground_deployments <- prepare_deployments(
  ground_dep_raw,
  protocol = "Ground-based lure",
  analysis_end = ground_analysis_end
)

tree_deployments <- bind_rows(tree_r1_raw, tree_r2_raw, tree_r3_raw) %>%
  prepare_deployments(
    protocol = "Tree-mounted bait",
    analysis_end = tree_analysis_end
  )

all_deployments <- bind_rows(ground_deployments, tree_deployments)

effort_by_site <- all_deployments %>%
  group_by(protocol, site) %>%
  summarise(
    start_date = min(crop_start),
    site_end = max(crop_end),
    trap_nights = sum(effort_days),
    .groups = "drop"
  ) %>%
  mutate(protocol = factor(protocol, levels = protocol_levels))

# ---- 4. Import observations and create 30-minute events -------------------

ground_obs_raw <- read_csv(paths$ground_observations, show_col_types = FALSE)
tree_obs_raw <- read_csv(paths$tree_observations, show_col_types = FALSE)

ground_media <- prepare_media_observations(ground_obs_raw, ground_deployments)
tree_media <- prepare_media_observations(tree_obs_raw, tree_deployments)

independent_events <- bind_rows(ground_media, tree_media) %>%
  make_independent_events(independence_minutes) %>%
  mutate(protocol = factor(protocol, levels = protocol_levels))

# ---- 5. Mandatory provenance checks ---------------------------------------

checkpoint <- effort_by_site %>%
  left_join(
    independent_events %>% count(protocol, site, name = "events"),
    by = c("protocol", "site")
  ) %>%
  mutate(events = replace_na(events, 0L)) %>%
  group_by(protocol) %>%
  summarise(
    sites = n(),
    detected_sites = sum(events > 0),
    revisit_sites = sum(events >= 2),
    independent_events = sum(events),
    exact_trap_nights = sum(trap_nights),
    reported_trap_nights = round(sum(trap_nights)),
    .groups = "drop"
  )

print(checkpoint)

expected_checkpoint <- tibble(
  protocol = factor(protocol_levels, levels = protocol_levels),
  sites = c(27L, 23L),
  detected_sites = c(17L, 17L),
  revisit_sites = c(7L, 15L),
  independent_events = c(40L, 197L),
  reported_trap_nights = c(2828, 2023)
)

checkpoint_comparison <- checkpoint %>%
  dplyr::select(-exact_trap_nights) %>%
  left_join(
    expected_checkpoint,
    by = "protocol",
    suffix = c("_observed", "_expected")
  )

check_columns <- c(
  "sites", "detected_sites", "revisit_sites",
  "independent_events", "reported_trap_nights"
)

for (column in check_columns) {
  observed <- checkpoint_comparison[[paste0(column, "_observed")]]
  expected <- checkpoint_comparison[[paste0(column, "_expected")]]
  if (!identical(as.numeric(observed), as.numeric(expected))) {
    stop(
      "Checkpoint failed for ", column, ". Observed: ",
      paste(observed, collapse = ", "), "; expected: ",
      paste(expected, collapse = ", ")
    )
  }
}

message("All provenance checkpoints passed.")

write_csv(checkpoint, file.path(out_dir, "analysis_checkpoints.csv"))
write_csv(independent_events, file.path(out_dir, "independent_events_30min.csv"))
write_csv(effort_by_site, file.path(out_dir, "effort_by_site.csv"))

# ---- 6. Naive occupancy ----------------------------------------------------

site_metrics <- effort_by_site %>%
  left_join(
    independent_events %>% count(protocol, site, name = "marten_events"),
    by = c("protocol", "site")
  ) %>%
  mutate(
    marten_events = replace_na(marten_events, 0L),
    detected = marten_events > 0,
    events_per_100_trap_nights = 100 * marten_events / trap_nights
  )

occupancy_table <- with(site_metrics, table(protocol, detected))
fisher_occupancy <- fisher.test(occupancy_table)

capture.output(
  occupancy_table,
  fisher_occupancy,
  file = file.path(out_dir, "fisher_naive_occupancy.txt")
)

# ---- 7. Event-rate analyses ------------------------------------------------

negative_binomial_model <- glm.nb(
  marten_events ~ protocol + offset(log(trap_nights)),
  data = site_metrics
)

poisson_model <- glm(
  marten_events ~ protocol + offset(log(trap_nights)),
  family = poisson(link = "log"),
  data = site_metrics
)

poisson_dispersion <- sum(
  residuals(poisson_model, type = "pearson")^2
) / df.residual(poisson_model)

rate_ratio_results <- tidy(
  negative_binomial_model,
  exponentiate = TRUE,
  conf.int = TRUE
)

write_csv(site_metrics, file.path(out_dir, "site_metrics.csv"))
write_csv(rate_ratio_results, file.path(out_dir, "negative_binomial_rate_ratios.csv"))

capture.output(
  summary(negative_binomial_model),
  rate_ratio_results,
  paste("Poisson dispersion:", poisson_dispersion),
  AIC(poisson_model, negative_binomial_model),
  file = file.path(out_dir, "event_rate_models.txt")
)

# ---- 8. Time to first detection -------------------------------------------

first_detection <- independent_events %>%
  group_by(protocol, site) %>%
  summarise(first_detection = min(event_start), .groups = "drop") %>%
  left_join(effort_by_site, by = c("protocol", "site")) %>%
  mutate(
    days_to_first_detection = as.numeric(
      difftime(first_detection, start_date, units = "days")
    )
  )

# Reproduces the manuscript's detected-stations-only comparison.
first_detection_wilcoxon <- wilcox.test(
  days_to_first_detection ~ protocol,
  data = first_detection,
  exact = FALSE
)

first_detection_summary <- first_detection %>%
  group_by(protocol) %>%
  summarise(
    detected_sites = n(),
    mean_days = mean(days_to_first_detection),
    sd_days = sd(days_to_first_detection),
    median_days = median(days_to_first_detection),
    se_days = sd_days / sqrt(detected_sites),
    ci_lower = mean_days - qt(0.975, df = detected_sites - 1) * se_days,
    ci_upper = mean_days + qt(0.975, df = detected_sites - 1) * se_days,
    .groups = "drop"
  )

# Sensitivity analysis retaining non-detected sites as right-censored.
first_detection_survival_data <- effort_by_site %>%
  left_join(
    first_detection %>% dplyr::select(protocol, site, first_detection),
    by = c("protocol", "site")
  ) %>%
  mutate(
    detection_status = as.integer(!is.na(first_detection)),
    detection_followup = if_else(
      detection_status == 1L,
      as.numeric(difftime(first_detection, start_date, units = "days")),
      as.numeric(difftime(site_end, start_date, units = "days"))
    )
  )

first_detection_logrank <- survdiff(
  Surv(detection_followup, detection_status) ~ protocol,
  data = first_detection_survival_data
)

write_csv(first_detection_summary, file.path(out_dir, "time_to_first_detection_summary.csv"))
write_csv(first_detection_survival_data, file.path(out_dir, "time_to_first_detection_survival_data.csv"))

capture.output(
  first_detection_wilcoxon,
  first_detection_logrank,
  file = file.path(out_dir, "time_to_first_detection_tests.txt")
)

# ---- 9. Revisit dynamics and survival -------------------------------------

revisit_data <- independent_events %>%
  arrange(protocol, site, event_start) %>%
  group_by(protocol, site) %>%
  summarise(
    first_event = first(event_start),
    second_event = if (n() >= 2) nth(event_start, 2) else as.POSIXct(NA, tz = "UTC"),
    n_events = n(),
    .groups = "drop"
  ) %>%
  left_join(
    effort_by_site %>% dplyr::select(protocol, site, site_end),
    by = c("protocol", "site")
  ) %>%
  mutate(
    revisit_status = as.integer(!is.na(second_event)),
    time_to_revisit_days = if_else(
      revisit_status == 1L,
      as.numeric(difftime(second_event, first_event, units = "days")),
      as.numeric(difftime(site_end, first_event, units = "days"))
    )
  )

revisit_survival <- Surv(
  time = revisit_data$time_to_revisit_days,
  event = revisit_data$revisit_status
)

revisit_km <- survfit(revisit_survival ~ protocol, data = revisit_data)
revisit_logrank <- survdiff(revisit_survival ~ protocol, data = revisit_data)
revisit_cox <- coxph(revisit_survival ~ protocol, data = revisit_data)
revisit_ph_test <- cox.zph(revisit_cox)

rmst_input <- revisit_data %>%
  mutate(arm = as.integer(protocol == "Tree-mounted bait"))

rmst_30 <- rmst2(
  time = rmst_input$time_to_revisit_days,
  status = rmst_input$revisit_status,
  arm = rmst_input$arm,
  tau = 30
)

rmst_50 <- rmst2(
  time = rmst_input$time_to_revisit_days,
  status = rmst_input$revisit_status,
  arm = rmst_input$arm,
  tau = 50
)

write_csv(revisit_data, file.path(out_dir, "time_to_first_revisit_by_site.csv"))

capture.output(
  table(revisit_data$protocol, revisit_data$revisit_status),
  revisit_logrank,
  summary(revisit_cox),
  revisit_ph_test,
  "RMST at 30 days:",
  rmst_30,
  "RMST at 50 days:",
  rmst_50,
  file = file.path(out_dir, "revisit_survival_results.txt")
)

revisit_dynamics <- independent_events %>%
  group_by(protocol, site) %>%
  mutate(
    days_since_first = as.numeric(
      difftime(event_start, min(event_start), units = "days")
    ),
    day_bin = floor(days_since_first)
  ) %>%
  ungroup() %>%
  filter(day_bin >= 0, day_bin <= 60) %>%
  count(protocol, day_bin, name = "events") %>%
  complete(protocol, day_bin = 0:60, fill = list(events = 0)) %>%
  left_join(
    first_detection %>% count(protocol, name = "detected_sites"),
    by = "protocol"
  ) %>%
  group_by(protocol) %>%
  arrange(day_bin, .by_group = TRUE) %>%
  mutate(
    mean_events_per_detected_site = events / detected_sites,
    cumulative_events_per_detected_site = cumsum(events) / detected_sites
  ) %>%
  ungroup()

write_csv(revisit_dynamics, file.path(out_dir, "revisit_dynamics.csv"))

# ---- 10. Manuscript summary table -----------------------------------------

manuscript_table <- site_metrics %>%
  group_by(protocol) %>%
  summarise(
    total_stations = n(),
    stations_with_pine_marten = sum(detected),
    naive_occupancy = mean(detected),
    total_trap_nights = round(sum(trap_nights)),
    independent_events = sum(marten_events),
    events_per_100_trap_nights = 100 * independent_events / total_trap_nights,
    .groups = "drop"
  )

write_csv(manuscript_table, file.path(out_dir, "Table1_survey_summary.csv"))

# ---- 11. Figures -----------------------------------------------------------

set.seed(20260426)

rate_model_p <- unname(
  summary(negative_binomial_model)$coefficients[2, "Pr(>|z|)"]
)

p_detection_rate <- ggplot(
  site_metrics,
  aes(x = protocol, y = events_per_100_trap_nights, colour = protocol, fill = protocol)
) +
  geom_boxplot(width = 0.5, alpha = 0.18, outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(width = 0.10, size = 1.5, alpha = 0.85, shape = 18) +
  annotate(
    "text",
    x = 1.5,
    y = max(site_metrics$events_per_100_trap_nights) * 0.95,
    label = format_p(rate_model_p),
    colour = "grey30"
  ) +
  scale_colour_manual(values = protocol_colours) +
  scale_fill_manual(values = protocol_colours) +
  labs(x = NULL, y = "Pine marten events per 100 trap nights") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", panel.grid = element_blank())

p_revisit_daily <- ggplot(
  revisit_dynamics,
  aes(x = day_bin, y = mean_events_per_detected_site, colour = protocol)
) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = protocol_colours) +
  labs(
    x = "Days since first detection",
    y = "Mean events per detected site",
    colour = NULL,
    title = "A"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_revisit_cumulative <- ggplot(
  revisit_dynamics,
  aes(x = day_bin, y = cumulative_events_per_detected_site, colour = protocol)
) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = protocol_colours) +
  labs(
    x = "Days since first detection",
    y = "Cumulative events per detected site",
    colour = NULL,
    title = "B"
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))

p_revisit_dynamics <- p_revisit_daily + p_revisit_cumulative

p_first_detection <- ggplot(
  first_detection,
  aes(x = protocol, y = days_to_first_detection, colour = protocol, fill = protocol)
) +
  geom_boxplot(width = 0.5, alpha = 0.18, outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(width = 0.10, size = 1.5, alpha = 0.85, shape = 18) +
  annotate(
    "text",
    x = 1.5,
    y = max(first_detection$days_to_first_detection) * 0.95,
    label = format_p(first_detection_wilcoxon$p.value),
    colour = "grey30"
  ) +
  scale_colour_manual(values = protocol_colours) +
  scale_fill_manual(values = protocol_colours) +
  labs(x = NULL, y = "Days to first pine marten detection") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", panel.grid = element_blank())

p_revisit_survival <- ggsurvplot(
  revisit_km,
  data = revisit_data,
  risk.table = TRUE,
  conf.int = FALSE,
  censor = TRUE,
  palette = unname(protocol_colours),
  legend.title = "",
  legend.labs = protocol_levels,
  xlab = "Days since first detection",
  ylab = "Probability of no revisit",
  break.time.by = 10,
  ggtheme = theme_bw(base_size = 14) + theme(panel.grid = element_blank())
)
p_revisit_survival$plot <- p_revisit_survival$plot +
  theme(legend.title = element_blank())
ggsave(
  file.path(fig_dir, "Figure_detection_rate.png"),
  p_detection_rate,
  width = 6.2, height = 4.8, dpi = 600, bg = "white"
)
ggsave(
  file.path(fig_dir, "Figure_detection_rate.pdf"),
  p_detection_rate,
  width = 6.2, height = 4.8, bg = "white"
)
ggsave(
  file.path(fig_dir, "Figure_revisit_dynamics.png"),
  p_revisit_dynamics,
  width = 7.1, height = 4.2, dpi = 600, bg = "white"
)
ggsave(
  file.path(fig_dir, "Figure_revisit_dynamics.pdf"),
  p_revisit_dynamics,
  width = 7.1, height = 4.2, bg = "white"
)
ggsave(
  file.path(fig_dir, "FigureS1_time_to_first_detection.png"),
  p_first_detection,
  width = 6.2, height = 4.8, dpi = 600, bg = "white"
)
ggsave(
  file.path(fig_dir, "Figure_time_to_first_revisit.png"),
  p_revisit_survival$plot,
  width = 6.8, height = 5.2, dpi = 600, bg = "white"
)
ggsave(
  file.path(fig_dir, "Figure_time_to_first_revisit.pdf"),
  p_revisit_survival$plot,
  width = 6.8, height = 5.2, bg = "white"
)

# ---- 12. Reproducibility record and final summary -------------------------

saveRDS(
  list(
    negative_binomial_model = negative_binomial_model,
    poisson_model = poisson_model,
    first_detection_logrank = first_detection_logrank,
    revisit_km = revisit_km,
    revisit_cox = revisit_cox,
    rmst_30 = rmst_30,
    rmst_50 = rmst_50
  ),
  file.path(out_dir, "fitted_models.rds")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("\n============================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")
print(checkpoint)
cat("\nFisher test for naive occupancy:\n")
print(fisher_occupancy)
cat("\nNegative-binomial rate ratios:\n")
print(rate_ratio_results)
cat("\nDetected-site time-to-first-detection test:\n")
print(first_detection_wilcoxon)
cat("\nCensored time-to-first-detection log-rank test:\n")
print(first_detection_logrank)
cat("\nRevisit log-rank test:\n")
print(revisit_logrank)
cat("\nCox proportional-hazards test:\n")
print(revisit_ph_test)
cat("\nOutputs saved in: ", normalizePath(out_dir), "\n", sep = "")

