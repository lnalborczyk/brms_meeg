###########################################################
# Summarise the simulation results                        #
# ------------------------------------------------------- #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@cnrs.fr                    #
# Last updated on November 4, 2026                        #
###########################################################

library(tidyverse)
library(MetBrewer)
library(patchwork)
library(gt)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

results_dir <- Sys.getenv("RESULTS_DIR", unset = "errors_results")
combined_dir <- file.path(results_dir, "combined")
fig_dir <- file.path(results_dir, "figures")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Optionally restrict inputs to comma-separated HPC job IDs
job_ids_env <- trimws(Sys.getenv("RESULT_JOB_IDS", unset = "") )
job_ids_to_include <- if (nzchar(job_ids_env) ) {
    
    trimws(strsplit(job_ids_env, ",", fixed = TRUE)[[1]])
    
    } else {
        
        character(0)
        
}

# Sampling: 250 distinct simulations in each effect_scenario x SNR condition.
# Four conditions give 1,000 simulations. Selection is deterministic: ascending
# sim_index, then simulation_id. Every method, threshold, and cluster belonging
# to a selected simulation is retained (no row-wise truncation).
simulations_per_condition <- 250L
expected_conditions <- 4L
cpt_variant_to_keep <- "cpt_auto"
selected_bgam_method <- "brms_odds_abs_qi95"
selected_bgam_baseline <- "posterior_cluster_threshold_abs_qi95"
default_bgam_odds_threshold <- 10
default_alpha_threshold <- 0.05

# Evaluation window, in seconds, and temporal sampling resolution
evaluation_start <- 0
evaluation_end <- 0.5
sampling_rate <- 500
evaluation_grid <- seq(evaluation_start, evaluation_end, by = 1 / sampling_rate)

method_labels <- c(
    raw_p = "Raw p-value",
    pval_bh = "FDR BH95",
    pval_by = "FDR BY01",
    pval_holm = "Holm",
    cluster_mass = "Cluster mass",
    cluster_tfce = "TFCE",
    cpt = "Change point",
    brms_odds_abs_qi95 = "BGAM (our method)"
    )

method_colours <- setNames(
    met.brewer(name = "Johnson", n = length(method_labels) ),
    unname(method_labels)
    )

# -----------------------------------------------------------------------------
# Input and summary helpers
# -----------------------------------------------------------------------------

extract_job_id <- function (path) {
    
    stringr::str_match(basename(path), "_hpc_([0-9]+)_[0-9]+\\.rds$")[, 2]
    
}

filter_files_by_job_id <- function (files, job_ids) {
    
    if (length(job_ids) == 0L) return (files)
    files[extract_job_id(files) %in% job_ids]
    
}

read_rds_set <- function (files) {
    
    purrr::map_dfr(files, function (file) {
        readRDS(file) |>
            as_tibble() |>
            mutate(
                source_file = basename(file),
                source_job_id = extract_job_id(file),
                source_mtime = file.info(file)$mtime
                )
        }
        )
}

assert_columns <- function (data, required, name) {
    
    missing <- setdiff(required, names(data) )
    
    if (length(missing) > 0L) {
        
        stop (name, " is missing columns: ", paste(missing, collapse = ", "))
        
    }
    
    invisible(data)
    
}

valid_error_row <- function (error) {
    
    is.na(error) | trimws(as.character(error) ) == ""
    
}

safe_mean <- function (x) {
    
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
    
}

safe_median <- function (x) {
    
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else median(x)
    
}

safe_mae <- function (x) safe_mean(abs(x) )

safe_rmse <- function (x) {
    
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else sqrt(mean(x^2) )
    
}

safe_sd <- function (x) {
    
    x <- x[is.finite(x)]
    if (length(x) < 2L) NA_real_ else sd(x)
    
}

mark_best <- function (x, minimise = TRUE) {
    
    values <- x[is.finite(x)]
    if (length(values) == 0L) return (rep(FALSE, length(x) ) )
    target <- if (minimise) min(values) else max(values)
    is.finite(x) & x == target
    
}

# Pick one source file per simulation. For results, prefer files containing at
# least one non-error row; otherwise select the latest source by modification
# time. All cluster/method rows from the chosen file remain together.
keep_latest_source <- function (data, prefer_valid = FALSE) {
    
    if (!"simulation_id" %in% names(data) ) return (data)

    source_choice <- data |>
        filter(!is.na(simulation_id)) |>
        mutate(row_valid = if (prefer_valid) valid_error_row(error) else TRUE) |>
        group_by(simulation_id, source_file, source_job_id, source_mtime) |>
        summarise(has_valid_rows = any(row_valid), .groups = "drop") |>
        arrange(
            simulation_id,
            desc(has_valid_rows),
            desc(source_mtime),
            desc(source_job_id),
            source_file
            ) |>
        group_by(simulation_id) |>
        slice_head(n = 1L) |>
        ungroup() |>
        select(simulation_id, source_file)

    inner_join(data, source_choice, by = c("simulation_id", "source_file") )
    
}

save_plot <- function (plot, filename, width, height) {
    
    ggsave(
        filename = file.path(fig_dir, filename),
        plot = plot,
        width = width,
        height = height,
        dpi = 300,
        bg = "white"
        )
    
}

# -----------------------------------------------------------------------------
# Load inputs and select the first 250 simulations per condition
# -----------------------------------------------------------------------------

sim_files <- list.files(
    results_dir,
    pattern = "^errors_results_extended_hpc_.*\\.rds$",
    full.names = TRUE
    ) |>
    filter_files_by_job_id(job_ids_to_include)

true_cluster_files <- list.files(
    results_dir,
    pattern = "^true_clusters_extended_hpc_.*\\.rds$",
    full.names = TRUE
    ) |>
    filter_files_by_job_id(job_ids_to_include)

if (length(sim_files) == 0L || length(true_cluster_files) == 0L) {
    
    stop ("Simulation or true-cluster RDS files were not found in: ", results_dir)
    
}

sim_results_raw <- read_rds_set(sim_files)
if (!"error" %in% names(sim_results_raw) ) sim_results_raw$error <- NA_character_

assert_columns(
    sim_results_raw,
    c("simulation_id", "sim_index", "effect_scenario", "snr", "method",
      "cluster_id", "onset", "offset"),
    "Simulation results"
    )

sim_results_latest <- keep_latest_source(sim_results_raw, prefer_valid = TRUE)

# One record per simulation, independent of the number of detected clusters,
# method variants, or thresholds associated with that simulation.
simulation_index <- sim_results_latest |>
    filter(!is.na(simulation_id), !is.na(effect_scenario), !is.na(snr) ) |>
    distinct(simulation_id, sim_index, effect_scenario, snr)

# A simulation ID should describe one condition and one simulation index.
repeated_ids <- simulation_index |>
    count(simulation_id, name = "n") |>
    filter(n != 1L)

if (nrow(repeated_ids) > 0L) {
    
    stop ("Some simulation IDs map to multiple conditions or indices.")
    
}

condition_counts <- simulation_index |> count(effect_scenario, snr, name = "n_available")

if (nrow(condition_counts) != expected_conditions) {
    
    print(condition_counts)
    
    stop (
        "Expected ", expected_conditions, " effect_scenario x SNR conditions, ",
        "but found ", nrow(condition_counts), ". Check the input data or ",
        "adjust the selection settings."
        )
    
}

if (any(condition_counts$n_available < simulations_per_condition) ) {
    
    print(condition_counts)
    
    stop (
        "At least one condition has fewer than ",
        simulations_per_condition, " simulations."
        )
    
}

selected_simulations <- simulation_index |>
    group_by(effect_scenario, snr) |>
    arrange(sim_index, simulation_id, .by_group = TRUE) |>
    slice_head(n = simulations_per_condition) |>
    ungroup() |>
    arrange(effect_scenario, snr, sim_index, simulation_id)

stopifnot(nrow(selected_simulations) == expected_conditions * simulations_per_condition)

# Save the included IDs, making the 1,000-simulation subset reproducible.
write_csv(
    selected_simulations,
    file.path(combined_dir, "selected_simulation_ids.csv")
    )

message ("Selected ", nrow(selected_simulations), " simulations:")
print(selected_simulations |> count(effect_scenario, snr) )
sim_results_latest <- sim_results_latest |> semi_join(selected_simulations, by = "simulation_id")
rm(sim_results_raw)

true_clusters_raw <- read_rds_set(true_cluster_files)

assert_columns(
    true_clusters_raw,
    c("simulation_id", "true_onset", "true_offset"),
    "True-cluster results"
    )

true_clusters <- true_clusters_raw |>
    keep_latest_source() |>
    semi_join(selected_simulations, by = "simulation_id") |>
    filter(
        is.finite(true_onset),
        is.finite(true_offset),
        true_offset >= true_onset
        ) |>
    select(-source_file, -source_job_id, -source_mtime) |>
    distinct()

rm(true_clusters_raw)

missing_truth <- setdiff(selected_simulations$simulation_id, unique(true_clusters$simulation_id) )

if (length(missing_truth) > 0L) {
    
    stop (
        "Missing true intervals for ", length(missing_truth),
        " selected simulations. First IDs: ",
        paste(head(missing_truth, 10L), collapse = ", ")
        )
    
}

# -----------------------------------------------------------------------------
# Normalise methods and select one BGAM baseline definition
# -----------------------------------------------------------------------------

standardise_methods <- function (data) {
    
    # For BGAM, odds_threshold is the inferential cutoff. 'threshold' is only
    # a legacy fallback. effect_threshold is the simulation-specific SESOI.
    for (field in c("threshold", "odds_threshold", "effect_threshold") ) {
        
        if (!field %in% names(data)) data[[field]] <- NA_real_
        
    }
    
    if (!"baseline_threshold_name" %in% names(data) ) {
        
        data$baseline_threshold_name <- NA_character_
        
    }

    data |>
        mutate(
            method_raw = as.character(method),
            threshold = suppressWarnings(as.numeric(as.character(threshold) ) ),
            odds_threshold = suppressWarnings(as.numeric(as.character(odds_threshold) ) ),
            effect_threshold = suppressWarnings(as.numeric(as.character(effect_threshold) ) ),
            baseline_threshold_name = as.character(baseline_threshold_name)
            ) |>
        filter(
            !method_raw %in% c("cpt_auto", "cpt_fixed_one_cluster") |
                method_raw == cpt_variant_to_keep
            ) |>
        mutate(
            method_key = if_else(
                method_raw == cpt_variant_to_keep,
                "cpt",
                method_raw
                )
            ) |>
        filter(method_key %in% names(method_labels) ) |>
        filter(
            method_key != selected_bgam_method |
                baseline_threshold_name == selected_bgam_baseline
            ) |>
        mutate(
            baseline_threshold_name = if_else(
                method_key == selected_bgam_method,
                baseline_threshold_name,
                NA_character_
                ),
            odds_threshold = if_else(
                method_key == selected_bgam_method,
                coalesce(odds_threshold, threshold),
                NA_real_
                ),
            effect_threshold = if_else(
                method_key == selected_bgam_method,
                effect_threshold,
                NA_real_
                ),
            threshold_value = if_else(
                method_key == selected_bgam_method,
                odds_threshold,
                threshold
                ),
            threshold_label = case_when(
                method_key == selected_bgam_method ~ paste0(
                    "Absolute 95% baseline interval; odds >= ", odds_threshold
                    ),
                is.finite(threshold) ~ paste0("alpha = ", threshold),
                TRUE ~ "No fixed threshold"
                ),
            is_default_threshold = case_when(
                method_key == selected_bgam_method ~
                    is.finite(odds_threshold) &
                    near(odds_threshold, default_bgam_odds_threshold),
                is.finite(threshold) ~
                    near(threshold, default_alpha_threshold),
                TRUE ~ FALSE
                ),
            method = factor(
                method_key,
                levels = names(method_labels),
                labels = unname(method_labels)
                ),
            method_label = as.character(method),
            effect_scenario_label = factor(
                recode(
                    effect_scenario,
                    one_cluster = "One true cluster",
                    two_clusters = "Two true clusters",
                    .default = effect_scenario
                    ),
                levels = c("One true cluster", "Two true clusters")
                )
            )
    
}

available_bgam <- sim_results_latest |> filter(as.character(method) == selected_bgam_method, valid_error_row(error) )

if (nrow(available_bgam) == 0L ||
    !"baseline_threshold_name" %in% names(available_bgam) ||
    !selected_bgam_baseline %in% available_bgam$baseline_threshold_name) {
    
    stop ("Selected BGAM method/baseline is absent from usable simulation results.")
    
}

sim_results <- sim_results_latest |>
    standardise_methods() |>
    mutate(valid_interval = is.finite(onset) & is.finite(offset) & offset >= onset)

available_odds <- sim_results |>
    filter(method_key == selected_bgam_method, valid_error_row(error) ) |>
    distinct(odds_threshold) |>
    pull(odds_threshold)

if (!any(near(available_odds, default_bgam_odds_threshold), na.rm = TRUE) ) {
    
    stop (
        "No usable BGAM results at posterior odds >= ",
        default_bgam_odds_threshold, ". Available odds: ",
        paste(sort(available_odds), collapse = ", ")
        )
    
}

message (
    "BGAM: ", selected_bgam_method,
    "; baseline: ", selected_bgam_baseline,
    "; default odds >= ", default_bgam_odds_threshold
    )

# -----------------------------------------------------------------------------
# Temporal interval evaluation
# -----------------------------------------------------------------------------

intervals_to_mask <- function (onset, offset, grid) {
    
    mask <- rep(FALSE, length(grid))
    valid <- which(is.finite(onset) & is.finite(offset) & offset >= onset)
    
    for (i in valid) mask <- mask | (grid >= onset[i] & grid <= offset[i])
    
    return (mask)
    
}

prepare_detected_intervals <- function (data) {
    
    data |>
        filter(valid_interval, is.finite(onset), is.finite(offset), offset >= onset) |>
        distinct(cluster_id, onset, offset) |>
        arrange(onset, offset)
    
}

prepare_true_intervals <- function (data) {
    
    data |>
        filter(is.finite(true_onset), is.finite(true_offset), true_offset >= true_onset) |>
        arrange(true_onset, true_offset)
    
}

evaluate_global_interval_set <- function (detected, truth, grid) {
    
    detected <- prepare_detected_intervals(detected)
    truth <- prepare_true_intervals(truth)

    detected_mask <- intervals_to_mask(detected$onset, detected$offset, grid)
    true_mask <- intervals_to_mask(truth$true_onset, truth$true_offset, grid)
    intersection <- sum(detected_mask & true_mask)
    union <- sum(detected_mask | true_mask)

    tibble(
        temporal_iou = if (union > 0L) intersection / union else NA_real_,
        temporal_precision = if (any(detected_mask)) {
            intersection / sum(detected_mask)
        } else {
            0
        },
        temporal_recall = if (any(true_mask)) {
            intersection / sum(true_mask)
        } else {
            NA_real_
        }
    )
    
}

make_true_cluster_windows <- function (truth, start, end) {
    
    truth <- prepare_true_intervals(truth)
    if (nrow(truth) == 0L) return(tibble() )

    if (nrow(truth) == 1L) {
        
        window_start <- start
        window_end <- end
        
    } else {
        
        # Partition the evaluation epoch at midpoints between true clusters.
        midpoints <- (head(truth$true_offset, -1L) + tail(truth$true_onset, -1L)) / 2
        window_start <- c(start, midpoints)
        window_end <- c(midpoints, end)
        
    }

    truth |>
        mutate(
            true_cluster_id = row_number(),
            window_start = window_start,
            window_end = window_end
            )
    
}

evaluate_local_true_clusters <- function (detected, truth, grid, start, end) {
    
    detected <- prepare_detected_intervals(detected)
    windows <- make_true_cluster_windows(truth, start, end)
    
    if (nrow(windows) == 0L) return (tibble() )

    map_dfr(seq_len(nrow(windows) ), function (i) {
        current <- windows[i, ]
        local_grid <- grid[grid >= current$window_start & grid <= current$window_end]

        detected_local <- detected |>
            filter(onset <= current$window_end, offset >= current$window_start) |>
            mutate(
                local_onset = pmax(onset, current$window_start),
                local_offset = pmin(offset, current$window_end)
                ) |>
            filter(local_offset >= local_onset)

        detected_mask <- intervals_to_mask(
            detected_local$local_onset, detected_local$local_offset, local_grid
            )
        
        true_mask <- intervals_to_mask(
            current$true_onset, current$true_offset, local_grid
            )

        # As in the original script, a finite boundary error requires a local
        # detection, not necessarily overlap with the true cluster. Recovery
        # is separately quantified by the overlap indicator below.
        estimated_onset <- if (nrow(detected_local) > 0L) {
            
            min(detected_local$local_onset)
            
        } else NA_real_
        
        estimated_offset <- if (nrow(detected_local) > 0L) { max(detected_local$local_offset)
            
        } else NA_real_

        tibble(
            true_cluster_id = current$true_cluster_id,
            true_onset = current$true_onset,
            true_offset = current$true_offset,
            recovered_true_cluster = any(detected_mask & true_mask),
            estimated_onset = estimated_onset,
            estimated_offset = estimated_offset,
            onset_error = estimated_onset - current$true_onset,
            offset_error = estimated_offset - current$true_offset
            )
        })
    
}

summary_grouping <- c(
    "effect_scenario", "effect_scenario_label", "snr",
    "method_key", "method", "method_label",
    "baseline_threshold_name", "odds_threshold",
    "threshold_value", "threshold_label", "is_default_threshold"
    )

analysis_groups <- sim_results |>
    filter(!is.na(simulation_id), !is.na(method_key), valid_error_row(error) ) |>
    group_by(simulation_id, sim_index, across(all_of(summary_grouping) ) )

global_interval_metrics <- analysis_groups |>
    group_modify(function (.x, .y) {
        truth <- true_clusters |> filter(simulation_id == .y$simulation_id[[1]])
        evaluate_global_interval_set(.x, truth, evaluation_grid) |>
            mutate(
                metric_window_start = evaluation_start,
                metric_window_end = evaluation_end,
                effect_threshold = safe_median(.x$effect_threshold)
                )
        }) |>
    ungroup()

local_cluster_metrics <- analysis_groups |>
    group_modify(function (.x, .y) {
        truth <- true_clusters |> filter(simulation_id == .y$simulation_id[[1]])
        evaluate_local_true_clusters(.x, truth, evaluation_grid, evaluation_start, evaluation_end) |>
            mutate(
                metric_window_start = evaluation_start,
                metric_window_end = evaluation_end
                )
        }) |>
    ungroup()

# -----------------------------------------------------------------------------
# Aggregate temporal overlap and boundary-localisation performance
# -----------------------------------------------------------------------------

# Temporal IoU is averaged across simulations. Boundary bias and error
# dispersion are calculated across finite local detections, including those
# that do not overlap the true cluster. Recovery is reported separately.
# Negative/positive signed bias denotes early/late estimated boundaries.
global_summary <- global_interval_metrics |>
    group_by(across(all_of(summary_grouping) ) ) |>
    summarise(
        n_simulations = n_distinct(simulation_id),
        median_effect_threshold = safe_median(effect_threshold),
        iou = safe_mean(temporal_iou),
        median_iou = safe_median(temporal_iou),
        temporal_recall = safe_median(temporal_recall),
        .groups = "drop"
        )

boundary_summary <- local_cluster_metrics |>
    mutate(
        onset_error_ms = onset_error * 1000,
        offset_error_ms = offset_error * 1000,
        has_local_detection = is.finite(estimated_onset) & is.finite(estimated_offset)
        ) |>
    group_by(across(all_of(summary_grouping) ) ) |>
    summarise(
        n_true_clusters_evaluated = n(),
        n_recovered_true_clusters = sum(recovered_true_cluster, na.rm = TRUE),
        true_cluster_recovery_rate = mean(recovered_true_cluster, na.rm = TRUE),
        n_true_clusters_with_local_detection = sum(has_local_detection),
        local_detection_rate = mean(has_local_detection),
        n_finite_onset_errors = sum(is.finite(onset_error_ms)),
        n_finite_offset_errors = sum(is.finite(offset_error_ms)),
        onset_bias_ms = safe_mean(onset_error_ms),
        onset_mae_ms = safe_mae(onset_error_ms),
        onset_rmse_ms = safe_rmse(onset_error_ms),
        onset_sd_ms = safe_sd(onset_error_ms),
        offset_bias_ms = safe_mean(offset_error_ms),
        offset_mae_ms = safe_mae(offset_error_ms),
        offset_rmse_ms = safe_rmse(offset_error_ms),
        offset_sd_ms = safe_sd(offset_error_ms),
        .groups = "drop"
        )

performance_summary <- global_summary |>
    left_join(boundary_summary, by = summary_grouping) |>
    mutate(threshold = threshold_value) |>
    arrange(effect_scenario, snr, method, threshold_value)

performance_summary_default <- performance_summary |> filter(is_default_threshold)

# Include simulation counts in exported data so missing detections can be
# distinguished from computational failures or incomplete method coverage.
write_csv(
    performance_summary,
    file.path(combined_dir, "simulation_performance_summary_all_thresholds.csv")
    )

write_csv(
    performance_summary_default,
    file.path(combined_dir, "simulation_performance_summary_default_thresholds.csv")
    )

saveRDS(
    performance_summary,
    file.path(combined_dir, "simulation_performance_summary_all_thresholds.rds")
    )

saveRDS(
    performance_summary_default,
    file.path(combined_dir, "simulation_performance_summary_default_thresholds.rds")
    )


# A compact audit of default-configuration coverage for all four conditions.
method_coverage <- performance_summary_default |>
    transmute(
        effect_scenario, snr, method_key, n_simulations,
        n_true_clusters_evaluated, n_recovered_true_clusters,
        n_true_clusters_with_local_detection,
        n_finite_onset_errors, n_finite_offset_errors
        )

write_csv(method_coverage, file.path(combined_dir, "method_coverage.csv") )

if (any(method_coverage$n_simulations < simulations_per_condition) ) {
    
    warning (
        "At least one default method has fewer than ",
        simulations_per_condition,
        " valid simulations in a condition. See method_coverage.csv."
        )
    
}

# -----------------------------------------------------------------------------
# Manuscript tables
# -----------------------------------------------------------------------------

manuscript_performance_table <- performance_summary_default |>
    transmute(
        scenario = effect_scenario_label,
        snr, method = method_label, threshold = threshold_label,
        n_simulations, recovery_rate = true_cluster_recovery_rate,
        iou, median_iou,
        onset_bias = onset_bias_ms,
        onset_mae = onset_mae_ms,
        onset_rmse = onset_rmse_ms,
        onset_sd = onset_sd_ms,
        offset_bias = offset_bias_ms,
        offset_mae = offset_mae_ms,
        offset_rmse = offset_rmse_ms,
        offset_sd = offset_sd_ms
        ) |>
    arrange(scenario, snr, method)

manuscript_performance_table_long <- manuscript_performance_table |>
    select(
        scenario, snr, method, iou,
        onset_bias, onset_mae, onset_rmse, onset_sd,
        offset_bias, offset_mae, offset_rmse, offset_sd
        ) |>
    pivot_longer(
        cols = matches("^(onset|offset)_(bias|mae|rmse|sd)$"),
        names_to = c("onset_offset", ".value"),
        names_pattern = "(onset|offset)_(bias|mae|rmse|sd)"
        ) |>
    mutate(
        onset_offset = factor(
            onset_offset, levels = c("onset", "offset"),
            labels = c("Onset", "Offset")
            )
        ) |>
    group_by(scenario, snr, onset_offset) |>
    mutate(
        bold_bias = mark_best(abs(bias) ),
        bold_mae = mark_best(mae),
        bold_rmse = mark_best(rmse),
        bold_sd = mark_best(sd),
        bold_iou = mark_best(iou, minimise = FALSE)
        ) |>
    ungroup() |>
    arrange(scenario, snr, onset_offset, mae)

write_csv(
    manuscript_performance_table,
    file.path(combined_dir, "manuscript_performance_table_default_thresholds.csv")
    )

write_csv(
    manuscript_performance_table_long,
    file.path(combined_dir, "manuscript_performance_table_long_default_thresholds.csv")
    )

saveRDS(
    manuscript_performance_table,
    file.path(combined_dir, "manuscript_performance_table_default_thresholds.rds")
    )

saveRDS(
    manuscript_performance_table_long,
    file.path(combined_dir, "manuscript_performance_table_long_default_thresholds.rds")
    )

# Print this object in Quarto, or use manuscript_performance_latex in LaTeX.
manuscript_performance_gt <- manuscript_performance_table_long |>
    mutate(snr = as.character(snr), scenario = as.character(scenario) ) |>
    group_by(scenario) |>
    gt() |>
    fmt_number(columns = c(iou, bias, mae, rmse, sd), decimals = 2) |>
    cols_align(align = "center", columns = everything() ) |>
    cols_align(align = "left", columns = method) |>
    cols_label(
        snr = "SNR", onset_offset = "Boundary", method = "Method",
        iou = "IoU", bias = "Bias", mae = "MAE", rmse = "RMSE", sd = "SD"
        ) |>
    tab_spanner(label = "Boundary localisation", columns = c(bias, mae, rmse, sd) ) |>
    tab_spanner(label = "Temporal overlap", columns = iou) |>
    tab_options(table.align = "center", table.font.size = 12) |>
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = bias, rows = bold_bias)
        ) |>
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = mae, rows = bold_mae)
        ) |>
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = rmse, rows = bold_rmse)
        ) |>
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = sd, rows = bold_sd)
        ) |>
    tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = iou, rows = bold_iou)
        ) |>
    cols_hide(columns = starts_with("bold_") )

manuscript_performance_latex <- as_latex(manuscript_performance_gt)

# -----------------------------------------------------------------------------
# Figure 1: temporal IoU and signed boundary biases, by scenario
# -----------------------------------------------------------------------------

compact_default <- performance_summary_default |>
    select(effect_scenario_label, snr, method, median_iou, onset_bias_ms, offset_bias_ms) |>
    pivot_longer(
        cols = c(median_iou, onset_bias_ms, offset_bias_ms),
        names_to = "metric", values_to = "value"
        ) |>
    mutate(
        metric = recode(
            metric,
            median_iou = "Temporal IoU",
            onset_bias_ms = "Onset bias (ms)",
            offset_bias_ms = "Offset bias (ms)"
            ),
        metric = factor(
            metric,
            levels = c("Temporal IoU", "Onset bias (ms)", "Offset bias (ms)")
            )
        )

zero_bias_lines <- tibble(
    metric = factor(
        c("Onset bias (ms)", "Offset bias (ms)"),
        levels = levels(compact_default$metric)
        )
    )

p_compact <- compact_default |>
    ggplot(aes(x = snr, y = value, group = method, colour = method) ) +
    geom_hline(
        data = zero_bias_lines,
        aes(yintercept = 0),
        linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE
        ) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 2, na.rm = TRUE) +
    facet_grid(metric ~ effect_scenario_label, scales = "free_y") +
    scale_colour_manual(values = method_colours, drop = FALSE) +
    labs(x = "SNR", y = NULL, colour = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", legend.box = "horizontal") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE))

save_plot(p_compact, "compact_localisation_bias_summary.png", 12.5, 10)

p_compact_poster <- p_compact +
    theme_bw(base_size = 25) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom"
        )

save_plot(p_compact_poster, "compact_localisation_bias_summary_poster.png", 12.5, 10)

# -----------------------------------------------------------------------------
# Figure 2: one-cluster temporal overlap and boundary bias
# -----------------------------------------------------------------------------

one_cluster_summary <- performance_summary_default |>
    filter(effect_scenario == "one_cluster") |>
    mutate(snr = readr::parse_number(as.character(snr) ) )

if (nrow(one_cluster_summary) == 0L) {
    
    stop ("No one-cluster results available for the one-cluster figure.")
    
}

# Order methods by average absolute boundary bias (display ordering only).
method_bias_order <- one_cluster_summary |>
    group_by(method) |>
    summarise(
        average_abs_bias = safe_mean(c(abs(onset_bias_ms), abs(offset_bias_ms) ) ),
        .groups = "drop"
        ) |>
    arrange(average_abs_bias) |>
    pull(method) |>
    as.character()

one_cluster_summary <- one_cluster_summary |>
    mutate(
        method = factor(as.character(method), levels = method_bias_order),
        snr_label = factor(
            snr,
            levels = sort(unique(snr) ),
            labels = paste0("SNR = ", sort(unique(snr) ) )
            )
        )

plot_method_colours <- method_colours[method_bias_order]

p_iou <- ggplot(
    one_cluster_summary,
    aes(x = snr, y = median_iou, colour = method, group = method)
    ) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 2.2, na.rm = TRUE) +
    scale_colour_manual(values = plot_method_colours, drop = TRUE) +
    scale_x_continuous(breaks = sort(unique(one_cluster_summary$snr) ) ) +
    coord_cartesian(ylim = c(0, 1) ) +
    labs(tag = "A", x = "SNR", y = "Temporal IoU with ground truth", colour = NULL) +
    guides(colour = guide_legend(ncol = 1, order = 1) )

p_boundary <- ggplot(
    one_cluster_summary,
    aes(y = forcats::fct_rev(method), colour = method)
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_segment(
        aes(x = onset_bias_ms, xend = offset_bias_ms, yend = forcats::fct_rev(method) ),
        linewidth = 0.8, na.rm = TRUE
        ) +
    geom_point(aes(x = onset_bias_ms, shape = "Onset"), size = 2.6, na.rm = TRUE) +
    geom_point(aes(x = offset_bias_ms, shape = "Offset"), size = 2.6, na.rm = TRUE) +
    facet_wrap(~ snr_label, nrow = 1) +
    scale_x_continuous(breaks = c(-50, 0, 50)) +
    scale_colour_manual(values = plot_method_colours, drop = TRUE) +
    scale_shape_manual(values = c(Onset = 16, Offset = 17)) +
    labs(tag = "B", x = "Bias (ms)", y = NULL, colour = NULL, shape = NULL) +
    guides(colour = "none", shape = guide_legend(ncol = 1, order = 2) )

# The legend sits beside panel A; both plots occupy the same left column.
figure_design <- "
AL
BS
"

p_two_part <- wrap_plots(
    A = free(p_iou, type = "label", side = "l"),
    B = p_boundary,
    L = guide_area(),
    S = plot_spacer(),
    design = figure_design,
    widths = c(0.8, 0.2),
    heights = c(2, 2),
    guides = "collect"
    ) &
    theme(legend.box.margin = margin(t = 12, r = 4, b = 0, l = 8))

save_plot(p_two_part, "one_cluster_temporal_performance.png", 12.5, 10)

p_two_part_poster <- p_two_part &
    theme(
        plot.tag = element_text(face = "bold", size = 14)
        )

save_plot(p_two_part_poster, "one_cluster_temporal_performance_poster.png", 9, 6)

# Save intermediate results for Quarto
saveRDS(
    local_cluster_metrics,
    file.path(combined_dir, "local_cluster_metrics.rds")
    )

saveRDS(
    global_interval_metrics,
    file.path(combined_dir, "global_interval_metrics.rds")
    )

message ("Analysis complete. Output folders: ", combined_dir, " ; ", fig_dir)
