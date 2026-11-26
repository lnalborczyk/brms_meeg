###########################################################
# Monte-Carlo simulation of onset/offset error properties #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@cnrs.fr                    #
# Last updated on October 13, 2026                        #
###########################################################

########################
# Handling errors
########################

options(
    error = function () {
        traceback(2)
        quit("no", status = 1, runLast = FALSE)
        }
    )

########################
# Packages
########################

library(changepoint)
library(reticulate)
library(assertthat)
library(doParallel)
library(tidybayes)
library(posterior)
library(foreach)
library(dplyr)
library(tidyr)
library(purrr)
library(glue)
library(brms)

###########################
# Source helper scripts
###########################

source ("code/utils.R")

simulate_eeg_source_file <- Sys.getenv(
    "SIMULATE_EEG_SOURCE",
    unset = "code/simulate_eeg.R"
    )

if (!file.exists(simulate_eeg_source_file) && file.exists("code/simulate_erp_data.R") ) {
    
    simulate_eeg_source_file <- "code/simulate_erp_data.R"
    
}

source (simulate_eeg_source_file)

meanpower <- unlist(read.table("code/meanpower.txt") )

########################
# Output directories
########################

dir.create("errors_results", showWarnings = FALSE, recursive = TRUE)

if (!nzchar(Sys.getenv("TMPDIR") ) ) {
    
    stop ("TMPDIR was not set by the SLURM wrapper.")
    
}

if (!dir.exists(Sys.getenv("TMPDIR") ) ) {
    
    stop ("TMPDIR does not exist: ", Sys.getenv("TMPDIR") )
    
}

message ("Using TMPDIR = ", Sys.getenv("TMPDIR") )

##########################
# HPC / local settings
##########################

args <- commandArgs(trailingOnly = TRUE)

array_id <- as.integer(
    Sys.getenv(
        "SLURM_ARRAY_TASK_ID",
        unset = ifelse(length(args) >= 1, args[1], "1")
        )
    )

n_array_jobs <- as.integer(Sys.getenv("N_ARRAY_JOBS", unset = "1") )

job_id <- Sys.getenv(
    "SLURM_ARRAY_JOB_ID",
    unset = Sys.getenv("SLURM_JOB_ID", unset = "local")
    )

task_id <- Sys.getenv(
    "SLURM_ARRAY_TASK_ID",
    unset = as.character(array_id)
    )

total_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "16") )

Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
    )

if (!nzchar(Sys.getenv("RETICULATE_PYTHON") ) ) {
    
    Sys.setenv(RETICULATE_PYTHON = Sys.which("python") )
    
}

# Set brms parameters
brms_chains <- as.integer(Sys.getenv("BRMS_CHAINS", unset = "4") )
brms_cores <- as.integer(Sys.getenv("BRMS_CORES", unset = "4") )
n_workers <- max(1, floor(total_cores / brms_cores) )

brms_warmup <- as.integer(Sys.getenv("BRMS_WARMUP", unset = "1000") )
brms_iter <- as.integer(Sys.getenv("BRMS_ITER", unset = "2000") )
brms_k <- as.integer(Sys.getenv("BRMS_K", unset = "20") )

# Set to NULL to use all posterior draws
n_epred_draws <- NULL

use_mne_cluster_methods <- TRUE

message ("array_id = ", array_id)
message ("n_array_jobs = ", n_array_jobs)
message ("job_id = ", job_id)
message ("task_id = ", task_id)
message ("total_cores = ", total_cores)
message ("brms_chains = ", brms_chains)
message ("brms_cores = ", brms_cores)
message ("initial n_workers = ", n_workers)
message ("total planned cores = ", n_workers * brms_cores)
message ("RETICULATE_PYTHON = ", Sys.getenv("RETICULATE_PYTHON"))

# Simulation parameters
n_trials <- 50
n_ppt <- 20
outvar <- 1
srate <- 500
epoch_duration <- 0.5

# Baseline used for both SESOI and posterior onset/offset thresholding
# It is also passed to simulate_eeg(), where it sets the scale used for SNR
baseline_window <- c(0, 0.1)

# Participant-level effect amplitude variability passed to simulate_eeg()
# These parameters scale the true cluster profile participant-by-participant
effect_amp_mean <- 1
effect_amp_sd <- 0.25
effect_amp_dist <- "truncated_normal"

# Participant-specific timing variability (seconds)
# Defaults correspond to 25 ms SD for onset and offset and can be overridden
# from SLURM with PARTICIPANT_ONSET_SD / PARTICIPANT_OFFSET_SD
participant_onset_sd <- as.numeric(
    Sys.getenv("PARTICIPANT_ONSET_SD", unset = "0.025")
    )

participant_offset_sd <- as.numeric(
    Sys.getenv("PARTICIPANT_OFFSET_SD", unset = "0.025")
    )

# Participant-level noise / overall variability passed to simulate_eeg()
noise_sd_mean <- 1
noise_sd_sd <- 0.5
noise_sd_dist <- "truncated_normal"

effect_scenarios <- list(
    one_cluster = tibble(
        true_cluster_id = 1L,
        true_onset = 0.160,
        true_offset = 0.320,
        relative_amplitude = 1
        )
    )

snr_values <- c(0.25, 0.50, 1.00, 2.00)

significance_thresholds <- c(0.05, 0.01, 0.005, 0.001)
threshold_values <- c(1, 3, 6, 10, 20, 50, 100)

cpt_auto_max_changepoints <- 10

n_sims_per_scenario <- as.integer(
    Sys.getenv("N_SIMS_PER_SCENARIO", unset = "250")
    )

n_sims <- n_sims_per_scenario

simulation_design_full <- tidyr::crossing(
    sim_index = seq_len(n_sims),
    effect_scenario = names(effect_scenarios),
    snr = snr_values
    ) |>
    mutate(
        simulation_id = sprintf(
            "s%05d_%s_snr%04.2f",
            sim_index,
            effect_scenario,
            snr
            ),
        simulation_id = gsub("\\.", "p", simulation_id),
        row_id = row_number(),
        sim_seed = 1000000L + row_id
        )

if (n_array_jobs <= n_sims) {
    # keep all SNR levels from the same
    # simulation replicate in the same array task.
    simulation_design <- simulation_design_full |>
        mutate(array_task = ((sim_index - 1L) %% n_array_jobs) + 1L)
    
} else {
    
    # Small smoke tests may use more array tasks than simulation replicates
    # In that case, distribute individual simulation cells across tasks
    simulation_design <- simulation_design_full |>
        mutate(array_task = ((row_id - 1L) %% n_array_jobs) + 1L)
    
}

simulation_design <- simulation_design |>
    filter(array_task == array_id) |>
    select(-row_id, -array_task)

if (nrow(simulation_design) == 0) {
    
    stop ("No simulation rows assigned to this array task.")
    
}

n_workers <- min(n_workers, nrow(simulation_design))

message ("n_sims per effect_scenario × snr = ", n_sims)
message ("Number of effect scenarios = ", length(effect_scenarios) )
message ("Number of SNR levels = ", length(snr_values) )
message ("Total simulation cells = ", nrow(simulation_design_full) )
message ("Simulation cells in this array task = ", nrow(simulation_design) )
message ("n_workers after truncation = ", n_workers)

simulation_design |>
    count(effect_scenario, snr) |>
    print(n = Inf)

# -------------------------------------------------------------------------
# Generic helper functions
# -------------------------------------------------------------------------

if (!exists("message_parallel") ) {
    
    message_parallel <- function (...) message (...)
    
}

get_scenario_name <- function (simulation_row) {
    
    as.character(simulation_row$effect_scenario[[1]])
    
}

safe_ttest <- function (x, mu = 0) {
    
    x <- x[is.finite(x)]

    if (length(x) < 2 || sd(x) == 0) {
        return (
            tibble(
                tval = NA_real_,
                pval = 1
                )
            )
        
    }

    out <- t.test(x = x, mu = mu)

    tibble(
        tval = unname(out$statistic)^2,
        pval = out$p.value
        )
    
}

safe_find_clusters <- function (
        data,
        threshold,
        method = "unknown",
        stop_on_error = as.logical(Sys.getenv("STOP_ON_CLUSTER_ERROR", unset = "FALSE") )
        ) {
    
    tryCatch (
        {
            if (!exists("find_clusters", mode = "function") ) {
                
                stop ("find_clusters() is not available.")
                
            }

            if (!all(c("time", "value") %in% names(data) ) ) {
                
                stop (
                    "Input data must contain columns `time` and `value`. Found: ",
                    paste(names(data), collapse = ", ")
                    )
                
            }

            data <- data |>
                mutate(
                    time = as.numeric(time),
                    value = as.numeric(value)
                    ) |>
                filter(is.finite(time), is.finite(value) ) |>
                arrange(time)

            if (nrow(data) == 0) {
                
                return (tibble() )
                
            }

            out <- find_clusters(data = data, threshold = threshold)

            if (is.null(out) || nrow(out) == 0) {
                
                tibble()
                
            } else {
                
                out
                
            }
        },
        error = function (e) {
            msg <- paste0(
                "find_clusters() failed for method = ",
                method,
                ", threshold = ",
                threshold,
                ": ",
                conditionMessage(e)
                )

            message_parallel(msg)

            if (stop_on_error) {
                
                stop (msg, call. = FALSE)
                
            }

            tibble()
            
        }
        
    )
    
}

empty_interval_result <- function (
        simulation_row,
        threshold,
        method,
        cluster_id = NA_integer_
        ) {
    
    scenario_name <- get_scenario_name(simulation_row)

    tibble(
        simulation_id = simulation_row$simulation_id,
        sim_index = simulation_row$sim_index,
        effect_scenario = scenario_name,
        snr = simulation_row$snr,
        true_n_clusters = nrow(effect_scenarios[[scenario_name]]),
        threshold = threshold,
        method = method,
        cluster_id = cluster_id,
        onset = NA_real_,
        offset = NA_real_,
        has_cluster = FALSE
        )
    
}

clusters_to_interval_result <- function (
        cluster_df,
        simulation_row,
        threshold,
        method
        ) {
    
    scenario_name <- get_scenario_name(simulation_row)

    if (is.null(cluster_df) || nrow(cluster_df) == 0) {
        
        return (
            empty_interval_result(
                simulation_row = simulation_row,
                threshold = threshold,
                method = method
                )
            )
        
    }

    cluster_df |>
        transmute(
            simulation_id = simulation_row$simulation_id,
            sim_index = simulation_row$sim_index,
            effect_scenario = scenario_name,
            snr = simulation_row$snr,
            true_n_clusters = nrow(effect_scenarios[[scenario_name]]),
            threshold = threshold,
            method = method,
            cluster_id = cluster_id,
            onset = cluster_onset,
            offset = cluster_offset,
            has_cluster = TRUE
            )
    
}

# -------------------------------------------------------------------------
# True-cluster metadata
# -------------------------------------------------------------------------

make_true_cluster_table <- function (simulation_row) {
    
    scenario_name <- get_scenario_name(simulation_row)

    effect_scenarios[[scenario_name]] |>
        transmute(
            simulation_id = simulation_row$simulation_id,
            sim_index = simulation_row$sim_index,
            effect_scenario = scenario_name,
            snr = simulation_row$snr,
            true_n_clusters = n(),
            true_cluster_id,
            true_onset,
            true_offset
            )
    
}

make_data_matrix <- function (ppt_data) {
    
    time_vec <- sort(unique(ppt_data$time) )
    time_chr <- sprintf("%.10f", time_vec)

    wide_df <- ppt_data |>
        mutate(time_chr = sprintf("%.10f", time) ) |>
        select(participant, time_chr, eeg_diff) |>
        pivot_wider(
            names_from = time_chr,
            values_from = eeg_diff
            ) |>
        arrange(participant)

    list(
        X = as.matrix(wide_df[, time_chr, drop = FALSE]),
        time_vec = time_vec
        )
    
}

# -------------------------------------------------------------------------
# Classical cluster helper functions
# -------------------------------------------------------------------------

run_cpt_fixed_one_cluster <- function (
        tests_results,
        simulation_row,
        alpha_level
        ) {
    
    tvals <- tests_results$tval
    tvals[!is.finite(tvals)] <- 0

    out <- tryCatch (
        {
            res <- cpt.meanvar(
                data = tvals,
                method = "BinSeg",
                Q = 2
                )

            cps <- sort(changepoint::cpts(res) )
            cps <- cps[is.finite(cps)]
            cps <- cps[cps > 1 & cps < nrow(tests_results)]

            if (length(cps) < 2) {
                
                return (tibble() )
                
            }

            tibble(
                cluster_id = 1L,
                cluster_onset = tests_results$time[cps[1]],
                cluster_offset = tests_results$time[cps[2]]
                )
            
        },
        error = function (e) {
            tibble()
            }
        )

    clusters_to_interval_result(
        cluster_df = out,
        simulation_row = simulation_row,
        threshold = alpha_level,
        method = "cpt_fixed_one_cluster"
        )
    
}

run_cpt_auto_clusters <- function (
        tests_results,
        simulation_row,
        alpha_level,
        n_participants,
        max_changepoints = 10
        ) {
    
    tvals <- tests_results$tval
    tvals[!is.finite(tvals)] <- 0

    out <- tryCatch (
        {
            res <- cpt.meanvar(
                data = tvals,
                method = "BinSeg",
                Q = max_changepoints,
                penalty = "MBIC"
                )

            cps <- sort(changepoint::cpts(res) )
            cps <- cps[is.finite(cps)]
            cps <- cps[cps > 1 & cps < nrow(tests_results)]

            segment_ends <- c(cps, nrow(tests_results) )
            segment_starts <- c(1, cps + 1)

            tval_threshold <- qt(
                p = 1 - alpha_level / 2,
                df = n_participants - 1
                )^2

            segments <- tibble(
                segment_id = seq_along(segment_starts),
                start_idx = segment_starts,
                end_idx = segment_ends
                ) |>
                mutate(
                    mean_tval = map2_dbl(
                        start_idx,
                        end_idx,
                        ~ mean(tvals[.x:.y], na.rm = TRUE)
                        ),
                    is_high = mean_tval > tval_threshold
                    ) |>
                filter(is_high)

            if (nrow(segments) == 0) {
                
                return (tibble() )
                
            }

            segments <- segments |>
                arrange(start_idx) |>
                mutate(
                    new_cluster = row_number() == 1 | start_idx > lag(end_idx) + 1,
                    cluster_id = cumsum(new_cluster)
                    ) |>
                summarise(
                    start_idx = min(start_idx),
                    end_idx = max(end_idx),
                    .by = cluster_id
                    )

            segments |>
                transmute(
                    cluster_id = as.integer(cluster_id),
                    cluster_onset = tests_results$time[start_idx],
                    cluster_offset = tests_results$time[end_idx]
                    )
        },
        error = function (e) {
            
            tibble()
            
        }
        
    )

    clusters_to_interval_result(
        cluster_df = out,
        simulation_row = simulation_row,
        threshold = alpha_level,
        method = "cpt_auto"
        )
    
}

# -------------------------------------------------------------------------
# BGAM posterior contrast and baseline-derived thresholds
# -------------------------------------------------------------------------

get_group_posterior_contrast <- function (meta_gam, n_epred_draws = NULL) {
    
    time_vec <- sort(unique(meta_gam$data$time) )
    condition_levels <- levels(meta_gam$data$condition)

    if (length(condition_levels) != 2L) {
        
        stop ("Expected exactly two condition levels.")
        
    }

    reference_condition <- condition_levels[1]
    target_condition <- condition_levels[2]
    
    # build a time x condition grid and add participant = NA only so that
    # brms can validate a model that contains participant-level terms.
    newdata_contrast <- tidyr::crossing(
        time = time_vec,
        condition = condition_levels
        ) |>
        mutate(
            condition = factor(
                condition,
                levels = condition_levels
                ),
            participant = NA,
            eeg_sem = 1
            )

    contrasts(newdata_contrast$condition) <- contrasts(meta_gam$data$condition)

    # Population-level expected predictions. re_formula = NA excludes
    # participant-specific varying effects and factor smooths.
    if (is.null(n_epred_draws) ) {
        
        epred_df <- tidybayes::add_epred_draws(
            newdata = newdata_contrast,
            object = meta_gam,
            re_formula = NA,
            incl_autocor = FALSE
            )
        
    } else {
        
        epred_df <- tidybayes::add_epred_draws(
            newdata = newdata_contrast,
            object = meta_gam,
            re_formula = NA,
            incl_autocor = FALSE,
            ndraws = n_epred_draws
            )
        
    }

    posterior_contrast <- epred_df |>
        data.frame() |>
        mutate(condition = as.character(condition) ) |>
        select(.draw, time, condition, .epred) |>
        pivot_wider(
            names_from = condition,
            values_from = .epred
            ) |>
        mutate(
            effect =
                .data[[target_condition]] -
                .data[[reference_condition]]
            ) |>
        select(.draw, time, effect)

    list(
        posterior_contrast = posterior_contrast,
        reference_condition = reference_condition,
        target_condition = target_condition,
        n_draws_eff = n_distinct(posterior_contrast$.draw)
        )
    
}

compute_baseline_thresholds <- function (
        posterior_contrast,
        baseline_window = c(0, 0.1)
        ) {
    
    baseline_post_preds <- posterior_contrast |>
        filter(
            time >= baseline_window[1],
            time < baseline_window[2]
            ) |>
        group_by(time) |>
        summarise(
            epred_median = median(effect, na.rm = TRUE),
            epred_q025 = quantile(effect, 0.025, na.rm = TRUE),
            epred_q975 = quantile(effect, 0.975, na.rm = TRUE),
            .groups = "drop"
            )

    if (nrow(baseline_post_preds) == 0) {
        
        stop ("No posterior contrast values found in baseline_window.")
        
    }

    meta_gam_sesoi_positive <- max(
        baseline_post_preds$epred_median,
        na.rm = TRUE
        )

    meta_gam_sesoi_absolute <- max(
        abs(baseline_post_preds$epred_median),
        na.rm = TRUE
        )

    # Threshold for draw-wise posterior onset/offset extraction:
    # maximum absolute endpoint of the pointwise 95% posterior interval during baseline.
    posterior_cluster_threshold_abs_qi95 <- max(
        abs(c(baseline_post_preds$epred_q025, baseline_post_preds$epred_q975) ),
        na.rm = TRUE
        )

    list(
        baseline_post_preds = baseline_post_preds,
        meta_gam_sesoi_positive = meta_gam_sesoi_positive,
        meta_gam_sesoi_absolute = meta_gam_sesoi_absolute,
        posterior_cluster_threshold_abs_qi95 = posterior_cluster_threshold_abs_qi95
        )
    
}

make_baseline_threshold_table <- function (
        baseline_thresholds,
        simulation_row
        ) {
    
    scenario_name <- get_scenario_name(simulation_row)

    tibble(
        simulation_id = simulation_row$simulation_id,
        sim_index = simulation_row$sim_index,
        effect_scenario = scenario_name,
        snr = simulation_row$snr,
        baseline_threshold_name = c(
            "meta_gam_sesoi_positive",
            "meta_gam_sesoi_absolute",
            "posterior_cluster_threshold_abs_qi95"
            ),
        effect_threshold = c(
            baseline_thresholds$meta_gam_sesoi_positive,
            baseline_thresholds$meta_gam_sesoi_absolute,
            baseline_thresholds$posterior_cluster_threshold_abs_qi95
            )
        )
    
}

compute_posterior_odds <- function (
        posterior_contrast,
        effect_threshold,
        n_draws_eff
        ) {
    
    posterior_contrast |>
        group_by(time) |>
        summarise(
            posterior_probability = mean(
                effect > effect_threshold,
                na.rm = TRUE
                ),
            .groups = "drop"
            ) |>
        mutate(
            # Avoid infinite odds when all posterior draws fall on one side
            # of the effect threshold.
            posterior_probability = pmin(
                pmax(posterior_probability, 1 / n_draws_eff),
                1 - 1 / n_draws_eff
                ),
            posterior_odds = posterior_probability / (1 - posterior_probability)
            )
    
}

get_bgam_method_name <- function (baseline_threshold_name) {
    
    switch (
        baseline_threshold_name,
        meta_gam_sesoi_positive = "brms_odds_sesoi_positive",
        meta_gam_sesoi_absolute = "brms_odds_sesoi_absolute",
        posterior_cluster_threshold_abs_qi95 = "brms_odds_abs_qi95",
        paste0("brms_odds_", baseline_threshold_name)
        )
    
}

run_bgam_odds_grid <- function (
        posterior_contrast,
        baseline_threshold_values,
        simulation_row,
        threshold_values,
        n_draws_eff
        ) {

    pmap_dfr(
        baseline_threshold_values |>
            select(
                baseline_threshold_name,
                effect_threshold
                ),
        function (baseline_threshold_name, effect_threshold) {

            method_name <- get_bgam_method_name(baseline_threshold_name)

            posterior_odds_df <- compute_posterior_odds(
                posterior_contrast = posterior_contrast,
                effect_threshold = effect_threshold,
                n_draws_eff = n_draws_eff
                )

            map_dfr(
                threshold_values,
                function (odds_threshold) {

                    safe_find_clusters(
                        data = posterior_odds_df |>
                            select(
                                time,
                                value = posterior_odds
                                ),
                        threshold = odds_threshold,
                        method = method_name
                        ) |>
                        clusters_to_interval_result(
                            simulation_row = simulation_row,
                            threshold = odds_threshold,
                            method = method_name
                            ) |>
                        mutate(
                            baseline_threshold_name = baseline_threshold_name,
                            effect_threshold = effect_threshold,
                            odds_threshold = odds_threshold
                            )
                    
                }
                
            )
            
        }
        
    )
    
}

# -------------------------------------------------------------------------
# Matching detected clusters to true clusters
# -------------------------------------------------------------------------

match_detected_to_true <- function (sim_results, true_clusters) {
    
    detected <- sim_results |>
        filter(has_cluster) |>
        select(
            simulation_id,
            sim_index,
            effect_scenario,
            snr,
            true_n_clusters,
            threshold,
            method,
            cluster_id,
            onset,
            offset
            )

    if (nrow(detected) == 0) {
        
        return (
            tibble(
                simulation_id = character(),
                sim_index = integer(),
                effect_scenario = character(),
                snr = numeric(),
                true_n_clusters = integer(),
                threshold = numeric(),
                method = character(),
                cluster_id = integer(),
                onset = numeric(),
                offset = numeric(),
                true_cluster_id = integer(),
                true_onset = numeric(),
                true_offset = numeric(),
                intersection = numeric(),
                union = numeric(),
                temporal_iou = numeric(),
                onset_error = numeric(),
                offset_error = numeric(),
                onset_abs_error = numeric(),
                offset_abs_error = numeric(),
                matched_true_cluster = logical()
                )
            )
        
    }

    detected |>
        left_join(
            true_clusters,
            by = c(
                "simulation_id",
                "sim_index",
                "effect_scenario",
                "snr",
                "true_n_clusters"
                ),
            relationship = "many-to-many"
            ) |>
        mutate(
            intersection = pmax(
                0,
                pmin(offset, true_offset) - pmax(onset, true_onset)
                ),
            union = pmax(offset, true_offset) - pmin(onset, true_onset),
            temporal_iou = if_else(union > 0, intersection / union, 0),
            onset_error = onset - true_onset,
            offset_error = offset - true_offset,
            onset_abs_error = abs(onset_error),
            offset_abs_error = abs(offset_error)
            ) |>
        group_by(simulation_id, threshold, method, cluster_id) |>
        arrange(
            desc(temporal_iou),
            onset_abs_error + offset_abs_error,
            .by_group = TRUE
            ) |>
        slice(1) |>
        ungroup() |>
        mutate(matched_true_cluster = temporal_iou > 0)
    
}

# -------------------------------------------------------------------------
# Run one simulation cell
# -------------------------------------------------------------------------

run_one_simulation_cell <- function (row_i) {
    
    simulation_row <- simulation_design[row_i, ]

    scenario_name <- get_scenario_name(simulation_row)
    true_clusters <- effect_scenarios[[scenario_name]]
    sim_seed <- as.integer(simulation_row$sim_seed[[1]])

    set.seed(sim_seed)

    message_parallel(
        sprintf(
            "\nRunning %s (%d/%d)...",
            simulation_row$simulation_id,
            row_i,
            nrow(simulation_design)
            )
        )

    sim_obj <- simulate_eeg(
        seed = sim_seed,
        n_trials = n_trials,
        n_ppt = n_ppt,
        outvar = outvar,
        srate = srate,
        epoch_duration = epoch_duration,
        meanpower = meanpower,
        true_clusters = true_clusters,
        snr = simulation_row$snr,
        effect_amp_mean = effect_amp_mean,
        effect_amp_sd = effect_amp_sd,
        effect_amp_dist = effect_amp_dist,
        participant_onset_sd = participant_onset_sd,
        participant_offset_sd = participant_offset_sd,
        noise_sd_mean = noise_sd_mean,
        noise_sd_sd = noise_sd_sd,
        noise_sd_dist = noise_sd_dist,
        baseline_window = baseline_window,
        reference_condition = "cond1",
        target_condition = "cond2",
        return_list = TRUE
        )

    raw_df <- sim_obj$raw_df

    participant_true_clusters <- sim_obj$participant_true_clusters |>
        mutate(
            simulation_id = simulation_row$simulation_id,
            sim_index = simulation_row$sim_index,
            effect_scenario = scenario_name,
            snr = simulation_row$snr,
            .before = 1
            )

    ppt_data <- raw_df |>
        summarise(
            eeg = mean(eeg),
            .by = c(participant, condition, time)
            ) |>
        pivot_wider(
            names_from = condition,
            values_from = eeg
            ) |>
        mutate(eeg_diff = cond2 - cond1)

    summary_df <- raw_df |>
        summarise(
            eeg_mean = mean(eeg),
            eeg_sd = sd(eeg),
            eeg_sem = eeg_sd / sqrt(n()),
            .by = c(participant, condition, time)
            ) |>
        mutate(
            participant = factor(participant),
            condition = factor(condition, levels = c("cond1", "cond2") )
            )

    contrasts(summary_df$condition) <- c(-0.5, 0.5)

    meta_gam <- brm(
        formula = eeg_mean | se(x = eeg_sem, sigma = TRUE) ~
            1 + condition +
            s(time, bs = "tp", k = brms_k, by = condition) +
            s(
                time,
                participant,
                bs = "fs",
                m = 1,
                k = 10,
                by = condition
                ),
        data = summary_df,
        family = gaussian(),
        warmup = brms_warmup,
        iter = brms_iter,
        chains = brms_chains,
        cores = brms_cores,
        refresh = 0,
        seed = sim_seed
        )

    contrast_obj <- get_group_posterior_contrast(
        meta_gam = meta_gam,
        n_epred_draws = n_epred_draws
        )

    posterior_contrast <- contrast_obj$posterior_contrast
    n_draws_eff <- contrast_obj$n_draws_eff

    baseline_thresholds <- compute_baseline_thresholds(
        posterior_contrast = posterior_contrast,
        baseline_window = baseline_window
        )

    baseline_threshold_values <- make_baseline_threshold_table(
        baseline_thresholds = baseline_thresholds,
        simulation_row = simulation_row
        )

    # For each of the three baseline-derived effect thresholds, compute
    # posterior odds P(effect > threshold) / P(effect <= threshold), then
    # identify clusters using the same threshold_values grid.
    #
    # This yields 3 baseline thresholds x 7 posterior-odds thresholds
    # = 21 BGAM cluster definitions per simulation.
    temp_bgam_results <- run_bgam_odds_grid(
        posterior_contrast = posterior_contrast,
        baseline_threshold_values = baseline_threshold_values,
        simulation_row = simulation_row,
        threshold_values = threshold_values,
        n_draws_eff = n_draws_eff
        )

    mat_obj <- make_data_matrix(ppt_data)
    data_matrix <- mat_obj$X
    time_vec <- mat_obj$time_vec
    n_participants_current <- length(unique(ppt_data$participant) )

    tests_results <- ppt_data |>
        group_by(time) |>
        summarise(
            safe_ttest(eeg_diff, mu = 0),
            .groups = "drop"
            ) |>
        mutate(
            pval_bh = p.adjust(p = pval, method = "BH"),
            pval_by = p.adjust(p = pval, method = "BY"),
            pval_holm = p.adjust(p = pval, method = "holm")
            )

    message_parallel(
        sprintf(
            "Diagnostics %s: min raw p = %.4g, min BH = %.4g, min BY = %.4g, min Holm = %.4g, max tval = %.3f",
            simulation_row$simulation_id,
            min(tests_results$pval, na.rm = TRUE),
            min(tests_results$pval_bh, na.rm = TRUE),
            min(tests_results$pval_by, na.rm = TRUE),
            min(tests_results$pval_holm, na.rm = TRUE),
            max(tests_results$tval, na.rm = TRUE)
            )
        )

    temp_other_results <- map_dfr(
        significance_thresholds,
        function (alpha_level) {
            onset_offset_raw_p <- safe_find_clusters(
                data = tests_results |>
                    mutate(value = -pval) |>
                    select(time, value),
                threshold = -alpha_level,
                method = "raw_p"
                ) |>
                clusters_to_interval_result(
                    simulation_row = simulation_row,
                    threshold = alpha_level,
                    method = "raw_p"
                    )

            onset_offset_pval_bh <- safe_find_clusters(
                data = tests_results |>
                    mutate(value = -pval_bh) |>
                    select(time, value),
                threshold = -alpha_level,
                method = "pval_bh"
                ) |>
                clusters_to_interval_result(
                    simulation_row = simulation_row,
                    threshold = alpha_level,
                    method = "pval_bh"
                    )

            onset_offset_pval_by <- safe_find_clusters(
                data = tests_results |>
                    mutate(value = -pval_by) |>
                    select(time, value),
                threshold = -alpha_level,
                method = "pval_by"
                ) |>
                clusters_to_interval_result(
                    simulation_row = simulation_row,
                    threshold = alpha_level,
                    method = "pval_by"
                    )

            onset_offset_pval_holm <- safe_find_clusters(
                data = tests_results |>
                    mutate(value = -pval_holm) |>
                    select(time, value),
                threshold = -alpha_level,
                method = "pval_holm"
                ) |>
                clusters_to_interval_result(
                    simulation_row = simulation_row,
                    threshold = alpha_level,
                    method = "pval_holm"
                    )

            onset_offset_cpt_fixed <- run_cpt_fixed_one_cluster(
                tests_results = tests_results,
                simulation_row = simulation_row,
                alpha_level = alpha_level
                )

            onset_offset_cpt_auto <- run_cpt_auto_clusters(
                tests_results = tests_results,
                simulation_row = simulation_row,
                alpha_level = alpha_level,
                n_participants = n_participants_current,
                max_changepoints = cpt_auto_max_changepoints
                )

            if (use_mne_cluster_methods) {
                
                p_values_cluster_mass <- freq_stats_cluster_matrix(
                    X = data_matrix,
                    timesteps = time_vec,
                    cluster_type = "mass",
                    alpha_level = alpha_level
                    )

                p_values_tfce <- freq_stats_cluster_matrix(
                    X = data_matrix,
                    timesteps = time_vec,
                    cluster_type = "tfce",
                    alpha_level = alpha_level
                    )

                onset_offset_cluster_mass <- safe_find_clusters(
                    data = p_values_cluster_mass |>
                        mutate(value = -pval) |>
                        select(time, value),
                    threshold = -alpha_level,
                    method = "cluster_mass"
                    ) |>
                    clusters_to_interval_result(
                        simulation_row = simulation_row,
                        threshold = alpha_level,
                        method = "cluster_mass"
                        )

                onset_offset_cluster_tfce <- safe_find_clusters(
                    data = p_values_tfce |>
                        mutate(value = -pval) |>
                        select(time, value),
                    threshold = -alpha_level,
                    method = "cluster_tfce"
                    ) |>
                    clusters_to_interval_result(
                        simulation_row = simulation_row,
                        threshold = alpha_level,
                        method = "cluster_tfce"
                        )
                
            } else {
                
                onset_offset_cluster_mass <- tibble()
                onset_offset_cluster_tfce <- tibble()
                
            }

            bind_rows(
                onset_offset_raw_p,
                onset_offset_pval_bh,
                onset_offset_pval_by,
                onset_offset_pval_holm,
                onset_offset_cpt_fixed,
                onset_offset_cpt_auto,
                onset_offset_cluster_mass,
                onset_offset_cluster_tfce
                )
            
        }
        
    )

    list(
        sim_results = bind_rows(
            temp_bgam_results,
            temp_other_results
            ),
        baseline_threshold_values = baseline_threshold_values,
        participant_true_clusters = participant_true_clusters
        )
    
}

# -------------------------------------------------------------------------
# Parallel execution
# -------------------------------------------------------------------------

reticulate_python <- Sys.getenv("RETICULATE_PYTHON")

cl <- makeCluster(n_workers)
registerDoParallel(cl)

clusterExport(
    cl,
    varlist = c(
        "reticulate_python",
        "use_mne_cluster_methods",
        "simulate_eeg_source_file"
        ),
    envir = environment()
    )

clusterEvalQ(
    cl,
    {
        if (nzchar(reticulate_python) ) {
            
            Sys.setenv(RETICULATE_PYTHON = reticulate_python)
            
        }

        Sys.setenv(
            OMP_NUM_THREADS = "1",
            MKL_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1",
            VECLIB_MAXIMUM_THREADS = "1",
            NUMEXPR_NUM_THREADS = "1"
            )

        library(changepoint)
        library(reticulate)
        library(assertthat)
        library(tidybayes)
        library(dplyr)
        library(tidyr)
        library(glue)
        library(brms)
        library(purrr)
        library(posterior)

        source("code/utils.R")
        source(simulate_eeg_source_file)

        message("Worker find_clusters loaded from: ", environmentName(environment(find_clusters) ) )

        if (!exists("simulate_eeg") ) {
            
            stop (
                "Could not find simulate_eeg() on worker. ",
                "Check SIMULATE_EEG_SOURCE or save the updated simulator as code/simulate_eeg.R."
                )
            
        }

        if (use_mne_cluster_methods) {
            
            reticulate::use_python(reticulate_python, required = TRUE)
            np <<- reticulate::import("numpy", delay_load = FALSE)
            mne <<- reticulate::import("mne", delay_load = FALSE)
            mne$set_log_level("WARNING")
            
        }

        NULL
        
    }
    
)

clusterExport(
    cl,
    varlist = c(
        "array_id",
        "n_array_jobs",
        "simulation_design",
        "effect_scenarios",
        "meanpower",
        "simulate_eeg_source_file",
        "n_trials",
        "n_ppt",
        "outvar",
        "srate",
        "epoch_duration",
        "baseline_window",
        "effect_amp_mean",
        "effect_amp_sd",
        "effect_amp_dist",
        "participant_onset_sd",
        "participant_offset_sd",
        "noise_sd_mean",
        "noise_sd_sd",
        "noise_sd_dist",
        "significance_thresholds",
        "threshold_values",
        "cpt_auto_max_changepoints",
        "brms_warmup",
        "brms_iter",
        "brms_chains",
        "brms_cores",
        "brms_k",
        "n_epred_draws",
        "use_mne_cluster_methods",
        "message_parallel",
        "get_scenario_name",
        "safe_ttest",
        "safe_find_clusters",
        "find_clusters",
        "freq_stats_cluster_matrix",
        "empty_interval_result",
        "clusters_to_interval_result",
        "make_true_cluster_table",
        "make_data_matrix",
        "run_cpt_fixed_one_cluster",
        "run_cpt_auto_clusters",
        "get_group_posterior_contrast",
        "compute_baseline_thresholds",
        "make_baseline_threshold_table",
        "compute_posterior_odds",
        "get_bgam_method_name",
        "run_bgam_odds_grid",
        "match_detected_to_true",
        "run_one_simulation_cell"
        ),
    envir = environment()
    )

combine_sim_outputs <- function (x, y) {
    
    list(
        sim_results = bind_rows(x$sim_results, y$sim_results),
        baseline_threshold_values = bind_rows(
            x$baseline_threshold_values,
            y$baseline_threshold_values
            ),
        participant_true_clusters = bind_rows(
            x$participant_true_clusters,
            y$participant_true_clusters
            )
        )
    
}

sim_output <- foreach(
    row_i = seq_len(nrow(simulation_design) ),
    .combine = combine_sim_outputs,
    .init = list(
        sim_results = tibble(),
        baseline_threshold_values = tibble(),
        participant_true_clusters = tibble()
        ),
    .errorhandling = "stop",
    .packages = c(
        "reticulate",
        "assertthat",
        "brms",
        "dplyr",
        "tidyr",
        "tidybayes",
        "changepoint",
        "purrr",
        "posterior"
        )
) %dopar% {
    tryCatch (
        {
            run_one_simulation_cell(row_i)
        },
        error = function (e) {
            simulation_row <- simulation_design[row_i, ]
            scenario_name <- as.character(simulation_row$effect_scenario[[1]])

            list(
                sim_results = tibble(
                    simulation_id = simulation_row$simulation_id,
                    sim_index = simulation_row$sim_index,
                    effect_scenario = scenario_name,
                    snr = simulation_row$snr,
                    true_n_clusters = nrow(effect_scenarios[[scenario_name]]),
                    threshold = NA_real_,
                    method = NA_character_,
                    cluster_id = NA_integer_,
                    onset = NA_real_,
                    offset = NA_real_,
                    has_cluster = FALSE,
                    error = e$message
                    ),
                baseline_threshold_values = tibble(),
                participant_true_clusters = tibble()
                )
        }
    )
}

stopCluster(cl)

sim_results <- sim_output$sim_results
baseline_threshold_values <- sim_output$baseline_threshold_values
participant_true_clusters <- sim_output$participant_true_clusters

# -------------------------------------------------------------------------
# Integrity checks
# -------------------------------------------------------------------------

expected_simulation_ids <- simulation_design$simulation_id
observed_simulation_ids <- unique(sim_results$simulation_id)

missing_simulation_ids <- setdiff(
    expected_simulation_ids,
    observed_simulation_ids
    )

failed_simulation_ids <- character()

if ("error" %in% names(sim_results) ) {
    
    failed_simulation_ids <- sim_results |>
        dplyr::filter(!is.na(error)) |>
        dplyr::pull(simulation_id) |>
        unique()
    
}

if (length(missing_simulation_ids) > 0 ||
    length(failed_simulation_ids) > 0) {

    message (
        "ERROR: incomplete array task. Missing cells = ",
        length(missing_simulation_ids),
        "; failed cells = ",
        length(failed_simulation_ids)
        )

    if (length(missing_simulation_ids) > 0) {
        
        message (
            "Missing simulation IDs: ",
            paste(missing_simulation_ids, collapse = ", ")
            )
        
    }

    if (length(failed_simulation_ids) > 0) {
        
        message (
            "Failed simulation IDs: ",
            paste(failed_simulation_ids, collapse = ", ")
            )
        
    }

    quit(save = "no", status = 1, runLast = FALSE)
    
}

message ("Method counts in sim_results:")

print(
    sim_results |> count(method, has_cluster),
    n = Inf
    )

if ("error" %in% names(sim_results) ) {
    
    method_errors <- sim_results |>
        filter(!is.na(error)) |>
        count(simulation_id, method, error)

    if (nrow(method_errors) > 0) {
        
        message ("Errors in sim_results:")
        print(method_errors, n = Inf)
        
    }
    
}

# -------------------------------------------------------------------------
# Save raw interval estimates and true-cluster metadata
# -------------------------------------------------------------------------

true_cluster_table <- map_dfr(
    seq_len(nrow(simulation_design) ),
    function (row_i) {
        make_true_cluster_table(simulation_design[row_i, ])
    }
)

result_prefix <- glue("errors_results_extended_hpc_{job_id}_{task_id}")

saveRDS(
    object = sim_results,
    file = glue("errors_results/{result_prefix}.rds")
    )

saveRDS(
    object = true_cluster_table,
    file = glue("errors_results/true_clusters_extended_hpc_{job_id}_{task_id}.rds")
    )

sim_results_long <- sim_results |>
    select(
        simulation_id,
        sim_index,
        effect_scenario,
        snr,
        true_n_clusters,
        threshold,
        method,
        baseline_threshold_name,
        effect_threshold,
        odds_threshold,
        cluster_id,
        has_cluster,
        onset,
        offset
        ) |>
    pivot_longer(
        cols = c(onset, offset),
        names_to = "onset_offset",
        values_to = "value"
        )

saveRDS(
    object = sim_results_long,
    file = glue("errors_results/errors_results_extended_long_hpc_{job_id}_{task_id}.rds")
    )

matched_results <- match_detected_to_true(
    sim_results = sim_results,
    true_clusters = true_cluster_table
    )

saveRDS(
    object = matched_results,
    file = glue("errors_results/errors_results_extended_matched_hpc_{job_id}_{task_id}.rds")
    )

summary_results <- matched_results |>
    summarise(
        n_detected = n(),
        prop_matched_to_true = mean(matched_true_cluster, na.rm = TRUE),
        median_temporal_iou = median(temporal_iou, na.rm = TRUE),
        median_onset_abs_error_ms = median(onset_abs_error * 1000, na.rm = TRUE),
        median_offset_abs_error_ms = median(offset_abs_error * 1000, na.rm = TRUE),
        .by = c(effect_scenario, snr, method, threshold)
        ) |>
    arrange(effect_scenario, snr, method, threshold)

saveRDS(
    object = summary_results,
    file = glue("errors_results/errors_results_extended_summary_hpc_{job_id}_{task_id}.rds")
    )

write.csv(
    x = summary_results,
    file = glue("errors_results/errors_results_extended_summary_hpc_{job_id}_{task_id}.csv"),
    row.names = FALSE
    )

# -------------------------------------------------------------------------
# Save baseline-derived BGAM thresholds and posterior-odds clusters
# -------------------------------------------------------------------------

baseline_threshold_clusters <- sim_results |>
    filter(!is.na(baseline_threshold_name) )

saveRDS(
    object = baseline_threshold_values,
    file = glue("errors_results/baseline_threshold_values_hpc_{job_id}_{task_id}.rds")
    )

write.csv(
    x = baseline_threshold_values,
    file = glue("errors_results/baseline_threshold_values_hpc_{job_id}_{task_id}.csv"),
    row.names = FALSE
    )

saveRDS(
    object = baseline_threshold_clusters,
    file = glue("errors_results/baseline_threshold_clusters_hpc_{job_id}_{task_id}.rds")
    )

write.csv(
    x = baseline_threshold_clusters,
    file = glue("errors_results/baseline_threshold_clusters_hpc_{job_id}_{task_id}.csv"),
    row.names = FALSE
    )

# Participant-specific ground-truth onset/offset values used in the simulator.
saveRDS(
    object = participant_true_clusters,
    file = glue("errors_results/participant_true_clusters_hpc_{job_id}_{task_id}.rds")
    )

write.csv(
    x = participant_true_clusters,
    file = glue("errors_results/participant_true_clusters_hpc_{job_id}_{task_id}.csv"),
    row.names = FALSE
    )

message ("Finished array task ", task_id)
message ("Saved files with job_id = ", job_id)
