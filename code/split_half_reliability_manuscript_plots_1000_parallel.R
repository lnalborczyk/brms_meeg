###########################################################
# Visualise the reliability results                       #
# ------------------------------------------------------- #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@cnrs.fr                    #
# Last updated on November 4, 2026                        #
###########################################################

library(tidyverse)
library(tidybayes)
library(patchwork)
library(MetBrewer)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

results_dir <- "reliability_results"
estimate_pattern <- "^split_half_estimates_.*\\.rds$"

time_step <- 0.005
include_offset <- TRUE
drop_zero_iou <- FALSE

# Select the BGAM posterior-odds threshold to display
bgam_sesoi_method_to_keep <- "baseline_max_abs_qi95"

# Available alternatives
bgam_sesoi_method_preference <- c(
    "baseline_max_abs_qi95",
    "baseline_max_positive_qi95",
    "baseline_max_abs_median",
    "baseline_max_positive_median"
    )

# Posterior odds threshold
bgam_odds_threshold_to_keep <- 10

method_levels <- c(
    "raw_p", "pval_bh", "pval_by", "pval_holm",
    "cluster_mass", "cluster_tfce", "cpt", "brms"
    )

method_labels <- c(
    "Raw p-value", "FDR BH95", "FDR BY01", "Holm",
    "Cluster mass", "TFCE", "Change point", "BGAM (our method)"
    )

method_label_lookup <- setNames(method_labels, method_levels)

method_cols <- setNames(
    met.brewer(name = "Johnson", n = length(method_labels) ),
    method_labels
    )

agreement_cols <- c(
    "Same number" = "grey35",
    "More in A" = "grey60",
    "More in B" = "grey82"
    )

# -----------------------------------------------------------------------------
# Import files
# -----------------------------------------------------------------------------

estimate_files <- list.files(
    path = results_dir,
    pattern = estimate_pattern,
    full.names = TRUE
    )

if (length(estimate_files) == 0) {
    
    stop ("No split_half_estimates_*.rds files found in: ", results_dir)
    
}

split_estimates <- purrr::map_dfr(
    estimate_files,
    \(f) {
        readRDS(f) |>
            mutate(source_file = basename(f) )
        }
    )

# Verify that the completed parallel run contains the expected split repetitions
n_loaded_splits <- dplyr::n_distinct(split_estimates$split_number)

message (
    "Loaded ",
    format(n_loaded_splits, big.mark = ","),
    " unique split repetitions from ",
    length(estimate_files),
    " estimate file(s)."
    )

if (n_loaded_splits < 1000L) {
    
    warning (
        "Only ", n_loaded_splits,
        " unique splits were found. The figure will use the available completed results."
        )
    
}

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

parse_bgam_method_name <- function (method) {
    
    method <- as.character(method)
    
    # Accept threshold encodings such as 10, 10p0, and 10.0.
    parsed <- stringr::str_match(
    method,
    paste0(
        "^brms_",
        "(baseline_max_abs_qi95|baseline_max_positive_qi95|baseline_max_abs_median|baseline_max_positive_median)",
        "_odds_([0-9]+(?:[p.][0-9]+)?)$"
        )
    )
    
    return (
        tibble(
            parsed_is_bgam = !is.na(parsed[, 1]),
            parsed_bgam_sesoi_method = parsed[, 2],
            parsed_bgam_odds_threshold = as.numeric(
                stringr::str_replace(parsed[, 3], "p", ".")
                )
            )
        )
    
}

first_existing_column <- function (data, candidates) {
    
    existing <- intersect(candidates, names(data) )
    if (length(existing) == 0) NA_character_ else existing[[1]]
    
}

median_iqr <- function (x) {
    
    x <- x[is.finite(x)]
    
    if (length(x) == 0) {
        
        return (
            tibble(
                n = 0,
                median_value = NA_real_,
                iqr_low = NA_real_,
                iqr_high = NA_real_
                )
            )
        
    }
    
    return (
        tibble(
            n = length(x),
            median_value = median(x),
            iqr_low = unname(quantile(x, probs = 0.25, names = FALSE) ),
            iqr_high = unname(quantile(x, probs = 0.75, names = FALSE) )
            )
        )
    
}

cluster_to_timestep_indices <- function (
    onset,
    offset,
    time_step = 0.001,
    include_offset = TRUE
    ) {
    
    onset <- as.numeric(onset)
    offset <- as.numeric(offset)
    
    if (
        length(onset) != 1 ||
        length(offset) != 1 ||
        !is.finite(onset) ||
        !is.finite(offset) ||
        offset < onset
        ) {
        
        return (integer() )
        
    }

    start_index <- ceiling(onset / time_step)

    if (include_offset) {
        
        end_index <- floor(offset / time_step)
        
    } else {
        
        end_index <- ceiling(offset / time_step) - 1L
    
    }
    
    if (!is.finite(start_index) || !is.finite(end_index) || end_index < start_index) {
        
        return (integer() )
        
    }
    
    return (
        seq.int(
            from = as.integer(start_index),
            to = as.integer(end_index),
            by = 1L
            )
        )
    
}

clusters_to_timestep_set <- function (
    onset,
    offset,
    time_step = 0.001,
    include_offset = TRUE
    ) {
    
    if (length(onset) == 0 || length(offset) == 0) {
        
        return (integer() )
        
    }
    
    purrr::map2(
        onset,
        offset,
        \(on, off) {
            cluster_to_timestep_indices(
                onset = on,
                offset = off,
                time_step = time_step,
                include_offset = include_offset
                )
            }
        ) |>
        unlist(use.names = FALSE) |>
        unique() |>
        sort() |>
        as.integer()
    
}

compute_timestep_iou <- function (timesteps_A, timesteps_B) {

    timesteps_A <- as.integer(timesteps_A)
    timesteps_B <- as.integer(timesteps_B)
    
    timesteps_A <- timesteps_A[!is.na(timesteps_A)]
    timesteps_B <- timesteps_B[!is.na(timesteps_B)]
    
    n_union <- length(union(timesteps_A, timesteps_B) )
    
    if (n_union == 0) {
        
        return (NA_real_)
        
    }
    
    n_intersection <- length(intersect(timesteps_A, timesteps_B) )
    
    return (n_intersection / n_union)
    
    }

# -----------------------------------------------------------------------------
# Check required columns
# -----------------------------------------------------------------------------

required_cols <- c(
    "split_number",
    "group",
    "method",
    "onset",
    "offset",
    "has_cluster"
    )

missing_cols <- setdiff(required_cols, names(split_estimates) )

if (length(missing_cols) > 0) {
    
    stop (
        "Missing required columns in split_estimates: ",
        paste(missing_cols, collapse = ", ")
        )
    
}

# -----------------------------------------------------------------------------
# Clean imported estimates and retain only the selected methods
# -----------------------------------------------------------------------------

# BGAM metadata may be encoded either in the method name or in separate columns
# saved in the RDS files. Resolve both representations before filtering
sesoi_column <- first_existing_column(
    split_estimates,
    c("bgam_threshold_method", "bgam_sesoi_method", "sesoi_method")
    )

odds_column <- first_existing_column(
    split_estimates,
    c("bgam_odds_threshold", "post_prob_ratio_threshold")
    )

existing_sesoi <- if (is.na(sesoi_column) ) {
    
    rep(NA_character_, nrow(split_estimates) )
    
    } else {
        
        as.character(split_estimates[[sesoi_column]])
        
    }

existing_odds <- if (is.na(odds_column) ) {
    
    rep(NA_real_, nrow(split_estimates) )
    
    } else {
        
        suppressWarnings(as.numeric(as.character(split_estimates[[odds_column]]) ) )
        
}

parsed_bgam <- parse_bgam_method_name(split_estimates$method)

split_estimates_tagged <- bind_cols(
    split_estimates,
    parsed_bgam,
    tibble(
        existing_bgam_sesoi_method = existing_sesoi,
        existing_bgam_odds_threshold = existing_odds
        )
    ) |>
    mutate(
        method = as.character(method),
        is_bgam = parsed_is_bgam | stringr::str_detect(method, "^brms(?:_|$)"),
        resolved_bgam_sesoi_method = coalesce(
            parsed_bgam_sesoi_method,
            existing_bgam_sesoi_method
            ),
    resolved_bgam_odds_threshold = coalesce(
        parsed_bgam_odds_threshold,
        existing_bgam_odds_threshold
        )
    )

bgam_candidates_threshold <- split_estimates_tagged |>
    filter(
        is_bgam,
        is.finite(resolved_bgam_odds_threshold),
        dplyr::near(
            resolved_bgam_odds_threshold,
            bgam_odds_threshold_to_keep
            )
        ) |>
    distinct(
        method,
        resolved_bgam_sesoi_method,
        resolved_bgam_odds_threshold
        )

if (nrow(bgam_candidates_threshold) == 0) {

    available_bgam <- split_estimates_tagged |>
        filter(is_bgam) |>
        distinct(
            method,
            resolved_bgam_sesoi_method,
            resolved_bgam_odds_threshold
            ) |>
        arrange(resolved_bgam_odds_threshold, method)
    
    stop (
        "No BGAM variant with posterior-odds threshold ",
        bgam_odds_threshold_to_keep,
        " was found. Available BGAM specifications:\n",
        paste(capture.output(print(available_bgam, n = Inf) ), collapse = "\n")
        )
    
}

available_sesoi <- bgam_candidates_threshold |>
    pull(resolved_bgam_sesoi_method) |>
    unique() |>
    na.omit() |>
    as.character()

if (!is.null(bgam_sesoi_method_to_keep) ) {
    
    if (!bgam_sesoi_method_to_keep %in% available_sesoi) {
        
        stop (
            "BGAM threshold ",
            bgam_odds_threshold_to_keep,
            " was found, but not with SESOI method '",
            bgam_sesoi_method_to_keep,
            "'. Available SESOI methods: ",
            paste(available_sesoi, collapse = ", ")
            )
        
    }
    
    selected_bgam_sesoi_method <- bgam_sesoi_method_to_keep
    
    } else {

        preferred_available <- bgam_sesoi_method_preference[
            bgam_sesoi_method_preference %in% available_sesoi
            ]

        if (length(preferred_available) > 0) {
            
            selected_bgam_sesoi_method <- preferred_available[[1]]
            
            } else if (length(available_sesoi) == 1) {
                
                selected_bgam_sesoi_method <- available_sesoi[[1]]
                
                } else if (length(available_sesoi) == 0) {
                    
                    # Some exports encode the threshold but not the SESOI definition
                    selected_bgam_sesoi_method <- NA_character_
                    
                    } else {
                        
                        stop (
                            "Several BGAM SESOI variants were found at threshold ",
                            bgam_odds_threshold_to_keep,
                            ". Set bgam_sesoi_method_to_keep explicitly. Available variants: ",
                            paste(available_sesoi, collapse = ", ")
                            )
                        
                    }
        
    }

message (
    "Using BGAM threshold ",
    bgam_odds_threshold_to_keep,
    if (is.na(selected_bgam_sesoi_method) ) {
        
        " (SESOI definition not encoded in the input files)."
        
        } else {
            
            paste0(" with SESOI method: ", selected_bgam_sesoi_method)
            
            }
    )

split_estimates_clean <- split_estimates_tagged |>
    mutate(
        source_file = as.character(source_file),
        split_number = as.integer(split_number),
        group = as.character(group),
        group = case_when(
            group %in% c("A", "a", "1", "first", "split_1", "split_A") ~ "A",
            group %in% c("B", "b", "2", "second", "split_2", "split_B") ~ "B",
            TRUE ~ group
            ),
        has_cluster = as.logical(has_cluster),
        onset = as.numeric(onset),
        offset = as.numeric(offset),
        keep_selected_bgam = is_bgam & is.finite(resolved_bgam_odds_threshold) &
            dplyr::near(resolved_bgam_odds_threshold, bgam_odds_threshold_to_keep) &
            if (is.na(selected_bgam_sesoi_method) ) {
                is.na(resolved_bgam_sesoi_method)
                } else {
                    resolved_bgam_sesoi_method == selected_bgam_sesoi_method
                    },
        method_plot = case_when(
            method %in% method_levels[method_levels != "brms"] ~ method,
            keep_selected_bgam ~ "brms",
            TRUE ~ NA_character_
            ),
        method_label = unname(method_label_lookup[method_plot]),
        method_sort_key = match(method_plot, method_levels)
        ) |>
    filter(
        group %in% c("A", "B"),
        !is.na(method_plot)
        )

# -----------------------------------------------------------------------------
# Robust cluster key
# -----------------------------------------------------------------------------

if (!"cluster_id" %in% names(split_estimates_clean) ) {
    
    split_estimates_clean <- split_estimates_clean |>
        group_by(source_file, split_number, method, group) |>
        mutate(
            cluster_id = NA_character_,
            cluster_key = paste0("row_", row_number() )
            ) |>
    ungroup()

    } else {

        split_estimates_clean <- split_estimates_clean |>
            group_by(source_file, split_number, method, group) |>
            mutate(
                cluster_id = as.character(cluster_id),
                cluster_key = if_else(
                    !is.na(cluster_id) & nzchar(cluster_id),
                    paste0("id_", cluster_id),
                    paste0("row_", row_number() )
                    )
                ) |>
            ungroup()
        
    }

# -----------------------------------------------------------------------------
# Build timestep sets for each split-half
# -----------------------------------------------------------------------------

split_grid <- split_estimates_clean |>
    distinct(
        source_file,
        split_number,
        method,
        method_plot,
        method_label,
        method_sort_key,
        group
        )

split_timestep_sets_detected <- split_estimates_clean |>
    filter(
        has_cluster,
        is.finite(onset),
        is.finite(offset),
        offset >= onset
        ) |>
    group_by(
        source_file,
        split_number,
        method,
        method_plot,
        method_label,
        method_sort_key,
        group
        ) |>
    summarise(
        has_detection = TRUE,
        n_clusters = as.integer(n_distinct(cluster_key) ),
        timestep_set = list(
            clusters_to_timestep_set(
                onset = onset,
                offset = offset,
                time_step = time_step,
                include_offset = include_offset
                )
            ),
        .groups = "drop"
        ) |>
    mutate(n_detected_timesteps = map_int(timestep_set, length) )

split_timestep_sets <- split_grid |>
    left_join(
        split_timestep_sets_detected,
        by = c(
            "source_file", "split_number", "method",
            "method_plot", "method_label", "method_sort_key", "group"
            )
        ) |>
    mutate(
        has_detection = replace_na(has_detection, FALSE),
    n_clusters = replace_na(n_clusters, 0L),
    timestep_set = map(
        timestep_set,
        \(x) {if (is.null(x) || length(x) == 0) integer() else as.integer(x)}
        ),
    n_detected_timesteps = map_int(timestep_set, length)
    )

# -----------------------------------------------------------------------------
# Pair split-halves and compute IoU
# -----------------------------------------------------------------------------

split_A <- split_timestep_sets |>
    filter(group == "A") |>
    select(
        source_file,
        split_number,
        method,
        method_plot,
        method_label,
        method_sort_key,
        has_detection_A = has_detection,
        n_clusters_A = n_clusters,
        n_detected_timesteps_A = n_detected_timesteps,
        timestep_set_A = timestep_set
        )

split_B <- split_timestep_sets |>
    filter(group == "B") |>
    select(
        source_file,
        split_number,
        method,
        has_detection_B = has_detection,
        n_clusters_B = n_clusters,
        n_detected_timesteps_B = n_detected_timesteps,
        timestep_set_B = timestep_set
        )

split_pairs_iou <- split_A |>
    full_join(
        split_B,
        by = c("source_file", "split_number", "method")
        ) |>
    mutate(
        has_detection_A = replace_na(has_detection_A, FALSE),
        has_detection_B = replace_na(has_detection_B, FALSE),
        n_clusters_A = replace_na(n_clusters_A, 0L),
        n_clusters_B = replace_na(n_clusters_B, 0L),
        n_detected_timesteps_A = replace_na(n_detected_timesteps_A, 0L),
        n_detected_timesteps_B = replace_na(n_detected_timesteps_B, 0L),
        timestep_set_A = map(
            timestep_set_A,
            \(x) {if (is.null(x) || length(x) == 0) integer() else as.integer(x)}
            ),
        timestep_set_B = map(
            timestep_set_B,
            \(x) {
                if (is.null(x) || length(x) == 0) integer() else as.integer(x)
                }
            ),
    both_detected = has_detection_A & has_detection_B,
    any_detected = has_detection_A | has_detection_B,
    neither_detected = !has_detection_A & !has_detection_B,
    detection_disagreement = xor(has_detection_A, has_detection_B),
    cluster_count_difference = n_clusters_A - n_clusters_B,
    abs_cluster_count_difference = abs(cluster_count_difference),
    same_cluster_count = n_clusters_A == n_clusters_B,
    more_clusters_A = n_clusters_A > n_clusters_B,
    more_clusters_B = n_clusters_B > n_clusters_A,
    n_intersection_timesteps = map2_int(
        timestep_set_A,
        timestep_set_B,
        \(x, y) length(intersect(x, y) )
        ),
    n_union_timesteps = map2_int(
        timestep_set_A,
        timestep_set_B,
        \(x, y) length(union(x, y) )
        ),
    timestep_iou = map2_dbl(
        timestep_set_A,
        timestep_set_B,
        compute_timestep_iou
        )
    )

# -----------------------------------------------------------------------------
# IoU summaries (median +/- IQR)
# -----------------------------------------------------------------------------

# Include all repetitions where at least one split-half detected something.
# If one half detected and the other did not, IoU = 0.
iou_long <- split_pairs_iou |>
    filter(any_detected, is.finite(timestep_iou) ) |>
    transmute(
        method_plot,
        method_label,
        method_sort_key,
        value = timestep_iou
        )

if (drop_zero_iou) {
    
    iou_long <- iou_long |> filter(value > 0)
  
}

iou_summary <- iou_long |>
    group_by(method_plot, method_label, method_sort_key) |>
    summarise(median_iqr(value), .groups = "drop")

# -----------------------------------------------------------------------------
# Cluster-count summaries
# -----------------------------------------------------------------------------

cluster_count_difference_summary <- split_pairs_iou |>
    group_by(method_plot, method_label, method_sort_key) |>
    summarise(
        n_pairs = n(),
        total_clusters_A = sum(n_clusters_A, na.rm = TRUE),
        total_clusters_B = sum(n_clusters_B, na.rm = TRUE),
        total_clusters_overall = sum(n_clusters_A + n_clusters_B, na.rm = TRUE),
        median_signed_difference = median(cluster_count_difference, na.rm = TRUE),
        iqr_low_signed_difference = unname(
            quantile(cluster_count_difference, probs = 0.25, na.rm = TRUE, names = FALSE)
            ),
        iqr_high_signed_difference = unname(
            quantile(cluster_count_difference, probs = 0.75, na.rm = TRUE, names = FALSE)
            ),
        n_same_cluster_count = sum(same_cluster_count, na.rm = TRUE),
        n_more_clusters_A = sum(more_clusters_A, na.rm = TRUE),
        n_more_clusters_B = sum(more_clusters_B, na.rm = TRUE),
        prop_same_cluster_count = mean(same_cluster_count, na.rm = TRUE),
        prop_more_clusters_A = mean(more_clusters_A, na.rm = TRUE),
        prop_more_clusters_B = mean(more_clusters_B, na.rm = TRUE),
        .groups = "drop"
        )

cluster_count_agreement_long <- cluster_count_difference_summary |>
    select(
        method_label,
        method_sort_key,
        prop_same_cluster_count,
        prop_more_clusters_A,
        prop_more_clusters_B
        ) |>
    pivot_longer(
        cols = starts_with("prop_"),
        names_to = "comparison",
        values_to = "proportion"
        ) |>
    mutate(
        comparison = recode(
          comparison,
          prop_same_cluster_count = "Same number",
          prop_more_clusters_A = "More in A",
          prop_more_clusters_B = "More in B"
          ),
    comparison = factor(
          comparison,
          levels = c("Same number", "More in A", "More in B")
          )
    )

cluster_difference_long <- split_pairs_iou |>
    transmute(
        method_plot,
        method_label,
        method_sort_key,
        cluster_count_difference
        )

# One observation per individual split-half. Pooling A and B here means that
# the mean below is the average number of clusters detected in one split-half,
# rather than the total number detected across a paired repetition.
cluster_counts_long <- split_pairs_iou |>
    select(
        source_file,
        split_number,
        method_plot,
        method_label,
        method_sort_key,
        n_clusters_A,
        n_clusters_B
        ) |>
    pivot_longer(
        cols = c(n_clusters_A, n_clusters_B),
        names_to = "split_half",
        names_pattern = "n_clusters_([AB])",
        values_to = "n_clusters"
        ) |>
    mutate(
        split_half = factor(split_half, levels = c("A", "B") ),
        n_clusters = as.integer(n_clusters)
        )

average_cluster_count_summary <- cluster_counts_long |>
    group_by(method_plot, method_label, method_sort_key) |>
    summarise(
        n_split_halves = sum(!is.na(n_clusters) ),
        mean_clusters = mean(n_clusters, na.rm = TRUE),
        sd_clusters = sd(n_clusters, na.rm = TRUE),
        sem_clusters = sd_clusters / sqrt(n_split_halves),
        sem_low_clusters = pmax(0, mean_clusters - sem_clusters),
        sem_high_clusters = mean_clusters + sem_clusters,
        .groups = "drop"
        )

# -----------------------------------------------------------------------------
# Plot-specific ordering and data formatting
# -----------------------------------------------------------------------------

# Match the ordering used in the original figures
iou_method_order <- iou_summary |>
    arrange(desc(median_value), method_sort_key) |>
    pull(method_label) |>
    as.character()

cluster_count_method_order <- cluster_count_difference_summary |>
    arrange(desc(total_clusters_overall), method_sort_key) |>
    pull(method_label) |>
    as.character()

average_cluster_method_order <- average_cluster_count_summary |>
    arrange(desc(mean_clusters), method_sort_key) |>
    pull(method_label) |>
    as.character()

iou_long_plot <- iou_long |>
    mutate(method_label = factor(method_label, levels = iou_method_order) )

iou_summary_plot <- iou_summary |>
    mutate(method_label = factor(method_label, levels = iou_method_order) )

cluster_count_agreement_plot <- cluster_count_agreement_long |>
    mutate(method_label = factor(method_label, levels = cluster_count_method_order) )

cluster_difference_long_plot <- cluster_difference_long |>
  mutate(method_label = factor(method_label, levels = cluster_count_method_order) )

cluster_difference_summary_plot <- cluster_count_difference_summary |>
  mutate(method_label = factor(method_label, levels = cluster_count_method_order) )

cluster_counts_long_plot <- cluster_counts_long |>
  mutate(method_label = factor(method_label, levels = average_cluster_method_order) )

average_cluster_count_summary_plot <- average_cluster_count_summary |>
  mutate(method_label = factor(method_label, levels = average_cluster_method_order) )

# -----------------------------------------------------------------------------
# Plot A: temporal overlap reliability (median +/- IQR)
# -----------------------------------------------------------------------------

p_timestep_iou_classical_best_bgam <- ggplot(
    iou_long_plot,
    aes(
        x = value,
        y = method_label,
        fill = method_label
        )
    ) +
    stat_halfeye(
        geom = "slab",
        side = "top",
        scale = 0.5,
        normalize = "groups",
        colour = NA,
        slab_alpha = 0.8,
        show.legend = FALSE
        ) +
    geom_segment(
        data = iou_summary_plot,
        aes(
            x = iqr_low,
            xend = iqr_high,
            y = method_label,
            yend = method_label
            ),
        inherit.aes = FALSE,
        linewidth = 0.75,
        colour = "black"
        ) +
    geom_point(
        data = iou_summary_plot,
        aes(
            x = median_value,
            y = method_label
            ),
        inherit.aes = FALSE,
        size = 2,
        colour = "black"
        ) +
    geom_label(
        data = iou_summary_plot,
        aes(
            x = median_value,
            y = method_label,
            label = sprintf("%.2f", median_value)
            ),
        inherit.aes = FALSE,
        size = 3,
        label.size = 0.75
        ) +
    scale_fill_manual(values = method_cols, breaks = method_labels, drop = FALSE) +
    scale_x_continuous(
        limits = c(0, 1.05),
        breaks = seq(0, 1, by = 0.25),
        expand = expansion(mult = c(0.02, 0.02) )
        ) +
    labs(
    x = "Split-half intersection-over-union (IoU)",
    y = NULL,
    title = ""
    )

# -----------------------------------------------------------------------------
# Plot B: agreement in cluster count
# -----------------------------------------------------------------------------

p_cluster_count_agreement <- ggplot(
    cluster_count_agreement_plot  %>% dplyr::filter(method_label != "Change point"),
    aes(
    x = proportion,
    y = method_label,
    fill = comparison
    ) ) +
    geom_col() +
    scale_fill_manual(values = agreement_cols, drop = FALSE) +
    scale_x_continuous(
        limits = c(0, 1),
        labels = scales::percent_format(accuracy = 1),
        expand = expansion(mult = c(0, 0.01) )
        ) +
    labs(
        x = "Proportion of split-half repetitions",
        y = NULL,
        fill = NULL,
        title = ""
        ) +
    theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# Plot C: signed cluster-count differences (median +/- IQR)
# -----------------------------------------------------------------------------

p_cluster_count_difference <- ggplot(
    cluster_difference_long_plot,
    aes(
        x = cluster_count_difference,
        y = method_label,
        fill = method_label
        ) ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.6,
        colour = "grey40"
        ) +
    stat_halfeye(
        geom = "slab",
        side = "top",
        scale = 0.5,
        normalize = "groups",
        colour = NA,
        slab_alpha = 0.8,
        show.legend = FALSE
        ) +
    geom_segment(
        data = cluster_difference_summary_plot,
        aes(
            x = iqr_low_signed_difference,
            xend = iqr_high_signed_difference,
            y = method_label,
            yend = method_label
            ),
        inherit.aes = FALSE,
        linewidth = 0.75,
        colour = "black"
        ) +
    geom_point(
    data = cluster_difference_summary_plot,
    aes(
        x = median_signed_difference,
        y = method_label
        ),
    inherit.aes = FALSE,
    size = 2,
    colour = "black"
    ) +
    scale_fill_manual(values = method_cols, breaks = method_labels, drop = FALSE) +
    labs(
        x = "Difference in cluster count (A - B)",
        y = NULL,
        title = ""
        )

# -----------------------------------------------------------------------------
# Plot D: average number of clusters detected per split-half
# -----------------------------------------------------------------------------

p_average_clusters_per_split_half <- ggplot(
    cluster_counts_long_plot %>% dplyr::filter(method_plot != "cpt"),
    aes(
        x = n_clusters,
        y = method_label,
        fill = method_label
        ) ) +
    stat_halfeye(
        geom = "slab",
        side = "top",
        scale = 0.5,
        normalize = "groups",
        colour = NA,
        slab_alpha = 0.8,
        show.legend = FALSE
        ) +
    geom_segment(
        data = average_cluster_count_summary_plot %>% dplyr::filter(method_plot != "cpt"),
        aes(
            x = sem_low_clusters,
            xend = sem_high_clusters,
            y = method_label,
            yend = method_label
            ),
        inherit.aes = FALSE,
        linewidth = 0.75,
        colour = "black"
        ) +
    geom_label(
    data = average_cluster_count_summary_plot %>% dplyr::filter(method_plot != "cpt"),
    aes(
        x = mean_clusters,
        y = method_label,
        label = sprintf("%.2f", mean_clusters)
        ),
    inherit.aes = FALSE,
    size = 3,
    label.size = 0.75
    ) +
    scale_fill_manual(values = method_cols, breaks = method_labels, drop = FALSE) +
    scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 2, 5, 10, 20, 50),
    expand = expansion(mult = c(0.05, 0.05) )
    ) +
    labs(
        x = "Number of detected clusters per split-half",
        y = NULL,
        title = ""
        )

cluster_variability_summary <- cluster_counts_long |>
    group_by(
        method_plot,
        method_label,
        method_sort_key
        ) |>
    summarise(n_split_halves = sum(!is.na(n_clusters) ),
    mean_clusters = mean(n_clusters, na.rm = TRUE),
    variance_clusters = var(n_clusters, na.rm = TRUE),
    sd_clusters = sd(n_clusters, na.rm = TRUE),
    # Variance-to-mean ratio: recommended general count-variability index
    dispersion_index = if_else(
        mean_clusters > 0,
        variance_clusters / mean_clusters,
        NA_real_
        ),
    # Relative SD; less stable when the mean is close to zero
    coefficient_variation = if_else(
        mean_clusters > 0,
        sd_clusters / mean_clusters,
        NA_real_
        ),
    # Robust measure of absolute variability
    mad_clusters = mad(
        n_clusters,
        center = median(n_clusters, na.rm = TRUE),
        constant = 1,
        na.rm = TRUE
        ),
    .groups = "drop"
    ) |>
    arrange(dispersion_index)

dispersion_method_order <- cluster_variability_summary |>
    arrange(desc(dispersion_index), method_sort_key) |>
    pull(method_label) |>
    as.character()

cluster_variability_summary_plot <- cluster_variability_summary |>
    mutate(
    method_label = factor(
        method_label,
        levels = dispersion_method_order
        )
    )

p_cluster_dispersion <- ggplot(
    cluster_variability_summary_plot %>% dplyr::filter(method_plot != "cpt"),
    aes(
    x = dispersion_index,
    y = method_label,
    fill = method_label
    ) ) +
    geom_col(
        width = 0.65,
        show.legend = FALSE
        ) +
    scale_fill_manual(
        values = method_cols,
        drop = FALSE
        ) +
    labs(
        x = "Dispersion index (variance / mean)",
        y = NULL
        )

# -----------------------------------------------------------------------------
# Combined patchwork figures for the manuscript
# -----------------------------------------------------------------------------

# Three-panel compact version
p_split_half_reliability_combined <- (
    p_timestep_iou_classical_best_bgam /
        (p_cluster_count_agreement | p_cluster_dispersion)
    ) +
    plot_layout(heights = c(1.1, 1) ) +
    plot_annotation(tag_levels = "A")

# Two-panel version
p_split_half_reliability_combined <- (
    p_timestep_iou_classical_best_bgam | p_cluster_count_agreement
    ) +
    plot_layout(heights = c(1, 1) ) +
    plot_annotation(tag_levels = "A")

p_split_half_reliability_combined <-
    p_split_half_reliability_combined &
    theme(
        plot.tag = element_text(
            face = "bold",
            size = 14
            )
        )

# -----------------------------------------------------------------------------
# Save manuscript figures
# -----------------------------------------------------------------------------

analysis_dir <- file.path(
    results_dir,
    "manuscript_figures"
    )

dir.create(
    analysis_dir,
    recursive = TRUE,
    showWarnings = FALSE
    )

ggsave(
    filename = file.path(
        analysis_dir,
        "split_half_reliability_1000_parallel.png"
        ),
    plot = p_split_half_reliability_combined,
    width = 14,
    height = 10,
    dpi = 300,
    bg = "white"
    )
