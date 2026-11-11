####################################################################
# Split-half reliability analysis                                  #
# Loca (non-SLURM) version comparing 4 baseline-derived thresholds #  
# ---------------------------------------------------------------- #
# Written by Ladislas Nalborczyk                                   #
# Contact: ladislas.nalborczyk@cnrs.fr                             #
# Last updated on November 4, 2026                                 #
####################################################################

# -------------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------------

library(changepoint)
library(tidybayes)
library(posterior)
library(dplyr)
library(tidyr)
library(purrr)
library(glue)
library(brms)

# reticulate is needed only for cluster-mass / TFCE
use_mne_cluster_methods <- TRUE

if (use_mne_cluster_methods) {
    
    library (reticulate)
    
}

# -------------------------------------------------------------------------
# Simulation settings
# -------------------------------------------------------------------------

# Change default brms backend
options(brms.backend = "cmdstanr")

# Number of half splits
n_splits_total <- 1000L

# Parallelisation across split-halves
# At most n_workers BGAMMs are fitted simultaneously
n_workers <- 3L

# Resume from completed per-fit checkpoints after interruption/restart
resume_existing <- TRUE
split_seed <- 666L

chance_level <- 0.5
alpha_level <- 0.05

# Baseline used to define the four BGAM thresholds
baseline_window <- c(-Inf, -0.1)

# BGAM posterior-odds thresholds
post_prob_ratio_thresholds <- c(1, 3, 6, 10, 20, 50, 100)

# Model settings
brms_chains <- 4L
brms_cores <- 4L

# Sanity check: 3 workers x 4 Stan chains = 12 active CPU processes.
total_requested_cores <- n_workers * brms_cores

if (parallel::detectCores() < total_requested_cores) {
    
    warning (
        "Requested ", total_requested_cores,
        " concurrent cores but parallel::detectCores() reports only ",
        parallel::detectCores(), "."
        )
    
}

# Number of iterations per chain
brms_iter <- 3000L
brms_warmup <- 1000L

# Population and participant-varying smooth basis dimensions
population_k <- 30L
participant_k <- 10L

# NULL = use all posterior draws
n_epred_draws <- NULL

# Data / output paths
data_file <- "data/decoding_results_reshaped.csv"
utils_file <- "code/utils.R"
results_dir <- "reliability_results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# One small RDS checkpoint per fitted split-half
checkpoint_dir <- file.path(results_dir, "task_checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

# Keep CmdStan-generated model files/executables outside temporary directories
# This also makes recompilation less likely across persistent worker sessions
cmdstan_cache_dir <- file.path(results_dir, "cmdstan_cache")
dir.create(cmdstan_cache_dir, recursive = TRUE, showWarnings = FALSE)
options(cmdstanr_write_stan_file_dir = normalizePath(cmdstan_cache_dir) )

# Prevent hidden BLAS/OpenMP oversubscription.
Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
    )

# -------------------------------------------------------------------------
# Project helpers / optional MNE setup
# -------------------------------------------------------------------------

if (!file.exists(utils_file) ) {
    
    stop ("Could not find ", utils_file, ". Run this script from the project root.")
    
}

if (!file.exists(data_file) ) {
    
    stop ("Could not find ", data_file, ". Run this script from the project root.")
    
}

source (utils_file)

if (!exists("find_clusters", mode = "function") ) {
    
    stop ("find_clusters() was not loaded from ", utils_file)
    
}

if (use_mne_cluster_methods) {

    if (!exists("freq_stats_cluster_matrix", mode = "function") ) {
        
        stop ("freq_stats_cluster_matrix() was not loaded from ", utils_file)
        
    }

    # Resolve the Python executable once in the main process.
    # Each PSOCK worker will start its own reticulate/Python session.
    reticulate_conda_binary <- path.expand("~/miniforge3/bin/conda")
    
    reticulate_python <- reticulate::conda_python(
        envname = "r-reticulate3",
        conda = reticulate_conda_binary
        )

    if (!file.exists(reticulate_python) ) {
        
        stop ("Could not resolve Python executable for r-reticulate3.")
        
    }

} else {

    reticulate_python <- NA_character_
}

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

safe_ttest <- function (x, mu = chance_level) {

    x <- x[is.finite(x)]

    if (length(x) < 2L || sd(x) == 0) {
        
        return (
            tibble(
                tval = NA_real_,
                pval = 1
                )
            )
        
    }

    out <- t.test(
        x = x,
        mu = mu
        )

    return (
        tibble(
            tval = unname(out$statistic)^2,
            pval = out$p.value
            )
        )
    
}

make_participant_split <- function (
    participants,
    split_number,
    seed = split_seed
    ) {

    set.seed(seed + split_number)

    shuffled <- sample(participants)
    n_a <- floor(length(shuffled) / 2)

    tibble(
        split_number = split_number,
        participant = shuffled,
        group = if_else(
            seq_along(shuffled) <= n_a,
            "A",
            "B"
            )
        )
    
}

jaccard_set <- function (a, b) {

    a <- unique(as.character(a) )
    b <- unique(as.character(b) )

    union_n <- length(union(a, b) )

    if (union_n == 0L) {
        
        return (NA_real_)
        
    }

    return (length(intersect(a, b) ) / union_n)
    
}

compute_partition_similarity <- function (split_df) {

    split_ids <- sort(unique(split_df$split_number) )
    participant_ids <- sort(unique(as.character(split_df$participant) ) )

    if (length(split_ids) < 2L) {
        
        return (
            tibble(
                split_1 = integer(),
                split_2 = integer(),
                jaccard_same_labels = double(),
                jaccard_swapped_labels = double(),
                partition_jaccard = double()
                )
            )
        
    }

    # Binary membership matrices: rows = splits, columns = participants
    membership_a <- matrix(
        0L,
        nrow = length(split_ids),
        ncol = length(participant_ids),
        dimnames = list(
            as.character(split_ids),
            participant_ids
            )
        )

    a_rows <- split_df |>
        filter(group == "A") |>
        transmute(
            split_index = match(split_number, split_ids),
            participant_index = match(
                as.character(participant),
                participant_ids
                )
            )

    membership_a[
        cbind(
            a_rows$split_index,
            a_rows$participant_index
            )
        ] <- 1L
    
    membership_b <- 1L - membership_a

    n_a <- rowSums(membership_a)
    n_b <- rowSums(membership_b)

    # All pairwise intersection sizes at once
    int_aa <- tcrossprod(membership_a)
    int_bb <- tcrossprod(membership_b)
    int_ab <- membership_a %*% t(membership_b)
    int_ba <- membership_b %*% t(membership_a)

    pair_index <- which(
        upper.tri(int_aa),
        arr.ind = TRUE
        )

    i <- pair_index[, 1]
    j <- pair_index[, 2]

    int_aa_v <- int_aa[cbind(i, j)]
    int_bb_v <- int_bb[cbind(i, j)]
    int_ab_v <- int_ab[cbind(i, j)]
    int_ba_v <- int_ba[cbind(i, j)]

    jac_aa <- int_aa_v / (n_a[i] + n_a[j] - int_aa_v)

    jac_bb <- int_bb_v / (n_b[i] + n_b[j] - int_bb_v)

    jac_ab <- int_ab_v / (n_a[i] + n_b[j] - int_ab_v)

    jac_ba <- int_ba_v / (n_b[i] + n_a[j] - int_ba_v)

    same_labels <- (jac_aa + jac_bb) / 2
    swapped_labels <- (jac_ab + jac_ba) / 2

    tibble(
        split_1 = split_ids[i],
        split_2 = split_ids[j],
        jaccard_same_labels = same_labels,
        jaccard_swapped_labels = swapped_labels,
        # A/B labels are arbitrary
        partition_jaccard = pmax(
            same_labels,
            swapped_labels
            )
        )
    
}

make_data_matrix <- function (summary_df_split) {

    time_vec <- sort(unique(summary_df_split$time) )
    time_chr <- sprintf("%.10f", time_vec)

    wide_df <- summary_df_split |>
        mutate(time_chr = sprintf("%.10f", time) ) |>
        select(
            participant,
            time_chr,
            auc_mean
            ) |>
        pivot_wider(
            names_from = time_chr,
            values_from = auc_mean
            ) |>
        arrange(participant)

    list(
        X = as.matrix(wide_df[, time_chr, drop = FALSE]),
        time_vec = time_vec
        )
    
}

safe_find_clusters <- function (data, threshold) {

    tryCatch (
        {
            out <- find_clusters(
                data = data,
                threshold = threshold
            )

            if (is.null(out) || nrow(out) == 0L) {
                tibble()
            } else {
                out
            }
        },
        error = function (e) {

            warning (
                "find_clusters() failed: ",
                conditionMessage(e)
                )

            tibble()
            
        }
        
    )
    
}

format_threshold_for_method <- function (x) {

    gsub(
        pattern = "\\.",
        replacement = "p",
        x = formatC(
            x,
            digits = 12,
            format = "fg",
            flag = "#"
            )
        )
    
}

make_bgam_method_name <- function (
    threshold_method,
    odds_threshold
    ) {

    paste0(
        "brms_",
        threshold_method,
        "_odds_",
        format_threshold_for_method(odds_threshold)
        )
    
}

as_cluster_table <- function (
    cluster_df,
    split_number,
    group,
    method,
    method_family = NA_character_,
    threshold_method = NA_character_,
    effect_threshold = NA_real_,
    post_prob_ratio_threshold = NA_real_
    ) {

    if (is.null(cluster_df) || nrow(cluster_df) == 0L) {

        return (
            tibble(
                split_number = split_number,
                group = group,
                method = method,
                method_family = method_family,
                sesoi_method = threshold_method,
                bgam_threshold_method = threshold_method,
                post_prob_ratio_threshold = post_prob_ratio_threshold,
                bgam_odds_threshold = post_prob_ratio_threshold,
                cluster_id = NA_integer_,
                onset = NA_real_,
                offset = NA_real_,
                has_cluster = FALSE,
                sesoi = effect_threshold,
                effect_threshold = effect_threshold
                )
            )
        
    }

    cluster_df |>
        transmute(
            split_number = split_number,
            group = group,
            method = method,
            method_family = method_family,
            sesoi_method = threshold_method,
            bgam_threshold_method = threshold_method,
            post_prob_ratio_threshold = post_prob_ratio_threshold,
            bgam_odds_threshold = post_prob_ratio_threshold,
            cluster_id = as.integer(cluster_id),
            onset = cluster_onset,
            offset = cluster_offset,
            has_cluster = TRUE,
            sesoi = effect_threshold,
            effect_threshold = effect_threshold
            )
    
}

get_population_epred_draws <- function (
    gam_split,
    n_epred_draws = NULL
    ) {

    # Population-level predictions following the neurogam strategy:
    # include participant in newdata for validation but exclude
    # participant-level terms using re_formula = NA.
    newdata <- tibble(
        time = sort(unique(gam_split$data$time) ),
        participant = NA,
        auc_se = 1
        )

    if (is.null(n_epred_draws) ) {

        draws <- tidybayes::add_epred_draws(
            newdata = newdata,
            object = gam_split,
            re_formula = NA
            )

    } else {

        draws <- tidybayes::add_epred_draws(
            newdata = newdata,
            object = gam_split,
            re_formula = NA,
            ndraws = n_epred_draws
            )
        
    }

    draws |>
        data.frame() |>
        transmute(
            .draw = .draw,
            time = time,
            auc = .epred,
            effect = .epred - chance_level
            )
    
}

compute_bgam_thresholds <- function (
    epred_draws,
    baseline_window
    ) {

    baseline_draws <- epred_draws |>
        filter(
            time >= baseline_window[1],
            time < baseline_window[2]
            )

    if (nrow(baseline_draws) == 0L) {
        
        stop (
            "No posterior predictions found in baseline window: ",
            paste(baseline_window, collapse = ", ")
            )
        
    }

    baseline_summary <- baseline_draws |>
        summarise(
            effect_median = median(
                effect,
                na.rm = TRUE
                ),
            effect_q025 = quantile(
                effect,
                0.025,
                na.rm = TRUE
                ),
            effect_q975 = quantile(
                effect,
                0.975,
                na.rm = TRUE
                ),
            .by = time
            )

    thresholds <- tibble(
        threshold_method = c(
            "baseline_max_abs_qi95",
            "baseline_max_positive_qi95",
            "baseline_max_abs_median",
            "baseline_max_positive_median"
            ),
        threshold_label = c(
            "Max absolute 95% CrI endpoint",
            "Max positive 95% CrI endpoint",
            "Max absolute posterior median",
            "Max positive posterior median"
            ),
        effect_threshold = c(
            # 1. Max absolute endpoint of the pointwise 95% CrI
            max(
                abs(
                    c(
                        baseline_summary$effect_q025,
                        baseline_summary$effect_q975
                        )
                    ),
                na.rm = TRUE
                ),
            # 2. Max positive endpoint of the pointwise 95% CrI
            max(
                0,
                baseline_summary$effect_q975,
                na.rm = TRUE
                ),
            # 3. Max absolute posterior-median effect
            max(
                abs(baseline_summary$effect_median),
                na.rm = TRUE
                ),
            # 4. Max positive posterior-median effect.
            max(
                0,
                baseline_summary$effect_median,
                na.rm = TRUE
                )
            )
        )

    return (
        list(
            baseline_summary = baseline_summary,
            thresholds = thresholds
            )
        )
    
}

compute_posterior_odds <- function (
    epred_draws,
    effect_threshold
    ) {

    n_draws_eff <- n_distinct(epred_draws$.draw)
    eps <- 1 / n_draws_eff

    epred_draws |>
        summarise(
            posterior_probability = mean(
                effect > effect_threshold,
                na.rm = TRUE
                ),
            .by = time
            ) |>
        mutate(
            posterior_probability = pmin(
                pmax(
                    posterior_probability,
                    eps
                    ),
                1 - eps
                ),
            prob_ratio = posterior_probability / (1 - posterior_probability)
            )
    
}

run_bgam_threshold_grid <- function (
    epred_draws,
    threshold_table,
    split_number,
    group
    ) {

    pmap_dfr(
        threshold_table,
        function (
            threshold_method,
            threshold_label,
            effect_threshold
            ) {

            odds_df <- compute_posterior_odds(
                epred_draws = epred_draws,
                effect_threshold = effect_threshold
                )

            map_dfr(
                post_prob_ratio_thresholds,
                function (odds_threshold) {

                    method_name <-
                        make_bgam_method_name(
                            threshold_method =
                                threshold_method,
                            odds_threshold =
                                odds_threshold
                            )

                    safe_find_clusters(
                        data = odds_df |>
                            select(
                                time,
                                value = prob_ratio
                                ),
                        threshold = odds_threshold
                        ) |>
                        as_cluster_table(
                            split_number =
                                split_number,
                            group = group,
                            method = method_name,
                            method_family = "brms",
                            threshold_method = threshold_method,
                            effect_threshold = effect_threshold,
                            post_prob_ratio_threshold = odds_threshold
                            )
                    
                }
                
            )
            
        }
        
    )
    
}

# -------------------------------------------------------------------------
# Import and summarise empirical data
# -------------------------------------------------------------------------

raw_df <- read.csv(file = data_file) |>
    mutate(
        participant_id = cur_group_id(),
        .by = participant
        ) |>
    filter(participant_id < 33) |>
    select(-participant_id)

summary_df <- raw_df |>
    summarise(
        auc_mean = mean(
            auc,
            na.rm = TRUE
            ),
        auc_sd = sd(
            auc,
            na.rm = TRUE
            ),
        auc_se = auc_sd / sqrt(n() ),
        .by = c(
            participant,
            time
            )
        ) |>
    mutate( participant = factor(participant) )

participants <- sort(unique(as.character(summary_df$participant) ) )

message ("Participants included: ", length(participants) )

# -------------------------------------------------------------------------
# Generate and SAVE all split compositions
# -------------------------------------------------------------------------

split_df <- map_dfr(
    seq_len(n_splits_total),
    function (current_split_number) {

        make_participant_split(
            participants = participants,
            split_number = current_split_number,
            seed = split_seed
            )
        
    }
    
)

split_counts <- split_df |>
    count(
        split_number,
        group
        ) |>
    pivot_wider(
        names_from = group,
        values_from = n
        )

print(split_counts)

split_composition_summary <- split_df |>
    arrange(
        split_number,
        group,
        participant
        ) |>
    summarise(
        n_participants = n(),
        participants = paste(
            participant,
            collapse = ";"
            ),
        .by = c(split_number, group)
        )

split_partition_similarity <- compute_partition_similarity(split_df)

saveRDS(
    split_df,
    file.path(
        results_dir,
        "split_composition_long.rds"
        )
    )

write.csv(
    split_df,
    file.path(
        results_dir,
        "split_composition_long.csv"
        ),
    row.names = FALSE
    )

write.csv(
    split_composition_summary,
    file.path(
        results_dir,
        "split_composition_summary.csv"
        ),
    row.names = FALSE
    )

saveRDS(
    split_partition_similarity,
    file.path(
        results_dir,
        "split_partition_similarity.rds"
        )
    )

write.csv(
    split_partition_similarity,
    file.path(
        results_dir,
        "split_partition_similarity.csv"
        ),
    row.names = FALSE
    )

# -------------------------------------------------------------------------
# Fit one split-half
# -------------------------------------------------------------------------

fit_one_split_half <- function (
    current_split,
    current_group
    ) {

    message ("\n============================================================")

    message (
        "Running split ",
        current_split,
        ", group ",
        current_group
        )

    participants_split <- split_df |>
        filter(
            split_number == current_split,
            group == current_group
            ) |>
        pull(participant)

    summary_df_split <- summary_df |>
        filter(as.character(participant) %in% participants_split) |>
        mutate(participant = droplevels(participant) )

    # -----------------------------------------------------
    # BGAMM
    #
    # Population smooth + participant-varying smooths,
    # with known SEs AND an additional residual sigma
    # -----------------------------------------------------

    gam_split <- brm(
        formula =
            auc_mean | se(x = auc_se, sigma = TRUE) ~ 1 +
            (1 | participant) +
            s(time, bs = "cr", k = population_k) +
            s(time, participant, bs = "fs", m = 1, k = participant_k),
        data = summary_df_split,
        family = gaussian(),
        warmup = brms_warmup,
        iter = brms_iter,
        chains = brms_chains,
        cores = brms_cores,
        refresh = 0,
        silent = 2,
        seed =
            100000L +
            current_split * 10L +
            if_else(current_group == "A", 1L, 2L)
        )

    # Population-level posterior expected trajectory
    epred_draws <- get_population_epred_draws(
        gam_split = gam_split,
        n_epred_draws = n_epred_draws
        )

    threshold_obj <- compute_bgam_thresholds(
        epred_draws = epred_draws,
        baseline_window = baseline_window
        )

    threshold_table <-
        threshold_obj$thresholds |>
        mutate(
            split_number = current_split,
            group = current_group,
            .before = 1
            )

    onset_offset_brms <-
        run_bgam_threshold_grid(
            epred_draws = epred_draws,
            threshold_table = threshold_obj$thresholds,
            split_number = current_split,
            group = current_group
            )

    # -----------------------------------------------------
    # Mass univariate t-tests
    # -----------------------------------------------------

    tests_results <- summary_df_split |>
        summarise(
            safe_ttest(
                auc_mean,
                mu = chance_level
                ),
            .by = time
            ) |>
        mutate(
            pval_bh = p.adjust(
                pval,
                method = "BH"
                ),
            pval_by = p.adjust(
                pval,
                method = "BY"
                ),
            pval_holm = p.adjust(
                pval,
                method = "holm"
                )
            )

    onset_offset_raw_p <- safe_find_clusters(
        data = tests_results |>
            mutate(value = -pval) |>
            select(time, value),
        threshold = -alpha_level
        ) |>
        as_cluster_table(
            split_number = current_split,
            group = current_group, method = "raw_p"
            )
    
    onset_offset_pval_bh <- safe_find_clusters(
        data = tests_results |>
            mutate(value = -pval_bh) |>
            select(time, value),
        threshold = -alpha_level
        ) |>
        as_cluster_table(
            split_number = current_split,
            group = current_group,
            method = "pval_bh"
            )

    onset_offset_pval_by <- safe_find_clusters(
        data = tests_results |>
            mutate(value = -pval_by) |>
            select(time, value),
        threshold = -alpha_level
        ) |>
        as_cluster_table(
            split_number = current_split,
            group = current_group,
            method = "pval_by"
            )

    onset_offset_pval_holm <-
        safe_find_clusters(
            data = tests_results |>
                mutate(value = -pval_holm) |>
                select(time, value),
            threshold = -alpha_level
            ) |>
        as_cluster_table(
            split_number = current_split,
            group = current_group,
            method = "pval_holm"
            )

    # -----------------------------------------------------
    # Change-point comparator
    # -----------------------------------------------------

    onset_offset_cpt <- tryCatch (
        {

            tvals <- tests_results$tval
            tvals[!is.finite(tvals)] <- 0

            res <- cpt.meanvar(
                data = tvals,
                method = "BinSeg",
                Q = 2
                )

            cpts <- changepoint::cpts(res)

            if (length(cpts) < 2L) {

                tibble()

            } else {

                tibble(
                    cluster_id = 1L,
                    cluster_onset = tests_results$time[cpts[1]],
                    cluster_offset = tests_results$time[cpts[2]]
                    )
                
            }
        },
        
        error = function (e) {
            
            tibble()
            
        }
        
        ) |>
        
        as_cluster_table(
            split_number = current_split,
            group = current_group,
            method = "cpt"
            )

    # -----------------------------------------------------
    # Cluster mass / TFCE
    # -----------------------------------------------------

    if (use_mne_cluster_methods) {

        mat_obj <- make_data_matrix(summary_df_split)

        data_matrix <- mat_obj$X
        time_vec <- mat_obj$time_vec

        p_values_cluster_mass <-
            freq_stats_cluster_matrix(
                X = data_matrix - chance_level,
                timesteps = time_vec,
                cluster_type = "mass",
                alpha_level = alpha_level
                )

        p_values_tfce <-
            freq_stats_cluster_matrix(
                X = data_matrix - chance_level,
                timesteps = time_vec,
                cluster_type = "tfce",
                alpha_level = alpha_level
                )

        onset_offset_cluster_mass <-
            safe_find_clusters(
                data =
                    p_values_cluster_mass |>
                    mutate(value = -pval) |>
                    select(time, value),
                threshold = -alpha_level
                ) |>
            as_cluster_table(
                split_number = current_split,
                group = current_group,
                method = "cluster_mass"
                )

        onset_offset_cluster_tfce <-
            safe_find_clusters(
                data =
                    p_values_tfce |>
                    mutate(value = -pval) |>
                    select(time, value),
                threshold = -alpha_level
                ) |>
            as_cluster_table(
                split_number = current_split,
                group = current_group,
                method = "cluster_tfce"
                )

    } else {

        onset_offset_cluster_mass <- tibble()
        onset_offset_cluster_tfce <- tibble()
        
    }

    return (
        list(
            estimates = bind_rows(
                onset_offset_brms,
                onset_offset_raw_p,
                onset_offset_pval_bh,
                onset_offset_pval_by,
                onset_offset_pval_holm,
                onset_offset_cpt,
                onset_offset_cluster_mass,
                onset_offset_cluster_tfce
            ),
            bgam_thresholds = threshold_table
            )
        )
    
}

# -------------------------------------------------------------------------
# Checkpointed wrapper for one split-half task
# -------------------------------------------------------------------------

run_one_task <- function (i) {

    current_split <- split_group_combinations$split_number[i]
    current_group <- split_group_combinations$group[i]

    checkpoint_file <- file.path(
        checkpoint_dir,
        sprintf(
            "split_%04d_group_%s.rds",
            current_split,
            current_group
            )
        )
    
    if (resume_existing && file.exists(checkpoint_file) ) {

        message (
            "Loading checkpoint: split ",
            current_split,
            ", group ",
            current_group
            )

        return (readRDS(checkpoint_file) )
        
    }
    
    out <- tryCatch (
        {
            result <- fit_one_split_half(
                current_split = current_split,
                current_group = current_group
                )

            # Save only successful fits. Failed fits will be retried when
            # the script is run again
            saveRDS(
                result,
                checkpoint_file,
                compress = FALSE
                )

            resul
            t
        },
        
        error = function (e) {

            warning (
                "Failed split ",
                current_split,
                ", group ",
                current_group,
                ": ",
                conditionMessage(e)
                )

            list(
                estimates = tibble(
                    split_number = current_split,
                    group = current_group,
                    method = NA_character_,
                    method_family = NA_character_,
                    sesoi_method = NA_character_,
                    bgam_threshold_method = NA_character_,
                    post_prob_ratio_threshold = NA_real_,
                    bgam_odds_threshold = NA_real_,
                    cluster_id = NA_integer_,
                    onset = NA_real_,
                    offset = NA_real_,
                    has_cluster = FALSE,
                    sesoi = NA_real_,
                    effect_threshold = NA_real_,
                    error = conditionMessage(e)
                    ),
                bgam_thresholds = tibble()
                )
            
        }
        
    )

    return (out)
    
}

# -------------------------------------------------------------------------
# Three persistent PSOCK workers are used. Each worker runs one brms model
# at a time; that model runs four chains on four cores. Thus the maximum
# intended CPU use is 3 x 4 = 12 cores.
#
# parLapplyLB() performs load balancing: when one worker finishes a fit,
# it immediately receives the next unfinished split-half.
# -------------------------------------------------------------------------

split_group_combinations <- crossing(
    split_number = seq_len(n_splits_total),
    group = c("A", "B")
    )

message (
    "Split repetitions: ", n_splits_total,
    " -> brms fits: ", nrow(split_group_combinations)
    )

message (
    "Parallel layout: ", n_workers,
    " simultaneous fits x ", brms_cores,
    " cores per fit = ", total_requested_cores,
    " active cores maximum."
    )

project_dir <- normalizePath(".", mustWork = TRUE)

cl <- parallel::makePSOCKcluster(
    n_workers,
    outfile = ""
    )

# Minimal values needed to initialise worker sessions
parallel::clusterExport(
    cl,
    varlist = c(
        "project_dir",
        "utils_file",
        "use_mne_cluster_methods",
        "reticulate_python",
        "brms_cores",
        "cmdstan_cache_dir"
        ),
    envir = environment()
    )

parallel::clusterEvalQ(
    cl,
    {

        setwd(project_dir)

        library(changepoint)
        library(tidybayes)
        library(posterior)
        library(dplyr)
        library(tidyr)
        library(purrr)
        library(glue)
        library(brms)

        worker_cmdstan_cache <- file.path(
            cmdstan_cache_dir,
            paste0("worker_", Sys.getpid() )
            )

        dir.create(
            worker_cmdstan_cache,
            recursive = TRUE,
            showWarnings = FALSE
            )

        options(
            brms.backend = "cmdstanr",
            mc.cores = brms_cores,
            cmdstanr_write_stan_file_dir =
                normalizePath(
                    worker_cmdstan_cache,
                    mustWork = TRUE
                    )
            )

        Sys.setenv(
            OMP_NUM_THREADS = "1",
            MKL_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1",
            VECLIB_MAXIMUM_THREADS = "1",
            NUMEXPR_NUM_THREADS = "1"
            )

        source (utils_file)

        if (use_mne_cluster_methods) {

            library (reticulate)

            Sys.setenv(RETICULATE_PYTHON = reticulate_python)

            np <<- reticulate::import(
                "numpy",
                delay_load = FALSE
                )

            mne <<- reticulate::import(
                "mne",
                delay_load = FALSE
                )

            mne$set_log_level("WARNING")
            
        }

        NULL
        
    }
    
)

# Objects and functions required during each task.
worker_exports <- c(
    "summary_df",
    "split_df",
    "split_group_combinations",
    "chance_level",
    "alpha_level",
    "baseline_window",
    "post_prob_ratio_thresholds",
    "brms_chains",
    "brms_cores",
    "brms_iter",
    "brms_warmup",
    "population_k",
    "participant_k",
    "n_epred_draws",
    "use_mne_cluster_methods",
    "results_dir",
    "checkpoint_dir",
    "resume_existing",
    "safe_ttest",
    "make_data_matrix",
    "safe_find_clusters",
    "format_threshold_for_method",
    "make_bgam_method_name",
    "as_cluster_table",
    "get_population_epred_draws",
    "compute_bgam_thresholds",
    "compute_posterior_odds",
    "run_bgam_threshold_grid",
    "fit_one_split_half",
    "run_one_task"
    )

parallel::clusterExport(
    cl,
    varlist = worker_exports,
    envir = environment()
    )

task_indices <- seq_len(nrow(split_group_combinations) )

local_outputs <- tryCatch (
    parallel::parLapplyLB(
        cl,
        X = task_indices,
        fun = run_one_task
        ),
    finally = {
        try (
            parallel::stopCluster(cl),
            silent = TRUE
            )
        
    }
    
)

cl <- NULL

reliability_estimates <- map_dfr(local_outputs, "estimates")
bgam_thresholds <- map_dfr(local_outputs, "bgam_thresholds")

# -------------------------------------------------------------------------
# Save outputs
# -------------------------------------------------------------------------

saveRDS(
    reliability_estimates,
    file.path(
        results_dir,
        "split_half_estimates_local.rds"
        )
    )

write.csv(
    reliability_estimates,
    file.path(
        results_dir,
        "split_half_estimates_local.csv"
        ),
    row.names = FALSE
    )

saveRDS(
    bgam_thresholds,
    file.path(
        results_dir,
        "bgam_threshold_values_local.rds"
        )
    )

write.csv(
    bgam_thresholds,
    file.path(
        results_dir,
        "bgam_threshold_values_local.csv"
        ),
    row.names = FALSE
    )

bgam_cluster_estimates <-
    reliability_estimates |>
    filter(method_family == "brms")

saveRDS(
    bgam_cluster_estimates,
    file.path(
        results_dir,
        "split_half_bgam_cluster_grid_local.rds"
        )
    )

write.csv(
    bgam_cluster_estimates,
    file.path(
        results_dir,
        "split_half_bgam_cluster_grid_local.csv"
        ),
    row.names = FALSE
    )

message ("\nFinished.")
message ("Results directory: ", results_dir)
