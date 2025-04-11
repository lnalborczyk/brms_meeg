###########################################################
# Split-half reliability of onset/offset estimates        #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on April 11, 2025                          #
###########################################################

# importing R packages
library(changepoint)
library(reticulate)
library(tidybayes)
library(dplyr)
library(tidyr)
library(brms)

# for setting up the cluster
library(doParallel)
library(foreach)

# timing the simulations
library(tictoc)
tic()

# total number of total CPU cores available on the HPC node
total_cores <- 12

# number of cores per brms model
cores_per_model <- 4

# number of parallel simulations to run at once
n_parallel <- total_cores / cores_per_model

# registering parallel backend
cl <- makeCluster(n_parallel)
registerDoParallel(cl)

# importing home-made helper functions
source("code/utils.R")

# importing the reshaped MEG decoding results
raw_df <- read.csv(file = "data/decoding_results_reshaped.csv") |>
    # removing the last participant (to get a round number of participants)
    mutate(participant_id = cur_group_id(), .by = participant) |>
    filter(participant_id < 33) |>
    select(-participant_id)

# averaging across participants
summary_df <- raw_df |>
    summarise(
        auc_mean = mean(auc),
        auc_sd = sd(auc),
        auc_se = sd(auc) / sqrt(n() ),
        .by = c(participant, time)
        )

# computing the smallest effect size of interest (SESOI)
# as the SD of decoding performance during the baseline
full_sesoi <- summary_df %>%
    filter(time < 0) %>%
    summarise(auc_mean = mean(auc_mean), .by = time) %>%
    summarise(sesoi = sd(auc_mean) ) %>%
    pull(sesoi)

# getting unique participants
unique_participants <- unique(raw_df$participant)

# getting the number of participants
n_participants <- as.numeric(length(unique_participants) )

# defining the number of splits
n_splits <- 1000

# number of posterior samples to use to compute the posterior probability ratio
# n_post_samples <- 1e4

# creating a dataframe where each participant appears in every split
split_df <- crossing(participant = unique_participants, split_number = 1:n_splits)

# assigning alternating groups (A or B) within each split for each participant
split_df <- split_df |>
    group_by(split_number) |>
    # shuffling groups for each split
    mutate(group = sample(rep(c("A", "B"), length.out = n() ) ) ) |>
    ungroup()

# defining the decoding chance level
chance_level <- 0.5

# defining the threshold on the posterior probability ratio
post_prob_ratio_threshold <- 20

# defining the significance level
alpha_level <- 0.05

# for each simulation
reliability_results <- tryCatch({foreach(
    i = seq_len(n_splits), .combine = bind_rows,
    .packages = c("reticulate", "brms", "dplyr", "tidyr", "changepoint", "tidybayes")
    ) %dopar% {
        
        # importing Python modules locally
        message_parallel("\nImporting Python modules...")
        use_condaenv("r-reticulate3", conda = "~/miniforge3/bin/conda", required = TRUE)
        
        # importing Python modules on the HPC cluster
        # py_require(packages = c("numpy", "mne"), python_version = "3.12.0")
        np <- import("numpy")
        mne <- import("mne")
        # mne$set_log_level("WARNING") # or "ERROR"
        
        # printing progress
        message_parallel(sprintf("\nRunning simulation %d out of %d...", i, n_splits) )
        
        # now for each split, assess the onset and offset
        participants_split <- split_df |>
            filter(split_number == i & group == "A") |>
            pull(participant)
        
        # splitting data
        raw_df_split <- raw_df |> filter(participant %in% participants_split)
        summary_df_split <- summary_df |> filter(participant %in% participants_split)
        
        # computing the smallest effect size of interest (SESOI)
        # as the SD of decoding performance during the baseline
        sesoi <- summary_df_split |>
            filter(time < 0) |>
            summarise(auc_mean = mean(auc_mean), .by = time) |>
            summarise(sesoi = sd(auc_mean) ) |>
            pull(sesoi)
        
        # fitting the GAM
        gam_split <- brm(
            auc_mean ~ s(time, bs = "cr", k = 50),
            data = summary_df_split,
            family = Beta(),
            warmup = 2000,
            iter = 5000,
            chains = cores_per_model,
            cores = cores_per_model
            # backend = "cmdstanr",
            # stan_model_args = list(stanc_options = list("O1") )
            )
        
        # computing the posterior probability
        message_parallel("\nModel fitted...")
        prob_y_above_0 <- data.frame(time = unique(gam_split$data$time) )|>
            # using a subset of posterior samples
            # add_epred_draws(object = gam_split, ndraws = n_post_samples) |>
            # or using all posterior samples
            add_epred_draws(object = gam_split)|>
            # converting to dataframe
            data.frame()|>
            group_by(time)|>
            summarise(m = mean(.epred > (0 + chance_level + sesoi) ) )|>
            mutate(prob_ratio = m / (1 - m) )|>
            ungroup()|>
            # ensuring there is no 0 or +Inf values
            mutate(prob_ratio = ifelse(is.infinite(prob_ratio), ndraws(gam_split), prob_ratio) )|>
            mutate(prob_ratio = ifelse(prob_ratio == 0, 1 / ndraws(gam_split), prob_ratio) )
            # mutate(prob_ratio = ifelse(is.infinite(prob_ratio), n_post_samples, prob_ratio) )|>
            # mutate(prob_ratio = ifelse(prob_ratio == 0, 1 / n_post_samples, prob_ratio) ) 
        
        # finding clusters
        onset_offset_brms <- find_clusters(
            data = prob_y_above_0 |> select(time, value = prob_ratio),
            threshold = post_prob_ratio_threshold
            ) |>
            mutate(split_id = formatC(x = i, width = 3, flag = 0) ) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "brms") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        # converting data to a matrix (for later use in MNE functions)
        data_matrix <- matrix(
            summary_df_split$auc_mean,
            ncol = length(unique(summary_df_split$time) ),
            byrow = TRUE
            )
        
        # massive univariate t-tests
        tests_results <- summary_df_split|>
            group_by(time) |>
            summarise(
                tval = t.test(x = auc_mean, mu = chance_level)$statistic^2,
                pval = t.test(x = auc_mean, mu = chance_level)$p.value
                ) |>
            mutate(
                pval_bh = p.adjust(p = pval, method = "BH"),
                pval_by = p.adjust(p = pval, method = "BY"),
                pval_holm = p.adjust(p = pval, method = "holm")
                ) |>
            ungroup()
            
        # using the binary segmentation method to identify onsets and offsets
        res <- cpt.meanvar(data = tests_results$tval, method = "BinSeg", Q = 2)
        onset_offset_cpt <- data.frame(
            split_id = formatC(x = i, width = 3, flag = 0),
            cluster_id = 1,
            onset = tests_results$time[res@cpts[1]],
            offset = tests_results$time[res@cpts[2]]
            ) |>
            mutate(method = "cpt") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        # finding clusters
        onset_offset_raw_p <- find_clusters(
            data = tests_results|> mutate(pval = pval * (-1) )|> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = formatC(x = i, width = 3, flag = 0) ) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "raw_p") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        onset_offset_pval_bh <- find_clusters(
            data = tests_results |> mutate(pval = pval_bh * (-1) ) |> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = formatC(x = i, width = 3, flag = 0) ) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "pval_bh") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        onset_offset_pval_by <- find_clusters(
            data = tests_results |> mutate(pval = pval_by * (-1) ) |> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = formatC(x = i, width = 3, flag = 0) ) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "pval_by") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        onset_offset_pval_holm <- find_clusters(
            data = tests_results |> mutate(pval = pval_holm * (-1) ) |> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = formatC(x = i, width = 3, flag = 0) ) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "pval_holm") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        # running the MNE cluster-based permutation
        p_values_cluster_mass <- freq_stats_cluster_matrix(
            X = data_matrix - chance_level,
            timesteps = unique(summary_df_split$time),
            cluster_type = "mass",
            alpha_level = alpha_level
            )
        p_values_tfce <- freq_stats_cluster_matrix(
            X = data_matrix - chance_level,
            timesteps = unique(summary_df_split$time),
            cluster_type = "tfce",
            alpha_level = alpha_level
            )
        
        onset_offset_cluster_mass <- find_clusters(
            data = p_values_cluster_mass |> mutate(pval = pval*(-1) ) |> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = formatC(x = i, width = 3, flag = 0) ) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "cluster_mass") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        onset_offset_cluster_tfce <- find_clusters(
            data = p_values_tfce |> mutate(pval = pval*(-1) ) |> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = formatC(x = i, width = 3, flag = 0) ) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "cluster_tfce") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        # binding all results together
        temp_reliability_results <- bind_rows(
            onset_offset_brms, onset_offset_raw_p, onset_offset_pval_bh,
            onset_offset_pval_by, onset_offset_pval_holm, onset_offset_cpt,
            onset_offset_cluster_mass, onset_offset_cluster_tfce
            )
        
        # returning the results
        temp_reliability_results
        
    }
    }, error = function(e) {
        message("Error in foreach: ", e$message)
    })

# importing Python modules
cat("\nSplits done. Importing Python modules...")
use_condaenv("r-reticulate3", conda = "~/miniforge3/bin/conda", required = TRUE)
# py_require(packages = c("numpy", "mne"), python_version = "3.12.0")
np <- import("numpy")
mne <- import("mne")

# fitting the final GAM
cat("\nNow fitting the model on the full dataset...")
full_gam <- brm(
    auc_mean ~ s(time, bs = "cr", k = 50),
    data = summary_df,
    family = Beta(),
    warmup = 2000,
    iter = 5000,
    chains = cores_per_model,
    cores = cores_per_model
    # backend = "cmdstanr",
    # stan_model_args = list(stanc_options = list("O1") )
    )

# sanity check
cat("\nNumber of posterior samples in the full model:", ndraws(full_gam) )

# computing the posterior probability
prob_y_above_0 <- data.frame(time = unique(full_gam$data$time) ) |>
    # retrieving the posterior samples
    # add_epred_draws(object = full_gam) |>
    # using a subset of posterior samples
    add_epred_draws(object = full_gam, ndraws = n_post_samples) |>
    # converting to dataframe
    data.frame() |>
    # computing mean posterior probability at the group level
    group_by(time) |>
    summarise(m = mean(.epred > (0 + chance_level + full_sesoi) ) ) |>
    mutate(prob_ratio = m / (1 - m) ) |>
    # ensuring there is no 0 or +Inf values
    # mutate(prob_ratio = ifelse(is.infinite(prob_ratio), ndraws(full_gam), prob_ratio) ) |>
    # mutate(prob_ratio = ifelse(prob_ratio == 0, 1 / ndraws(full_gam), prob_ratio) ) |>
    mutate(prob_ratio = ifelse(is.infinite(prob_ratio), n_post_samples, prob_ratio) ) |>
    mutate(prob_ratio = ifelse(prob_ratio == 0, 1 / n_post_samples, prob_ratio) ) |>
    ungroup()

# finding clusters
onset_offset_brms <- find_clusters(
    data = prob_y_above_0 |> select(time, value = prob_ratio),
    threshold = post_prob_ratio_threshold
    ) |>
    mutate(split_id = formatC(x = n_splits+1, width = 3, flag = 0) ) |>
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
    mutate(method = "brms") |>
    pivot_longer(cols = onset:offset) |>
    select(split_id, method, cluster_id, onset_offset = name, value) |>
    data.frame()

# converting data to a matrix (for later use in MNE functions)
data_matrix <- matrix(
    summary_df$auc_mean,
    ncol = length(unique(summary_df$time) ),
    byrow = TRUE
    )

# massive univariate t-tests
tests_results <- summary_df |>
    group_by(time) |>
    summarise(
        tval = t.test(x = auc_mean, mu = chance_level)$statistic^2,
        pval = t.test(x = auc_mean, mu = chance_level)$p.value
        ) |>
    mutate(
        pval_bh = p.adjust(p = pval, method = "BH"),
        pval_by = p.adjust(p = pval, method = "BY"),
        pval_holm = p.adjust(p = pval, method = "holm")
        ) |>
    ungroup()

# using the binary segmentation method to identify onsets and offsets
# res <- cpt.meanvar(data = tests_results$tval, method = "BinSeg", Q = 10)
# res@cpts
# tests_results$time[res@cpts]
# onset_offset_cluster_tfce
# res@cpttype
# res@cpts.full
res <- cpt.meanvar(data = tests_results$tval, method = "BinSeg", Q = 2)

onset_offset_cpt <- data.frame(
    split_id = formatC(x = n_splits+1, width = 3, flag = 0),
    cluster_id = 1,
    onset = tests_results$time[res@cpts[1]],
    offset = tests_results$time[res@cpts[2]]
    ) |>
    mutate(method = "cpt") |>
    pivot_longer(cols = onset:offset) |>
    select(split_id, method, cluster_id, onset_offset = name, value) |>
    data.frame()

# finding clusters
onset_offset_raw_p <- find_clusters(
    data = tests_results |> mutate(pval = pval * (-1) ) |> select(time, value = pval),
    threshold = -alpha_level
    ) |>
    mutate(split_id = formatC(x = n_splits+1, width = 3, flag = 0) ) |>
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
    mutate(method = "raw_p") |>
    pivot_longer(cols = onset:offset) |>
    select(split_id, method, cluster_id, onset_offset = name, value) |>
    data.frame()

onset_offset_pval_bh <- find_clusters(
    data = tests_results |> mutate(pval = pval_bh * (-1) ) |> select(time, value = pval),
    threshold = -alpha_level
    ) |>
    mutate(split_id = formatC(x = n_splits+1, width = 3, flag = 0) ) |>
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
    mutate(method = "pval_bh") |>
    pivot_longer(cols = onset:offset) |>
    select(split_id, method, cluster_id, onset_offset = name, value) |>
    data.frame()

onset_offset_pval_by <- find_clusters(
    data = tests_results |> mutate(pval = pval_by * (-1) ) |> select(time, value = pval),
    threshold = -alpha_level
    ) |>
    mutate(split_id = formatC(x = n_splits+1, width = 3, flag = 0) ) |>
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
    mutate(method = "pval_by") |>
    pivot_longer(cols = onset:offset) |>
    select(split_id, method, cluster_id, onset_offset = name, value) |>
    data.frame()

onset_offset_pval_holm <- find_clusters(
    data = tests_results |> mutate(pval = pval_holm * (-1) ) |> select(time, value = pval),
    threshold = -alpha_level
    ) |>
    mutate(split_id = formatC(x = n_splits+1, width = 3, flag = 0) ) |>
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
    mutate(method = "pval_holm") |>
    pivot_longer(cols = onset:offset) |>
    select(split_id, method, cluster_id, onset_offset = name, value) |>
    data.frame()

# running the MNE cluster-based permutation
p_values_cluster_mass <- freq_stats_cluster_matrix(
    X = data_matrix - chance_level,
    timesteps = unique(summary_df$time),
    cluster_type = "mass",
    alpha_level = alpha_level
    )
p_values_tfce <- freq_stats_cluster_matrix(
    X = data_matrix - chance_level,
    timesteps = unique(summary_df$time),
    cluster_type = "tfce",
    alpha_level = alpha_level
    )

onset_offset_cluster_mass <- find_clusters(
    data = p_values_cluster_mass |> mutate(pval = pval*(-1) ) |> select(time, value = pval),
    threshold = -alpha_level
    ) |>
    mutate(split_id = formatC(x = n_splits+1, width = 3, flag = 0) ) |>
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
    mutate(method = "cluster_mass") |>
    pivot_longer(cols = onset:offset) |>
    select(split_id, method, cluster_id, onset_offset = name, value) |>
    data.frame()

onset_offset_cluster_tfce <- find_clusters(
    data = p_values_tfce |> mutate(pval = pval*(-1) ) |> select(time, value = pval),
    threshold = -alpha_level
    ) |>
    mutate(split_id = formatC(x = n_splits+1, width = 3, flag = 0) ) |>
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
    mutate(method = "cluster_tfce") |>
    pivot_longer(cols = onset:offset) |>
    select(split_id, method, cluster_id, onset_offset = name, value) |>
    data.frame()

# binding all results together
temp_reliability_results <- bind_rows(
    onset_offset_brms, onset_offset_raw_p, onset_offset_pval_bh,
    onset_offset_pval_by, onset_offset_pval_holm, onset_offset_cpt,
    onset_offset_cluster_mass, onset_offset_cluster_tfce
    )

# and appending it to previous splits results
reliability_results <- bind_rows(reliability_results, temp_reliability_results)

# saving the results
saveRDS(object = reliability_results, file = "./results/reliability_results_1000splits.rds")

# stopping the cluster
stopCluster(cl)

# timing the simulations
toc()
