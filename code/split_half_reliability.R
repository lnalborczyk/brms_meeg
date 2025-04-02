###########################################################
# Split-half reliability of onset/offset estimates        #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on April 2, 2025                           #
###########################################################

# importing R packages
library(changepoint)
library(reticulate)
library(easystats)
library(tidyverse)
library(tidybayes)
library(scales)
library(brms)

# defining the default ggplot2 theme
theme_set(theme_light(base_size = 12, base_family = "Open Sans") )

# importing Python modules
use_condaenv("r-reticulate3", conda = "~/miniforge3/bin/conda", required = TRUE)
np <- import("numpy")
mne <- import("mne")

# reshaping the decoding results
# read.csv(file = "code/decoding_results.csv") %>%
#     filter(crossval_strategy == "4") %>%
#     select(participant, cv_fold = fold, time, auc = score) %>%
#     mutate(cv_fold = cv_fold + 1) %>%
#     mutate(participant = str_trunc(string = participant, width = 12, ellipsis = "") ) %>%
#     mutate(participant = str_replace(string = participant, pattern = "tickertape", replacement = "participant") ) %>%
#     data.frame() %>%
#     filter(time >= -0.2) %>%
#     write.csv(file = "code/decoding_results_reshaped.csv", row.names = FALSE)

# importing the reshaped MEG decoding results
raw_df <- read.csv(file = "code/decoding_results_reshaped.csv") %>%
    # removing the last participant (to get a round number of participants)
    mutate(participant_id = cur_group_id(), .by = participant) %>%
    filter(participant_id < 33) %>%
    select(-participant_id)

# averaging across participants
summary_df <- raw_df %>%
    summarise(
        auc_mean = mean(auc),
        auc_sd = sd(auc),
        auc_se = sd(auc) / sqrt(n() ),
        .by = c(participant, time)
        )

# plotting the group-level average decoding performance
# raw_df %>%
#     summarise(
#         auc_mean = mean(auc),
#         auc_sd = sd(auc),
#         auc_se = sd(auc) / sqrt(n() ),
#         .by = time
#         ) %>%
#     ggplot(aes(x = time, y = auc_mean) ) +
#     geom_hline(yintercept = 0.5, linetype = 2) +
#     geom_vline(xintercept = 0.0, linetype = 2) +
#     geom_ribbon(
#         aes(ymin = auc_mean - auc_se, ymax = auc_mean + auc_se),
#         alpha = 0.2,
#         ) +
#     geom_line() +
#     labs(x = "Time (s)", y = "Average decoding performance (AUC)")

# saving the plot
# ggsave(
#     filename = "figures/meg_decoding_results.png",
#     width = 12, height = 6, dpi = 300,
#     device = "png"
#     )

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
n_splits <- 100

# creating a dataframe where each participant appears in every split
split_df <- crossing(participant = unique_participants, split_number = 1:n_splits)

# assigning alternating groups (A or B) within each split for each participant
split_df <- split_df %>%
    group_by(split_number) %>%
    # shuffling groups for each split
    mutate(group = sample(rep(c("A", "B"), length.out = n() ) ) ) %>%
    ungroup()

# defining the decoding chance level
chance_level <- 0.5

# defining the threshold on the posterior probability ratio
post_prob_ratio_threshold <- 20

# defining the significance level
alpha_level <- 0.05

# initialising empty results
reliability_results <- data.frame()

# defining the function in R (it will be executed in Python)
freq_stats_cluster_matrix <- function (
        X, timesteps,
        cluster_type = c("mass", "tfce"),
        alpha_level = 0.05
        ) {
    
    # converting R matrix to NumPy array
    X_np <- np$array(X, dtype = "float64")
    
    # defining the type of cluster-based method to use
    if (cluster_type == "mass") {
        
        threshold <- py_none()
        
    } else {
        
        threshold <- dict(start = 0, step = 0.2)
        
    }
    
    # running the statistical test
    results <- mne$stats$spatio_temporal_cluster_1samp_test(
        X = X_np,
        out_type = "mask",
        threshold = threshold,
        n_permutations = as.integer(2^12),
        n_jobs = as.integer(4),
        verbose = TRUE
        )
    
    # extracting results
    T_obs_ <- results[[1]]
    clusters <- results[[2]]
    p_values <- results[[3]]
    
    # retrieving significant clusters
    if (cluster_type == "mass") {
        
        # using the first slice of the 3D array
        p_values_ <- np$transpose(np$ones_like(X_np[1, , drop = FALSE]) )
        
        # defining the (main) cluster onset and offset
        main_cluster_onset <- clusters[which(p_values < alpha_level)][[1]][[1]]$start+1
        main_cluster_offset <- clusters[which(p_values < alpha_level)][[1]][[1]]$stop
        
        # assigning p-values to clusters
        for (i in seq_along(clusters) ) {
            
            # retrieving the current cluster
            cluster_mask <- clusters[[i]][[1]]
            idx_start <- cluster_mask$start+1
            # idx_stop <- cluster_mask$stop
            idx_stop <- cluster_mask$stop+1
            
            # retrieving the p-value for this cluster
            pval <- p_values[[i]]
            
            # assigning this p-value to timesteps belonging
            # to the current cluster
            p_values_[idx_start:idx_stop] <- pval
            
        }
        
        # converting back to a matrix
        results_matrix <- as.matrix(np$squeeze(np$transpose(p_values_) ) )
        
    } else {
        
        # for TFCE
        main_cluster_onset <- min(which(p_values < alpha_level) )
        main_cluster_offset <- max(which(p_values < alpha_level) )
        
        # reshaping the results in a matrix
        results_matrix <- as.matrix(np$squeeze(np$transpose(p_values) ) )
        
    }
    
    # converting back to a dataframe if needed
    results_df <- data.frame(
        pval = results_matrix,
        time = timesteps
        ) %>%
        mutate(
            cluster_onset = time[main_cluster_onset],
            cluster_offset = time[main_cluster_offset]
            )
    
    # returning the results
    return (results_df)
    
}

# computing the p-map
# see https://github.com/easystats/bayestestR/blob/v0.15.2/R/p_map.R
# and https://github.com/easystats/bayestestR/blob/v0.15.2/R/p_map.R#L406
# bayestestR::p_map(x = gam_split, use_iterations = TRUE)
# p_map_custom <- function (post_samples, map_value, null_value = 0.5) {
#     
#     # computing density at null
#     density_null <- density_at(
#         posterior = post_samples, x = null_value, method = "kernel"
#         )
#     
#     if (is.na(density_null) ) density_null <- 0
#     
#     # computing density at MAP
#     density_map <- density_at(
#         posterior = post_samples, x = map_value, method = "kernel"
#         )
#     
#     # computing the odds (ratio)
#     odds <- density_map / density_null
#     
#     # returning these odds
#     return (odds)
#     
# }

# for each simulation
for (i in seq(from = 1, to = n_splits, by = 1) ) {
    
    # printing progress
    cat("\nAnalysing split number:", i, "out of", n_splits, "splits...\n")
    
    # now for each split, assess the onset and offset
    participants_split <- split_df %>%
        filter(split_number == i & group == "A") %>%
        pull(participant)
    
    # sanity check
    cat("\nCurrent split includes participants:", participants_split, "\n")
    
    # splitting data
    raw_df_split <- raw_df %>% filter(participant %in% participants_split)
    summary_df_split <- summary_df %>% filter(participant %in% participants_split)
    
    # default value for k
    # 50 looks good but it depends on the data preprocessing...
    # a sensible default seems to be half the number of timesteps
    # basis_dimension <- as.integer(round(n_distinct(raw_df_split$time) / 2) )
    
    # computing the smallest effect size of interest (SESOI)
    # as the SD of decoding performance during the baseline
    sesoi <- summary_df_split %>%
        filter(time < 0) %>%
        summarise(auc_mean = mean(auc_mean), .by = time) %>%
        summarise(sesoi = sd(auc_mean) ) %>%
        pull(sesoi)
    
    # computing the smallest effect size of interest (SESOI)
    # as the maximum absolute value of decoding performance during the baseline
    # sesoi <- summary_df_split %>%
    #     filter(time < 0) %>%
    #     summarise(auc_mean = mean(auc_mean), .by = time) %>%
    #     summarise(sesoi = max(abs(range(auc_mean) - chance_level) ) ) %>%
    #     pull(sesoi)
    
    # fitting the GAM
    gam_split <- brm(
        # auc_mean | se(auc_sd) ~ s(time, bs = "cr", k = 50) + (1 | participant),
        # auc_mean | se(auc_sd) ~ s(time, bs = "tp", k = 50) + (1 | participant),
        # data = summary_df_split,
        # auc_mean | se(auc_sd) ~ s(time) + s(time, participant, bs = "fs", m = 1),
        # auc_mean ~ s(time, bs = "cr", k = 100),
        # auc_mean ~ s(time, bs = "tp", k = basis_dimension),
        # auc_mean ~ s(time, bs = "cr", k = 50),
        # bf(auc_mean ~ s(time, bs = "cr", k = 50), phi ~ s(time, bs = "cr") ),
        # bf(auc_mean ~ s(time, bs = "cr", k = 50), phi ~ time),
        # auc_mean ~ s(time, bs = "cr", k = 50) + (1 + time | participant),
        auc_mean ~ s(time, bs = "cr", k = 50),
        data = summary_df_split,
        family = Beta(),
        warmup = 2000,
        iter = 5000,
        chains = 8,
        cores = 8
        # backend = "cmdstanr",
        # stan_model_args = list(stanc_options = list("O1") )
        )
    
    # fitting the multilevel GAM (much better PPCs, but very slow)...
    # gamm_split <- brm(
    #     # auc ~ s(time, bs = "cr", k = 20) +
    #     #     s(time, participant, bs = "fs", xt = "cr", k = 20, m = 1),
    #     # data = raw_df_split,
    #     auc_mean | se(auc_sd) ~ s(time, bs = "cr", k = 30) +
    #         s(time, participant, bs = "fs", xt = "cr", k = 30, m = 1),
    #     data = summary_df_split,
    #     family = Beta(),
    #     warmup = 2000,
    #     iter = 5000,
    #     chains = 8,
    #     cores = 8
    #     )
    
    # checking model's predictions against raw data
    plot(
        conditional_effects(x = gam_split),
        line_args = list(colour = "steelblue", fill = "steelblue", alpha = 0.3),
        points = FALSE, plot = FALSE
        )[[1]] +
        geom_hline(yintercept = 0.5, linetype = 2) +
        geom_vline(xintercept = 0.0, linetype = 2) +
        annotate(
            geom = "rect",
            xmin = -Inf, xmax = Inf,
            ymin = chance_level - sesoi,
            ymax = chance_level + sesoi,
            fill ="orangered",
            alpha = 0.2
            ) +
        geom_line(
            data = summary_df_split %>% summarise(auc_mean = mean(auc_mean), .by = time),
            aes(x = time, y = auc_mean), inherit.aes = FALSE
            ) +
        labs(x = "Time (s)", y = "Decoding accuracy (ROC AUC)")
    
    # saving the plot
    ggsave(
        filename = paste0(
            "figures/split_reliability/gam_split_predictions_split_",
            formatC(x = i, width = 3, flag = 0), ".png"
            ),
        width = 12, height = 6, dpi = 300,
        device = "png"
        )
    
    # checking the residuals
    # hist(residuals(gam_split) )
    # hist(residuals(gamm_split) )
    
    # checking the ACF in residuals
    # see https://jroy042.github.io/nonlinear/week4.html
    # and https://ge-chunyu.github.io/posts/2024-04-gamm/
    # and https://jacolienvanrij.com/Tutorials/GAMM.html
    # acf(resid(gam_split), lag.max = 10, main = "ACF")
    # acf(resid(gamm_split), lag.max = 10, main = "ACF")
    
    # retrieving the auto-correlation rho
    # from https://jacolienvanrij.com/Tutorials/GAMM.html
    # (valRho <- acf(resid(gamm_split), plot = FALSE)$acf[2])
    
    # posterior predictive checking
    # pp_check(gam_split, prefix = "ppd")
    # pp_check(gam_split)
    # pp_check(gam_split, type = "error_hist")
    # pp_check(gam_split, type = "stat_2d")
    # pp_check(gamm_split, type = "ecdf_overlay")
    
    # initialising dataframe to store brms results
    # temp_brms_results <- data.frame(
    #     split_id = formatC(x = i, width = 3, flag = 0),
    #     onset_brms = NA, offset_brms = NA
    #     )
    
    # computing the posterior probability
    prob_y_above_0 <- data.frame(time = unique(gam_split$data$time) ) %>%
        # gam_split$data %>%
        # using a subset of posterior samples
        # add_epred_draws(object = gam_split, ndraws = n_post_samples) %>%
        # or using all posterior samples
        add_epred_draws(object = gam_split) %>%
        # converting to dataframe
        data.frame() %>%
        # dplyr::select(participant, time, .epred, .draw) %>%
        # dplyr::select(time, .epred, .draw) %>%
        # pivot_wider(names_from = condition, values_from = .epred) %>%
        # mutate(epred_diff = cond2 - cond1) %>%
        # computing mean posterior probability at the participant level
        # group_by(participant, time) %>%
        # summarise(m = mean(.epred > (0 + chance_level) ) ) %>%
        # ungroup() %>%
        # computing mean posterior probability at the group level
        group_by(time) %>%
        summarise(
            # m = mean(.epred > (0 + chance_level) )
            m = mean(.epred > (0 + chance_level + sesoi) ),
            m_full = mean(.epred > (0 + chance_level + full_sesoi) )
            # map_est = map_estimate(x = .epred, method = "kernel")$MAP_Estimate,
            # d_at_map = density_at(posterior = .epred, x = map_est, method = "kernel"),
            # d_at_null = density_at(posterior = .epred, x = chance_level, method = "kernel"),
            # p_map = p_map_custom(post_samples = .epred, map_value = map_est, null_value = chance_level)
            ) %>%
        mutate(
            prob_ratio = m / (1 - m),
            prob_ratio_full = m_full / (1 - m_full)
            ) %>%
        # ensuring there is no 0 or +Inf values
        mutate(prob_ratio = ifelse(is.infinite(prob_ratio), ndraws(gam_split), prob_ratio) ) %>%
        mutate(prob_ratio = ifelse(prob_ratio == 0, 1 / ndraws(gam_split), prob_ratio) ) %>%
        mutate(prob_ratio_full = ifelse(is.infinite(prob_ratio_full), ndraws(gam_split), prob_ratio_full) ) %>%
        mutate(prob_ratio_full = ifelse(prob_ratio_full == 0, 1 / ndraws(gam_split), prob_ratio_full) ) %>%
        # mutate(d_at_null = ifelse(is.na(d_at_null), 0, d_at_null) ) %>%
        # mutate(p_map = d_at_map / d_at_null) %>%
        ungroup()
    
    # sanity visual checks
    prob_y_above_0 %>%
        ggplot(aes(x = time, y = prob_ratio) ) +
        geom_hline(yintercept = 1, linetype = "dashed") +
        geom_hline(yintercept = 3, linetype = "dashed", color = "darkred") +
        geom_hline(yintercept = 1/3, linetype = "dashed", color = "darkred") +
        geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
        geom_hline(yintercept = post_prob_ratio_threshold, linetype = "dashed", color = "darkgreen") +
        geom_hline(yintercept = 1/10, linetype = "dashed", color = "red") +
        geom_hline(yintercept = 100, linetype = "dashed", color = "orangered") +
        geom_hline(yintercept = 1/100, linetype = "dashed", color = "orangered") +
        geom_vline(xintercept = 0, linetype = 2) +
        geom_line() +
        scale_y_log10(labels = label_log(digits = 2) ) +
        labs(x = "Time (s)", y = "Posterior probability ratio")
    
    # saving the plot
    ggsave(
        filename = paste0(
            "figures/split_reliability/gam_split_post_prob_ratio_split_",
            formatC(x = i, width = 3, flag = 0), ".png"
            ),
        width = 12, height = 6, dpi = 300,
        device = "png"
        )
    
    # finding clusters with neurogam
    onset_offset_brms <- neurogam::find_clusters(
        data = prob_y_above_0 %>% select(time, value = prob_ratio),
        threshold = post_prob_ratio_threshold
        ) %>%
        mutate(split_id = formatC(x = i, width = 3, flag = 0) ) %>%
        select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
        mutate(method = "brms") %>%
        pivot_longer(cols = onset:offset) %>%
        select(split_id, method, cluster_id, onset_offset = name, value) %>%
        data.frame()
    
    # converting data to a matrix (for later use in MNE functions)
    data_matrix <- matrix(
        summary_df_split$auc_mean,
        ncol = length(unique(summary_df_split$time) ),
        byrow = TRUE
        )
    
    # massive univariate t-tests
    tests_results <- summary_df_split %>%
        group_by(time) %>%
        summarise(
            tval = t.test(x = auc_mean, mu = chance_level)$statistic^2,
            pval = t.test(x = auc_mean, mu = chance_level)$p.value
            ) %>%
        mutate(
            pval_bh = p.adjust(p = pval, method = "BH"),
            pval_by = p.adjust(p = pval, method = "BY"),
            pval_holm = p.adjust(p = pval, method = "holm")
            ) %>%
        ungroup()
        
    # using the binary segmentation method to identify onsets and offsets
    # res <- cpt.meanvar(data = tests_results$tval, method = "BinSeg", Q = 10)
    res <- cpt.meanvar(data = tests_results$tval, method = "BinSeg", Q = 2)
    onset_offset_cpt <- data.frame(
        split_id = formatC(x = i, width = 3, flag = 0),
        cluster_id = 1,
        onset = tests_results$time[res@cpts[1]],
        offset = tests_results$time[res@cpts[2]]
        ) %>%
        mutate(method = "cpt") %>%
        pivot_longer(cols = onset:offset) %>%
        select(split_id, method, cluster_id, onset_offset = name, value) %>%
        data.frame()
    
    # finding clusters
    onset_offset_raw_p <- neurogam::find_clusters(
        data = tests_results %>% mutate(pval = pval * (-1) ) %>% select(time, value = pval),
        threshold = -alpha_level
        ) %>%
        mutate(split_id = formatC(x = i, width = 3, flag = 0) ) %>%
        select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
        mutate(method = "raw_p") %>%
        pivot_longer(cols = onset:offset) %>%
        select(split_id, method, cluster_id, onset_offset = name, value) %>%
        data.frame()
    
    onset_offset_pval_bh <- neurogam::find_clusters(
        data = tests_results %>% mutate(pval = pval_bh * (-1) ) %>% select(time, value = pval),
        threshold = -alpha_level
        ) %>%
        mutate(split_id = formatC(x = i, width = 3, flag = 0) ) %>%
        select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
        mutate(method = "pval_bh") %>%
        pivot_longer(cols = onset:offset) %>%
        select(split_id, method, cluster_id, onset_offset = name, value) %>%
        data.frame()
    
    onset_offset_pval_by <- neurogam::find_clusters(
        data = tests_results %>% mutate(pval = pval_by * (-1) ) %>% select(time, value = pval),
        threshold = -alpha_level
        ) %>%
        mutate(split_id = formatC(x = i, width = 3, flag = 0) ) %>%
        select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
        mutate(method = "pval_by") %>%
        pivot_longer(cols = onset:offset) %>%
        select(split_id, method, cluster_id, onset_offset = name, value) %>%
        data.frame()
    
    onset_offset_pval_holm <- neurogam::find_clusters(
        data = tests_results %>% mutate(pval = pval_holm * (-1) ) %>% select(time, value = pval),
        threshold = -alpha_level
        ) %>%
        mutate(split_id = formatC(x = i, width = 3, flag = 0) ) %>%
        select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
        mutate(method = "pval_holm") %>%
        pivot_longer(cols = onset:offset) %>%
        select(split_id, method, cluster_id, onset_offset = name, value) %>%
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
    
    onset_offset_cluster_mass <- neurogam::find_clusters(
        data = p_values_cluster_mass %>% mutate(pval = pval*(-1) ) %>% select(time, value = pval),
        threshold = -alpha_level
        ) %>%
        mutate(split_id = formatC(x = i, width = 3, flag = 0) ) %>%
        select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
        mutate(method = "cluster_mass") %>%
        pivot_longer(cols = onset:offset) %>%
        select(split_id, method, cluster_id, onset_offset = name, value) %>%
        data.frame()
    
    onset_offset_cluster_tfce <- neurogam::find_clusters(
        data = p_values_tfce %>% mutate(pval = pval*(-1) ) %>% select(time, value = pval),
        threshold = -alpha_level
        ) %>%
        mutate(split_id = formatC(x = i, width = 3, flag = 0) ) %>%
        select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
        mutate(method = "cluster_tfce") %>%
        pivot_longer(cols = onset:offset) %>%
        select(split_id, method, cluster_id, onset_offset = name, value) %>%
        data.frame()
    
    # binding all results together
    temp_reliability_results <- bind_rows(
        onset_offset_brms, onset_offset_raw_p, onset_offset_pval_bh,
        onset_offset_pval_by, onset_offset_pval_holm, onset_offset_cpt,
        onset_offset_cluster_mass, onset_offset_cluster_tfce
        )
    
    # and appending it to previous splits results
    reliability_results <- bind_rows(reliability_results, temp_reliability_results)
    
}

# fitting the final GAM
full_gam <- brm(
    auc_mean ~ s(time, bs = "cr", k = 50),
    data = summary_df,
    family = Beta(),
    warmup = 2000,
    iter = 5000,
    chains = 8,
    cores = 8
    )

# computing the posterior probability
prob_y_above_0 <- data.frame(time = unique(full_gam$data$time) ) %>%
    # retrieving the posterior samples
    add_epred_draws(object = full_gam) %>%
    # converting to dataframe
    data.frame() %>%
    # computing mean posterior probability at the group level
    group_by(time) %>%
    summarise(m = mean(.epred > (0 + chance_level + full_sesoi) ) ) %>%
    mutate(prob_ratio = m / (1 - m) ) %>%
    # ensuring there is no 0 or +Inf values
    mutate(prob_ratio = ifelse(is.infinite(prob_ratio), ndraws(full_gam), prob_ratio) ) %>%
    mutate(prob_ratio = ifelse(prob_ratio == 0, 1 / ndraws(full_gam), prob_ratio) ) %>%
    ungroup()

# finding clusters with neurogam
onset_offset_brms <- neurogam::find_clusters(
    data = prob_y_above_0 %>% select(time, value = prob_ratio),
    threshold = post_prob_ratio_threshold
    ) %>%
    mutate(split_id = formatC(x = i+1, width = 3, flag = 0) ) %>%
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
    mutate(method = "brms") %>%
    pivot_longer(cols = onset:offset) %>%
    select(split_id, method, cluster_id, onset_offset = name, value) %>%
    data.frame()

# converting data to a matrix (for later use in MNE functions)
data_matrix <- matrix(
    summary_df$auc_mean,
    ncol = length(unique(summary_df$time) ),
    byrow = TRUE
    )

# massive univariate t-tests
tests_results <- summary_df %>%
    group_by(time) %>%
    summarise(
        tval = t.test(x = auc_mean, mu = chance_level)$statistic^2,
        pval = t.test(x = auc_mean, mu = chance_level)$p.value
        ) %>%
    mutate(
        pval_bh = p.adjust(p = pval, method = "BH"),
        pval_by = p.adjust(p = pval, method = "BY"),
        pval_holm = p.adjust(p = pval, method = "holm")
        ) %>%
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
    split_id = formatC(x = i+1, width = 3, flag = 0),
    cluster_id = 1,
    onset = tests_results$time[res@cpts[1]],
    offset = tests_results$time[res@cpts[2]]
    ) %>%
    mutate(method = "cpt") %>%
    pivot_longer(cols = onset:offset) %>%
    select(split_id, method, cluster_id, onset_offset = name, value) %>%
    data.frame()

# finding clusters
onset_offset_raw_p <- neurogam::find_clusters(
    data = tests_results %>% mutate(pval = pval * (-1) ) %>% select(time, value = pval),
    threshold = -alpha_level
    ) %>%
    mutate(split_id = formatC(x = i+1, width = 3, flag = 0) ) %>%
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
    mutate(method = "raw_p") %>%
    pivot_longer(cols = onset:offset) %>%
    select(split_id, method, cluster_id, onset_offset = name, value) %>%
    data.frame()

onset_offset_pval_bh <- neurogam::find_clusters(
    data = tests_results %>% mutate(pval = pval_bh * (-1) ) %>% select(time, value = pval),
    threshold = -alpha_level
    ) %>%
    mutate(split_id = formatC(x = i+1, width = 3, flag = 0) ) %>%
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
    mutate(method = "pval_bh") %>%
    pivot_longer(cols = onset:offset) %>%
    select(split_id, method, cluster_id, onset_offset = name, value) %>%
    data.frame()

onset_offset_pval_by <- neurogam::find_clusters(
    data = tests_results %>% mutate(pval = pval_by * (-1) ) %>% select(time, value = pval),
    threshold = -alpha_level
    ) %>%
    mutate(split_id = formatC(x = i+1, width = 3, flag = 0) ) %>%
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
    mutate(method = "pval_by") %>%
    pivot_longer(cols = onset:offset) %>%
    select(split_id, method, cluster_id, onset_offset = name, value) %>%
    data.frame()

onset_offset_pval_holm <- neurogam::find_clusters(
    data = tests_results %>% mutate(pval = pval_holm * (-1) ) %>% select(time, value = pval),
    threshold = -alpha_level
    ) %>%
    mutate(split_id = formatC(x = i+1, width = 3, flag = 0) ) %>%
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
    mutate(method = "pval_holm") %>%
    pivot_longer(cols = onset:offset) %>%
    select(split_id, method, cluster_id, onset_offset = name, value) %>%
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

onset_offset_cluster_mass <- neurogam::find_clusters(
    data = p_values_cluster_mass %>% mutate(pval = pval*(-1) ) %>% select(time, value = pval),
    threshold = -alpha_level
    ) %>%
    mutate(split_id = formatC(x = i+1, width = 3, flag = 0) ) %>%
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
    mutate(method = "cluster_mass") %>%
    pivot_longer(cols = onset:offset) %>%
    select(split_id, method, cluster_id, onset_offset = name, value) %>%
    data.frame()

onset_offset_cluster_tfce <- neurogam::find_clusters(
    data = p_values_tfce %>% mutate(pval = pval*(-1) ) %>% select(time, value = pval),
    threshold = -alpha_level
    ) %>%
    mutate(split_id = formatC(x = i+1, width = 3, flag = 0) ) %>%
    select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) %>%
    mutate(method = "cluster_tfce") %>%
    pivot_longer(cols = onset:offset) %>%
    select(split_id, method, cluster_id, onset_offset = name, value) %>%
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
saveRDS(object = reliability_results, file = "results/reliability_results.rds")
