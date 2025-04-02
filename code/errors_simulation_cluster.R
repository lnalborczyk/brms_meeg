###########################################################
# Monte-Carlo simulation of onset/offset error properties #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on March 31, 2025                          #
###########################################################

# importing R packages
library(changepoint)
library(reticulate)
library(tidyverse)
library(tidybayes)
library(brms)

# for setting up the cluster
library(doParallel)
library(foreach)

# total number of total CPU cores available on the HPC node
total_cores <- 64

# number of cores per brms model
cores_per_model <- 8

# number of parallel simulations to run at once
n_parallel <- total_cores / cores_per_model

# registering parallel backend
cl <- makeCluster(n_parallel)
registerDoParallel(cl)

# importing Python modules
# use_condaenv("/Users/ladislas/opt/anaconda3/envs/r-reticulate/bin/python", required = TRUE)
# use_condaenv("r-reticulate3", conda = "~/miniforge3/bin/conda", required = TRUE)
use_condaenv("r-reticulate", required = TRUE)
np <- import("numpy")
mne <- import("mne")

# importing R version of Matlab code from Yeung et al. (2004)
source("code/eeg_noise.R")

# importing the ERP template with true onset = 160 ms, F=81, and max at F=126
source("code/erp_template.R")

# importing helper functions
source("code/functions.R")

# to use with the eeg_noise function
meanpower <- unlist(read.table("code/meanpower.txt") )

# defining the number of simulations to be performed
nsims <- 1e4

# defining simulation parameters
n_trials <- 50 # number of trials
n_ppt <- 20 # number of participants
outvar <- 1 # noise variance
srate <- 500 # sampling rate in Hz
ronset <- seq(from = 150, to = 170, by = 2) # random onset for each participant

# defining the true onset and offset in seconds
true_onset <- 0.160
true_offset <- 0.342

# defining significance (alpha) thresholds
significance_thresholds <- c(0.05, 0.01, 0.005, 0.001)

# defining a function to generate a dataframe
generate_data <- function (n_trials, n_ppt, outvar, srate, ronset) {
    
    for (P in 1:n_ppt) { # for each participant
        
        # get random onset
        ponset <- sample(x = ronset, size = 1)
        
        # find starting point
        st <- which(Xf == ponset)
        
        # pad vector
        temp2 <- c(rep(0, st - 2), erp, rep(0, Nf - st - length(erp) + 2) )
        
        # initialising empty conditions
        cond1 <- matrix(0, nrow = n_trials, ncol = Nf)
        cond2 <- matrix(0, nrow = n_trials, ncol = Nf)
        
        for (trial in 1:n_trials) { # for each trial
            
            cond1[trial, ] <- temp1 + eeg_noise(frames = Nf, srate = srate, outvar = outvar, meanpower)
            cond2[trial, ] <- temp2 + eeg_noise(frames = Nf, srate = srate, outvar = outvar, meanpower)
            
        }
        
        # converting results to dataframe
        temp_df <- data.frame(
            x = rep(Xf, 2 * nrow(cond1) ),
            y = c(c(t(cond1) ), c(t(cond2) ) ),
            trial = c(rep(1:n_trials, each = length(Xf) ), rep(1:n_trials, each = length(Xf) ) ),
            condition = factor(rep(c("cond1", "cond2"), each = Nf * n_trials) ),
            participant = paste0("participant_", sprintf("%02d", P) )
        ) %>%
            select(participant, condition, trial, x, y)
        
        # and appending it to previous participants
        if (exists("raw_data") ) {
            
            raw_data <- bind_rows(raw_data, temp_df)
            
        } else {
            
            raw_data <- temp_df
            
        }
        
    }
    
    # converting time from ms to seconds
    raw_data <- raw_data %>%
        # converting time from ms to seconds
        mutate(x = x / 1000) %>%
        # renaming columns
        rename(time = x, eeg = y)
    
    # returning the simulated data
    return (raw_data)
    
}

# defining a function to find onset and offset in timeseries
find_onset_offset <- function (mask, timeseries) {
    
    # initialising onset and offset values
    onset <- NA
    offset <- NA
    
    # identifying the onset and offset
    masked_timesteps <- try(timeseries[which(mask)], silent = TRUE)
    onset <- head(x = masked_timesteps, n = 1)
    offset <- tail(x = masked_timesteps, n = 1)
    
    # returning the onset and offset
    return (c(onset, offset) )
    
}

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
            idx_stop <- cluster_mask$stop
            
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

# generate some data
# raw_df <- generate_data(
#     n_trials = n_trials, n_ppt = n_ppt, outvar = outvar,
#     srate = srate, ronset = ronset
#     )

# testing the find_onset_offset() function
# find_onset_offset(
#     mask = tests_results$pval[1:nrow(tests_results)] <= aath,
#     timeseries = tests_results$time[1:nrow(tests_results)]
#     )

# initialising empty simulation results
sim_results <- data.frame()

# for each simulation
# for (i in seq(from = 1, to = nsims, by = 1) ) {
sim_results_list <- foreach(i = 1:nsims, .packages = c("brms", "tidyverse", "data.table", "mne", "reticulate") ) %dopar% {
    
    # printing progress
    cat("\n-------- Simulation number:", i, "out of", nsims, "simulations...\n")
    
    # simulating some EEG data
    raw_df <- generate_data(
        n_trials = n_trials, n_ppt = n_ppt, outvar = outvar,
        srate = srate, ronset = ronset
        )
    
    # summarising raw_data per participant
    ppt_data <- raw_df %>%
        summarise(eeg = mean(eeg), .by = c(participant, condition, time) ) %>%
        pivot_wider(names_from = condition, values_from = eeg) %>%
        mutate(eeg_diff = cond2 - cond1)
    
    # summarising raw_data per participant (mean + SD)
    summary_df <- raw_df %>%
        summarise(
            eeg_mean = mean(eeg),
            eeg_sd = sd(eeg),
            .by = c(participant, condition, time)
            )
    
    # defining a contrast for condition
    contrasts(summary_df$condition) <- c(-0.5, 0.5)
    
    # fitting the GAMM
    meta_gam <- brm(
        # using by-participant SD of ERPs across trials
        eeg_mean | se(eeg_sd) ~
            condition + s(time, bs = "cr", k = 20, by = condition) +
            (1 | participant),
        data = summary_df,
        family = gaussian(),
        warmup = 2000,
        iter = 5000,
        chains = 8,
        cores = cores_per_model
        # backend = "cmdstanr",
        # opencl = opencl(c(0, 0) )
        # threads = threading(1),
        # stan_model_args = list(stanc_options = list("O1") )
        )
    
    # checking the model's predictions
    # conditional_effects(meta_gam)
    
    # defining values of eps and threshold to test
    # eps_values <- seq(from = 0, to = 0.1, by = 0.01)
    # threshold_values <- seq(from = 1, to = 30, by = 1)
    eps_values <- 0
    threshold_values <- c(1, 3, 6, 10, 20, 50, 100)
    
    # number of posterior samples to use
    n_post_samples <- 1e3
    
    # initialising eps_threshold_results dataframe
    eps_threshold_results <- crossing(
        eps = eps_values, threshold = threshold_values
        ) %>%
        mutate(
            onset_brms = NA, offset_brms = NA,
            simulation_id = formatC(x = i, width = 3, flag = 0)
            )
    
    # looping over different values of eps and threshold
    for (j in seq_len(nrow(eps_threshold_results) ) ) {
        
        # printing progress
        # cat("Assessing combination:", j, "out of", nrow(eps_threshold_results), "combinations...\n")
        
        # retrieving current eps and threshold values
        eps <- eps_threshold_results$eps[j]
        threshold <- eps_threshold_results$threshold[j]
        
        # computing probability metrics with the current eps value
        prob_y_above_0 <- meta_gam$data %>%
            # crossing(
            # time = unique(meta_gam$data$time),
            # condition = unique(meta_gam$data$condition),
            # participant = unique(meta_gam$data$participant),
            # ) %>%
            # meta_gam$data %>%
            add_epred_draws(object = meta_gam, ndraws = n_post_samples) %>%
            # add_epred_draws(object = meta_gam) %>%
            data.frame() %>%
            dplyr::select(participant, time, condition, .epred, .draw) %>%
            pivot_wider(names_from = condition, values_from = .epred) %>%
            mutate(epred_diff = cond2 - cond1) %>%
            # computing mean posterior probability at the participant level
            # group_by(participant, time) %>%
            # summarise(m = mean(epred_diff > eps) ) %>%
            # ungroup() %>%
            # computing mean posterior prob at the group level
            group_by(time) %>%
            summarise(m = mean(epred_diff > eps) ) %>%
            # summarise(m = mean(m) ) %>%
            mutate(prob_ratio = m / (1 - m) ) %>%
            ungroup()
        
        # sanity visual check
        # prob_y_above_0 %>%
        #     # ggplot(aes(x = time, y = m) ) +
        #     ggplot(aes(x = time, y = log(prob_ratio) ) ) +
        #     geom_line()
        
        # finding onset, offset, and peak for the current threshold
        exceeding_times <- prob_y_above_0 %>%
            dplyr::filter(prob_ratio > threshold) %>%
            summarise(
                cluster_onset = min(time, na.rm = TRUE),
                cluster_offset = max(time, na.rm = TRUE)
            ) %>%
            mutate(true_onset = true_onset, true_offset = true_offset)
        
        # storing the results in the dataframe
        if (nrow(exceeding_times) > 0) {
            
            # filling onset+-/offset values
            eps_threshold_results$onset_brms[j] <- exceeding_times$cluster_onset
            eps_threshold_results$offset_brms[j] <- exceeding_times$cluster_offset
            
        }
        
    }
    
    # converting data to a matrix (for later use in MNE functions)
    data_matrix <- matrix(
        ppt_data$eeg_diff,
        ncol = length(unique(ppt_data$time) ),
        byrow = TRUE
    )
    
    # massive univariate t-tests
    tests_results <- ppt_data %>%
        group_by(time) %>%
        summarise(
            tval = t.test(x = eeg_diff, mu = 0)$statistic^2,
            pval = t.test(x = eeg_diff, mu = 0)$p.value
            ) %>%
        mutate(
            pval_bh = p.adjust(p = pval, method = "BH"),
            pval_by = p.adjust(p = pval, method = "BY"),
            pval_holm = p.adjust(p = pval, method = "holm")
            ) %>%
        ungroup()
    
    # initialising empty dataframe
    temp_sim_results <- data.frame()
    
    # for each significance (alpha) threshold
    for (alpha_level in significance_thresholds) {
        
        # using the binary segmentation method to identify onsets and offsets
        res <- cpt.meanvar(data = tests_results$tval, method = "BinSeg", Q = 2)
        onset_offset_cpt <- c(tests_results$time[res@cpts[1]], tests_results$time[res@cpts[2]])
        
        # finding onsets
        onset_offset_raw_p <- find_onset_offset(
            tests_results$pval <= alpha_level,
            tests_results$time
            )
        onset_offset_bh <- find_onset_offset(
            tests_results$pval_bh <= alpha_level,
            tests_results$time
            )
        onset_offset_by <- find_onset_offset(
            tests_results$pval_by <= alpha_level,
            tests_results$time
            )
        onset_offset_holm <- find_onset_offset(
            tests_results$pval_holm <= alpha_level,
            tests_results$time
            )
        
        # running the MNE cluster-based permutation
        p_values_cluster_mass <- freq_stats_cluster_matrix(
            X = data_matrix,
            timesteps = unique(ppt_data$time),
            cluster_type = "mass",
            alpha_level = alpha_level
            )
        p_values_tfce <- freq_stats_cluster_matrix(
            X = data_matrix,
            timesteps = unique(ppt_data$time),
            cluster_type = "tfce",
            alpha_level = alpha_level
            )
        
        # putting everything together
        alpha_temp_sim_results <- data.frame(
            simulation_id = formatC(x = i, width = 3, flag = 0),
            alpha_level = alpha_level,
            onset_p = onset_offset_raw_p[1],
            offset_p = onset_offset_raw_p[2],
            onset_bh = onset_offset_bh[1],
            offset_bh = onset_offset_bh[2],
            onset_by = onset_offset_by[1],
            offset_by = onset_offset_by[2],
            onset_holm = onset_offset_holm[1],
            offset_holm = onset_offset_holm[2],
            onset_cpt = onset_offset_cpt[1],
            offset_cpt = onset_offset_cpt[2],
            onset_cluster_mass = unique(p_values_cluster_mass$cluster_onset),
            offset_cluster_mass = unique(p_values_cluster_mass$cluster_offset),
            onset_cluster_tfce = unique(p_values_tfce$cluster_onset),
            offset_cluster_tfce = unique(p_values_tfce$cluster_offset)
            )
        
        # appending to previous results
        temp_sim_results <- bind_rows(temp_sim_results, alpha_temp_sim_results)
        
    }
    
    # bind these results with the brms results
    eps_threshold_results <- relocate(eps_threshold_results, simulation_id, .before = 1)
    temp_sim_results <- left_join(eps_threshold_results, temp_sim_results, by = "simulation_id")
    
    # and appending it to previous simulation results
    # sim_results <- bind_rows(sim_results, temp_sim_results)
    
    # returning the results
    return (temp_sim_results)
    
}

# binding all simulation results
sim_results <- bind_rows(sim_results_list)

# saving the simulation results
saveRDS(object = sim_results, file = "results/sim_results.rds")

# stopping the cluster
stopCluster(cl)
