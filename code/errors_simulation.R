###########################################################
# Monte-Carlo simulation of onset/offset error properties #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on February 10, 2025                       #
###########################################################

library(changepoint)
library(reticulate)
library(tidyverse)
library(brms)

# retrieving the code from the Quarto document
# knitr::purl(input = "brms_meeg.qmd", output = "brms_meeg.R", documentation = 0)

# importing R version of Matlab code from Yeung et al. (2004)
source("eeg_noise.R")

# importing the ERP template with true onset = 160 ms, F=81, and max at F=126
source("erp_template.R")

# importing helper functions
source("functions.R")

# to use with the eeg_noise function
meanpower <- unlist(read.table("meanpower.txt") )

# defining the number of simulations to be performed
nsims <- 10

# defining simulation parameters
n_trials <- 50 # number of trials
n_ppt <- 20 # number of participants
outvar <- 1 # noise variance
srate <- 500 # sampling rate in Hz
ronset <- seq(from = 150, to = 170, by = 2) # random onset for each participant

# defining the true onset and offset in seconds
true_onset <- 0.160
true_offset <- 0.342

# defining an arbitrary alpha threshold
aath <- 0.05

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
        
        for (T in 1:n_trials) { # for each trial
            
            cond1[T, ] <- temp1 + eeg_noise(frames = Nf, srate = srate, outvar = outvar, meanpower)
            cond2[T, ] <- temp2 + eeg_noise(frames = Nf, srate = srate, outvar = outvar, meanpower)
            
        }
        
        # converting results to dataframe
        temp_df <- data.frame(
            x = rep(Xf, 2*nrow(cond1) ),
            y = c(c(t(cond1) ), c(t(cond2) ) ),
            trial = c(rep(1:n_trials, each = length(Xf) ), rep(1:n_trials, each = length(Xf) ) ),
            condition = factor(rep(c("cond1", "cond2"), each = Nf * n_trials) ),
            participant = paste0("participant_", sprintf("%02d", P) )
        ) %>%
            select(participant, condition, trial, x, y)
        
        # and appending it to previous participants
        if (exists("raw_df") ) {
            
            raw_df <- bind_rows(raw_df, temp_df)
            
        } else {
            
            raw_df <- temp_df
            
        }
        
    }
    
    # converting time from ms to seconds
    raw_df <- raw_df %>%
        # converting time from ms to seconds
        mutate(x = x / 1000) %>%
        # renaming columns
        rename(time = x, eeg = y)
    
    # returning the simulated data
    return (raw_df)
    
}

# initialising empty simulation results
sim_results <- data.frame()

# for each simulation
for (i in seq_range(x = nsims) ) {
    
    # simulating some EEG data
    raw_df <- generate_data()
    
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
            condition + s(time, bs = "cr", k = 10, by = condition) +
            (1 | participant),
        data = summary_df,
        family = gaussian(),
        warmup = 2000,
        iter = 5000,
        chains = 8,
        cores = 8
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
    
    # finding onsets
    onset_p <- find_onset(
        tests_results$pval[1:nrow(tests_results)] <= aath,
        tests_results$time[1:nrow(tests_results)]
    )
    onset_bh <- find_onset(
        tests_results$pval_bh[1:nrow(tests_results)] <= aath,
        tests_results$time[1:nrow(tests_results)]
    )
    onset_by <- find_onset(
        tests_results$pval_by[1:nrow(tests_results)] <= aath,
        tests_results$time[1:nrow(tests_results)]
    )
    onset_holm <- find_onset(
        tests_results$pval_holm[1:nrow(tests_results)] <= aath,
        tests_results$time[1:nrow(tests_results)]
    )
    
    # using the changepoint package to identify onsets and offsets
    res <- cpt.meanvar(data = tests_results$tval, method = "BinSeg", Q = 2)
    cpt_onset <- res@cpts[1]
    cpt_offset <- res@cpts[2]
    
    # defining values of eps and threshold to test
    eps_values <- seq(from = 0, to = 0.1, by = 0.01)
    threshold_values <- seq(from = 1, to = 30, by = 1)
    
    # number of posterior samples to use
    n_post_samples <- 1e3
    
    # initialising results dataframe
    results <- crossing(eps = eps_values, threshold = threshold_values) %>%
        mutate(
            estimated_onset = NA, estimated_peak = NA, estimated_offset = NA,
            error_onset = NA, error_peak = NA, error_offset = NA
        )
    
    # looping over different values of eps and threshold
    for (i in seq_len(nrow(results) ) ) {
        
        # printing progress
        cat("Assessing combination:", i, "out of", nrow(results), "combinations...")
        
        # retrieving current eps and threshold values
        eps <- results$eps[i]
        threshold <- results$threshold[i]
        
        # computing probability metrics with the current eps value
        prob_y_above_0 <- meta_gam$data %>%
            add_epred_draws(object = meta_gam, ndraws = n_post_samples) %>%
            data.frame() %>%
            dplyr::select(participant, time, condition, .epred, .draw) %>%
            pivot_wider(names_from = condition, values_from = .epred) %>%
            mutate(epred_diff = cond2 - cond1) %>%
            # computing mean posterior prob at the participant level
            group_by(participant, time) %>%
            summarise(m = mean(epred_diff > eps) ) %>%
            ungroup() %>%
            # computing mean posterior prob at the group level
            group_by(time) %>%
            summarise(m = mean(m) ) %>%
            mutate(prob_ratio = m / (1 - m) ) %>%
            ungroup()
        
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
            
            results$estimated_onset[i] <- exceeding_times$cluster_onset
            results$estimated_offset[i] <- exceeding_times$cluster_offset
            
            # computing errors
            results$error_onset[i] <- abs(
                exceeding_times$cluster_onset - exceeding_times$true_onset
            )
            results$error_offset[i] <- abs(
                exceeding_times$cluster_offset - exceeding_times$true_offset
            )
            
        }
        
    }
    
    
    # importing the numpy and mne python modules
    np <- import("numpy")
    mne <- import("mne")
    
    # defining the function in R (it will be executed in Python)
    freq_stats_gat_matrix <- function (X) {
        
        # converting R matrix to NumPy array
        X_np <- np$array(X, dtype = "float64")
        
        # running the statistical test
        results <- mne$stats$spatio_temporal_cluster_1samp_test(
            X_np,
            out_type = "mask",
            n_permutations = as.integer(1000),
            n_jobs = as.integer(4),
            verbose = TRUE
        )
        
        # extracting results
        T_obs_ <- results[[1]]
        clusters <- results[[2]]
        p_values <- results[[3]]
        
        # using the first slice of the 3D array
        p_values_ <- np$transpose(np$ones_like(X_np[1, , drop = FALSE]) )
        
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
        
        # converting result back to R format
        return (as.matrix(np$squeeze(np$transpose(p_values_) ) ) )
        
    }
    
    # converting decoding_data to a matrix
    data_matrix <- matrix(
        decoding_data$auc,
        ncol = length(unique(decoding_data$time) ),
        byrow = TRUE
        )
    
    # running the function
    chance_level <- 0.5
    p_values_results <- freq_stats_gat_matrix(data_matrix-chance_level)
    
    # converting back to a dataframe if needed
    p_values_df <- data.frame(
        pval = p_values_results,
        time = decoding_data$time
        )
    
}

# saving the simulation results
saveRDS(object = sim_results, file = "results/sim_results.rds")
