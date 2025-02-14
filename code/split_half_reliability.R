###########################################################
# Monte-Carlo simulation of onset/offset error properties #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on February 13, 2025                       #
###########################################################

# importing R packages
library(changepoint)
library(reticulate)
library(tidyverse)
library(tidybayes)
library(brms)

# importing Python modules
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
# nsims <- 1e4
nsims <- 20

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
        
        # defining the (main) cluster onsef and offset
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
