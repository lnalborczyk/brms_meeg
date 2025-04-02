###########################################################
# Split-half reliability of onset/offset estimates        #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on March 31, 2025                          #
###########################################################

# loading the reticulate package
library(reticulate)

# importing the numpy and mne python modules
use_condaenv("r-reticulate3", conda = "~/miniforge3/bin/conda", required = TRUE)
np <- import("numpy")
mne <- import("mne")

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

# converting data to a matrix (for later use in MNE functions)
data_matrix <- matrix(
    ppt_data$eeg_diff,
    ncol = length(unique(ppt_data$time) ),
    byrow = TRUE
    )

# running the MNE cluster-based permutation
# alpha_level <- 0.5
p_values_cluster_mass <- freq_stats_cluster_matrix(
    X = data_matrix,
    timesteps = unique(ppt_data$time),
    cluster_type = "mass",
    alpha_level = aath
    )

p_values_tfce <- freq_stats_cluster_matrix(
    X = data_matrix,
    timesteps = unique(ppt_data$time),
    cluster_type = "tfce",
    alpha_level = aath
    )

# retrieving the onset and offset
onset_cluster_mass <- unique(p_values_cluster_mass$cluster_onset)
offset_cluster_mass <- unique(p_values_cluster_mass$cluster_offset)
onset_cluster_mass_tfce <- unique(p_values_tfce$cluster_onset)
offset_cluster_mass_tfce <- unique(p_values_tfce$cluster_offset)
