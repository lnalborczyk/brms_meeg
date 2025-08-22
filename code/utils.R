###########################################################
# Defining some utility functions                         #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on April 4, 2025                           #
###########################################################

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

# defining a function to identify clusters
# from https://github.com/lnalborczyk/neurogam
find_clusters <- function (data, threshold = 10, above_threshold = TRUE) {
    
    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("threshold must be a numeric..." = is.numeric(threshold) )
    stopifnot("above_threshold must be a logical..." = is.logical(above_threshold) )
    
    required_columns <- c("time", "value")
    assertthat::assert_that(
        identical(sort(colnames(data) ), sort(required_columns) ),
        msg = paste(
            "Data must have exactly two columns named 'time' and 'value'. Found:",
            paste(colnames(data), collapse = ", ")
            )
        )
    
    if (isFALSE(above_threshold) ) {
        
        threshold <- -threshold
        data <- data |> dplyr::mutate(value = if (above_threshold) .data$value else -.data$value)
        
    }
    
    clusters <- data |>
        dplyr::mutate(
            above = .data$value >= threshold,
            change = dplyr::lag(.data$above, default = FALSE) != .data$above,
            cluster_id = cumsum(.data$change & .data$above)
            ) |>
        dplyr::filter(.data$above) |>
        dplyr::group_by(.data$cluster_id) |>
        dplyr::summarise(
            cluster_onset = dplyr::first(.data$time),
            cluster_offset = dplyr::last(.data$time),
            .groups = "drop"
            ) |>
        data.frame()
    
    if (nrow(clusters) == 0) {
        clusters <- data.frame(
            cluster_id = NA_integer_,
            cluster_onset = NA_real_,
            cluster_offset = NA_real_
            )
        }
    
    return (clusters)
    
}

# defining the function in R (it will be executed in Python)
freq_stats_cluster_matrix <- function (
        X, timesteps,
        cluster_type = c("mass", "tfce"),
        alpha_level = 0.05,
        verbose = FALSE
        ) {
    
    # converting R matrix to NumPy array
    X_np <- np$array(X, dtype = "float64")
    
    # sanity check
    if (verbose) {message(print(X_np) )}
    
    # defining the type of cluster-based method to use
    if (cluster_type == "mass") {
        
        threshold <- py_none()
        
    } else {
        
        threshold <- dict(start = 0, step = 0.2)
        
    }
    
    # running the statistical test
    # results <- mne$stats$spatio_temporal_cluster_1samp_test(
    #     X = X_np,
    #     out_type = "mask",
    #     threshold = threshold,
    #     n_permutations = as.integer(2^12),
    #     n_jobs = as.integer(4),
    #     verbose = TRUE
    #     )
    results <- tryCatch({
        mne$stats$spatio_temporal_cluster_1samp_test(
            X = X_np,
            out_type = "mask",
            threshold = threshold,
            n_permutations = as.integer(2^12),
            n_jobs = as.integer(1),
            verbose = verbose
            )
    }, error = function (e) {
        message("Python error: ", e$message)
        message(reticulate::py_last_error() )
        return(NULL)
    })
    
    if (is.null(results) ) {
        message(reticulate::py_last_error() )
        stop("Python function failed.")
    }
    
    # sanity check
    if (verbose) {message("Length of results: ", length(results) )}
    
    if (length(results) < 3) {
        stop("Python function did not return 3 results as expected.")
    }
    
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
        
        # old
        # main_cluster_offset <- clusters[which(p_values < alpha_level)][[1]][[1]]$stop+1
        
        # assigning p-values to clusters
        for (i in seq_along(clusters) ) {
            
            # retrieving the current cluster
            cluster_mask <- clusters[[i]][[1]]
            idx_start <- cluster_mask$start+1
            idx_stop <- cluster_mask$stop
            
            # old
            # idx_stop <- cluster_mask$stop+1
            
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

# defining a function which prints a message using shell echo; useful for printing messages from inside mclapply when running in Rstudio
# from https://stackoverflow.com/questions/17345837/printing-from-mclapply-in-r-studio
message_parallel <- function (...) {
    
    system(sprintf('echo "\n%s\n"', paste0(..., collapse = "") ) )
    
}
