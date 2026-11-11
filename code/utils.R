###########################################################
# Defining some utility functions                         #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@cnrs.fr                    #
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

freq_stats_cluster_matrix <- function(
    X,
    timesteps,
    cluster_type = c("mass", "tfce"),
    alpha_level = 0.05,
    n_permutations = 4096L,
    seed = 123L,
    tfce_start = 0,
    tfce_step = 0.2,
    verbose = FALSE
    ) {

    cluster_type <- match.arg(cluster_type)

    X <- as.matrix(X)
    storage.mode(X) <- "double"
    timesteps <- as.numeric(timesteps)

    if (length(dim(X)) != 2L) {
        stop("X must be a two-dimensional observations-by-time matrix.")
    }

    if (nrow(X) < 2L) {
        stop("X must contain at least two observations.")
    }

    if (ncol(X) != length(timesteps)) {
        stop(
            "ncol(X) = ", ncol(X),
            " but length(timesteps) = ", length(timesteps), "."
        )
    }

    if (any(!is.finite(X))) {
        stop("X contains non-finite values.")
    }

    np <- reticulate::import("numpy", convert = FALSE)
    mne <- reticulate::import("mne", convert = FALSE)

    X_np <- np$array(X, dtype = "float64")

    threshold_py <- if (cluster_type == "tfce") {
        reticulate::dict(
            start = as.numeric(tfce_start),
            step = as.numeric(tfce_step)
        )
    } else {
        reticulate::py_none()
    }

    results <- mne$stats$permutation_cluster_1samp_test(
        X = X_np,
        threshold = threshold_py,
        n_permutations = as.integer(n_permutations),
        tail = as.integer(0),
        out_type = "indices",
        n_jobs = as.integer(1),
        seed = as.integer(seed),
        verbose = verbose
    )

    # Python tuple:
    # T_obs, clusters, cluster_p_values, H0
    clusters_py <- reticulate::py_get_item(results, 1L)

    cluster_p_values <- reticulate::py_get_item(results, 2L) |>
        reticulate::py_to_r() |>
        as.numeric()

    # Default corrected p-value outside identified clusters
    pointwise_p_values <- rep(1, ncol(X))

    if (length(cluster_p_values) > 0L) {

        for (i in seq_along(cluster_p_values)) {

            cluster_py <- reticulate::py_get_item(
                clusters_py,
                as.integer(i - 1L)
            )

            indices_py <- reticulate::py_get_item(cluster_py, 0L)

            # Python is zero-indexed; R is one-indexed
            indices_r <- as.integer(
                reticulate::py_to_r(indices_py)
            ) + 1L

            indices_r <- indices_r[
                indices_r >= 1L &
                indices_r <= ncol(X)
            ]

            if (length(indices_r) > 0L) {
                pointwise_p_values[indices_r] <- pmin(
                    pointwise_p_values[indices_r],
                    cluster_p_values[i]
                )
            }
        }
    }

    if (
        length(pointwise_p_values) != length(timesteps) ||
        any(!is.finite(pointwise_p_values)) ||
        any(pointwise_p_values < 0 | pointwise_p_values > 1)
    ) {
        stop("Invalid corrected p-values returned by MNE.")
    }

    tibble::tibble(
        time = timesteps,
        pval = pointwise_p_values
    )
}

# defining a function which prints a message using shell echo; useful for printing messages from inside mclapply when running in Rstudio
# from https://stackoverflow.com/questions/17345837/printing-from-mclapply-in-r-studio
message_parallel <- function (...) {
    
    system(sprintf('echo "\n%s\n"', paste0(..., collapse = "") ) )
    
}
