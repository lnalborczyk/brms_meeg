###########################################################
# Generate synthetic EEG data                             #
# ------------------------------------------------------- #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@cnrs.fr                    #
# Last updated on June 25, 2026                           #
###########################################################

################################################################################
# Adapted and extended from previous code written by Guillaume Rousselet
# R version + control variance: GAR, University of Glasgow, Oct 2023
# From https://github.com/GRousselet/onsetsim/blob/main/code/eeg_noise.R
# And https://github.com/GRousselet/onsetsim/blob/main/code/erp_template.R
# Implemented by: Rafal Bogacz and Nick Yeung, Princeton Univesity, December 2002
# Adapted from Matlab code:
# https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator
################################################################################

###########################################################
# Flexible version with arbitrary true clusters, SNR,     #
# participant-level amplitude variation, and              #
# participant-level noise variability.                    #
###########################################################

library(tidyverse)

simulate_eeg <- function (
    seed = 666,
    n_trials = 50,
    n_ppt = 20,
    outvar = 1,
    srate = 500,
    epoch_duration = 0.6,
    meanpower = NULL,
    meanpower_path = "code/meanpower.txt",

    # True effect scenario
    true_clusters = tibble(
        true_cluster_id = c(1, 2),
        true_onset = c(0.160, 0.410),
        true_offset = c(0.320, 0.490),
        relative_amplitude = c(1.0, 1.0)
        ),

    # Signal-to-noise ratio
    snr = 1,

    # Participant-level effect amplitude variability
    effect_amp_mean = 1,
    effect_amp_sd = 0.25,
    effect_amp_dist = c("truncated_normal", "lognormal", "fixed"),

    # Participant-level noise / overall variability
    noise_sd_mean = 1,
    noise_sd_sd = 0.5,
    noise_sd_dist = c("truncated_normal", "lognormal", "fixed"),

    # Baseline used to estimate noise scale for SNR
    baseline_window = c(0, 0.1),
    
    reference_condition = "cond1",
    target_condition = "cond2",

    # Return simulation data with settings and summaries
    return_list = TRUE
    ) {

    if (!is.null(seed) ) set.seed(seed)

    effect_amp_dist <- match.arg(effect_amp_dist)
    noise_sd_dist <- match.arg(noise_sd_dist)

    if (length(snr) != 1) {
        
        stop ("snr must be a single numeric value. For multiple SNR values, call simulate_eeg() once per SNR.")
        
    }

    if (!all(c("true_onset", "true_offset", "relative_amplitude") %in% names(true_clusters))) {
        stop(
            "true_clusters must contain at least: ",
            "true_onset, true_offset, relative_amplitude."
        )
    }

    if (!"true_cluster_id" %in% names(true_clusters)) {
        true_clusters <- true_clusters |>
            mutate(true_cluster_id = row_number())
    }

    true_clusters <- true_clusters |>
        arrange(true_onset) |>
        mutate(true_peak = true_onset + (true_offset - true_onset) / 2)

    if (any(true_clusters$true_offset <= true_clusters$true_onset)) {
        stop("Each true_offset must be greater than true_onset.")
    }

    if (min(true_clusters$true_onset) < 0 ||
        max(true_clusters$true_offset) > epoch_duration) {
        stop("All true clusters must fall within [0, epoch_duration].")
    }

    # ---------------------------------------------------------------------
    # Draw participant-level parameters
    # ---------------------------------------------------------------------

    draw_values <- function (n, mean, sd, dist, strictly_positive = TRUE) {

        if (dist == "fixed") {

            values <- rep(mean, n)

        } else if (dist == "truncated_normal") {

            if (sd == 0) {
                values <- rep(mean, n)
            } else {
                values <- rnorm(n, mean = mean, sd = sd)

                if (strictly_positive) {
                    while (any(values <= 0)) {
                        values[values <= 0] <- rnorm(
                            sum(values <= 0),
                            mean = mean,
                            sd = sd
                        )
                    }
                } else {
                    while (any(values < 0)) {
                        values[values < 0] <- rnorm(
                            sum(values < 0),
                            mean = mean,
                            sd = sd
                        )
                    }
                }
            }

        } else if (dist == "lognormal") {

            if (mean <= 0) {
                stop("For lognormal values, mean must be > 0.")
            }

            if (sd == 0) {
                values <- rep(mean, n)
            } else {
                sdlog <- sqrt(log(1 + (sd^2 / mean^2)))
                meanlog <- log(mean) - 0.5 * sdlog^2
                values <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
            }
        }

        if (strictly_positive && any(values <= 0)) {
            stop("Generated non-positive values where strictly positive values were required.")
        }

        if (!strictly_positive && any(values < 0)) {
            stop("Generated negative effect amplitudes.")
        }

        values
    }

    participant_params <- tibble(
        participant = paste0("participant_", sprintf("%02d", seq_len(n_ppt))),

        effect_amp = draw_values(
            n = n_ppt,
            mean = effect_amp_mean,
            sd = effect_amp_sd,
            dist = effect_amp_dist,
            strictly_positive = FALSE
        ),

        noise_sd = draw_values(
            n = n_ppt,
            mean = noise_sd_mean,
            sd = noise_sd_sd,
            dist = noise_sd_dist,
            strictly_positive = TRUE
        )
    ) |>
        mutate(
            noise_var = outvar * noise_sd^2
        )

    # ---------------------------------------------------------------------
    # EEG noise function
    # ---------------------------------------------------------------------

    eeg_noise <- function (frames = 51, srate = 100, outvar = 1, meanpower) {

        sumsig <- 50
        out <- vector(mode = "numeric", length = frames)
        freq <- 0

        for (i in seq_len(sumsig)) {

            freq <- freq + (4 * runif(1))
            freq_index <- min(ceiling(freq), length(meanpower))

            freqamp <- meanpower[freq_index] / meanpower[1]
            phase <- runif(1) * 2 * pi

            out <- out +
                sin(seq_len(frames) / srate * 2 * pi * freq + phase) *
                freqamp
        }

        out <- out - mean(out)
        out <- out / sd(out)
        out <- out * sqrt(outvar)

        return (out)
        
    }

    # ---------------------------------------------------------------------
    # Read or define meanpower
    # ---------------------------------------------------------------------

    if (is.null(meanpower)) {

        if (!is.null(meanpower_path) && file.exists(meanpower_path)) {

            meanpower <- unlist(read.table(meanpower_path))

        } else {

            warning(
                "No meanpower vector or valid meanpower_path provided. ",
                "Using a simple 1/f fallback spectrum. ",
                "For exact reproducibility with the original code, provide meanpower.txt."
            )

            meanpower <- 1 / seq_len(125)
        }
    }

    # ---------------------------------------------------------------------
    # Time vector and true effect profile
    # ---------------------------------------------------------------------

    time_vec <- seq(0, epoch_duration, by = 1 / srate)
    n_time <- length(time_vec)

    make_effect_profile <- function(time, true_clusters) {

        profile <- rep(0, length(time))

        for (cc in seq_len(nrow(true_clusters))) {

            onset <- true_clusters$true_onset[cc]
            offset <- true_clusters$true_offset[cc]
            rel_amp <- true_clusters$relative_amplitude[cc]

            inside <- time >= onset & time <= offset

            if (any(inside)) {
                phase <- (time[inside] - onset) / (offset - onset)

                # Raised-cosine bump: zero at onset/offset, one at peak
                bump <- 0.5 - 0.5 * cos(2 * pi * phase)

                profile[inside] <- profile[inside] + rel_amp * bump
            }
        }

        profile
    }

    effect_template <- tibble(
        time = time_vec,
        true_unit_effect = make_effect_profile(time_vec, true_clusters)
    )

    # ---------------------------------------------------------------------
    # Generate noise-only cond1 and cond2
    # ---------------------------------------------------------------------

    raw_df <- map_dfr(seq_len(n_ppt), function(P) {

        participant_id <- participant_params$participant[P]
        ppt_noise_var <- participant_params$noise_var[P]

        cond1 <- matrix(0, nrow = n_trials, ncol = n_time)
        cond2 <- matrix(0, nrow = n_trials, ncol = n_time)

        for (T in seq_len(n_trials)) {

            cond1[T, ] <- eeg_noise(
                frames = n_time,
                srate = srate,
                outvar = ppt_noise_var,
                meanpower = meanpower
            )

            cond2[T, ] <- eeg_noise(
                frames = n_time,
                srate = srate,
                outvar = ppt_noise_var,
                meanpower = meanpower
            )
        }

        tibble(
            participant = participant_id,
            condition = factor(
                rep(c(reference_condition, target_condition), each = n_time * n_trials),
                levels = c(reference_condition, target_condition)
            ),
            trial = rep(rep(seq_len(n_trials), each = n_time), times = 2),
            time = rep(time_vec, times = 2 * n_trials),
            eeg = c(c(t(cond1)), c(t(cond2)))
        )
    }) |>
        left_join(participant_params, by = "participant")

    # ---------------------------------------------------------------------
    # Compute baseline noise scale for SNR
    # ---------------------------------------------------------------------

    baseline_df <- raw_df |>
        filter(
            time >= baseline_window[1],
            time < baseline_window[2]
        )

    if (nrow(baseline_df) < 10) {
        warning(
            "Fewer than 10 rows found in baseline_window. ",
            "Using all non-effect time points to estimate baseline SD."
        )

        baseline_df <- raw_df |>
            filter(effect_template$true_unit_effect[match(time, effect_template$time)] == 0)
    }

    baseline_sd <- sd(baseline_df$eeg, na.rm = TRUE)

    if (!is.finite(baseline_sd) || baseline_sd <= 0) {
        baseline_sd <- sd(raw_df$eeg, na.rm = TRUE)
    }

    if (!is.finite(baseline_sd) || baseline_sd <= 0) {
        stop("Could not estimate a valid baseline SD.")
    }

    # ---------------------------------------------------------------------
    # Inject true condition effect into target condition
    # ---------------------------------------------------------------------

    raw_df <- raw_df |>
        left_join(effect_template, by = "time") |>
        mutate(
            true_effect = snr * baseline_sd * effect_amp * true_unit_effect,
            eeg = if_else(
                condition == target_condition,
                eeg + true_effect,
                eeg
            ),
            condition = factor(
                condition,
                levels = c(reference_condition, target_condition)
            )
        )

    # ---------------------------------------------------------------------
    # Derived truth summaries
    # ---------------------------------------------------------------------

    true_effect_by_time <- raw_df |>
        distinct(participant, time, true_effect) |>
        summarise(
            true_effect_mean = mean(true_effect),
            true_effect_sd = sd(true_effect),
            true_effect_min = min(true_effect),
            true_effect_max = max(true_effect),
            .by = time
        ) |>
        arrange(time)

    truth <- true_clusters |>
        mutate(
            snr = snr,
            baseline_sd = baseline_sd,
            n_true_clusters = n()
        )

    # ---------------------------------------------------------------------
    # Return objects
    # ---------------------------------------------------------------------

    if (return_list) {

        return(
            list(
                raw_df = raw_df,
                participant_params = participant_params,
                truth = truth,
                effect_template = effect_template,
                true_effect_by_time = true_effect_by_time,
                settings = list(
                    seed = seed,
                    n_trials = n_trials,
                    n_ppt = n_ppt,
                    outvar = outvar,
                    srate = srate,
                    epoch_duration = epoch_duration,
                    snr = snr,
                    baseline_window = baseline_window,
                    baseline_sd = baseline_sd,
                    reference_condition = reference_condition,
                    target_condition = target_condition
                )
            )
        )

    } else {

        return (raw_df)
        
    }
    
}
