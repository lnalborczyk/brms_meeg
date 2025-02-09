###########################################################
# Monte-Carlo simulation of onset/offset error properties #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on February 9, 2025                        #
###########################################################

library(changepoint)
library(tidyverse)
library(brms)

# retrieving the code from the qdm document
# knitr::purl(input = "brms_meeg.qmd", output = "brms_meeg.R", documentation = 0)

# importing R version of Matlab code from Yeung et al. (2004)
source("eeg_noise.R")

# importing the ERP template with true onset = 160 ms, F=81, and max at F=126
source("erp_template.R")

# importing helper functions
source("functions.R")

# to use with the eeg_noise function
meanpower <- unlist(read.table("meanpower.txt") )

# defining simulation parameters
n_trials <- 50 # number of trials
n_ppt <- 20 # number of participants
outvar <- 1 # noise variance
srate <- 500 # sampling rate in Hz
ronset <- seq(150, 170, 2) # random onset for each participant

# removing raw_df
if (exists("raw_df") ) rm(raw_df)

for (P in 1:n_ppt) { # for each participant
    
    # get random onset
    ponset <- sample(ronset, 1)
    
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

# defining the true onset and offset in seconds
true_onset <- 0.160
true_offset <- 0.342

# averaging across participants
raw_df %>%
    pivot_wider(names_from = condition, values_from = eeg, values_fn = mean) %>%
    mutate(eeg_diff = cond2 - cond1) %>%
    summarise(
        eeg_mean = mean(eeg_diff),
        eeg_sem = sd(eeg_diff) / sqrt(n() ),
        .by = time
        ) %>%
    mutate(lower = eeg_mean - eeg_sem, upper = eeg_mean + eeg_sem) %>%
    # plotting the data
    ggplot(aes(x = time, y = eeg_mean) ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = true_onset, linetype = 2) +
    geom_ribbon(
        aes(ymin = lower, ymax = upper, colour = NULL),
        alpha = 0.25, show.legend = FALSE
        ) +
    geom_line(show.legend = FALSE) +
    labs(x = "Time (s)", y = "EEG difference (a.u.)")

# averaging across participants
ppt_df <- raw_df %>%
    group_by(participant, condition, time) %>%
    summarise(eeg = mean(eeg) ) %>%
    # mutate(eeg = as.numeric(scale(x = eeg, scale = FALSE) ) ) %>%
    ungroup()

# defining a contrast for condition
contrasts(ppt_df$condition) <- c(-0.5, 0.5)

# fitting the GAM
gam <- brm(
    # cubic regression splines with k-1 basis functions
    eeg ~ condition + s(time, bs = "cr", k = 10, by = condition),
    data = ppt_df,
    family = gaussian(),
    iter = 5000,
    chains = 4,
    cores = 4,
    file = "models/gam.rds"
    )

gp_model <- brm(
    # k refers to the number of basis functions for
    # computing Hilbert-space approximate GPs
    # if k = NA (default), exact GPs are computed
    # eeg ~ gp(time, by = condition),
    eeg ~ condition + gp(time, by = condition, k = 20, cov = "exp_quad"),
    data = ppt_df,
    family = gaussian(),
    # prior = gp_priors,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    backend = "cmdstanr",
    iter = 2000,
    chains = 4,
    cores = 4,
    file = "models/gp.rds"
)

# defining a function to plot posterior predictions
plot_post_preds <- function (model, data = NULL) {
    
    if (is.null(data) ) data <- model$data
    
    post_preds <- conditional_effects(
        x = model,
        effect = "time:condition",
        method = "posterior_epred",
        re_formula = NULL,
        prob = 0.99
    )[[1]]
    
    post_preds %>%
        ggplot(aes(x = time, y = eeg, colour = condition, fill = condition) ) +
        geom_vline(xintercept = true_onset, linetype = 2) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(
            data = data %>% summarise(eeg = mean(eeg), .by = c(time, condition) ),
            aes(fill = NULL),
            pch = 21, show.legend = FALSE
        ) +
        geom_line(aes(y = estimate__), show.legend = FALSE) +
        geom_ribbon(
            aes(ymin = lower__, ymax = upper__, colour = NULL),
            alpha = 0.25, show.legend = FALSE
        ) +
        labs(x = "Time (s)", y = "EEG signal (a.u.)")
    
}

# plotting the posterior predictions
plot_post_preds(model = gam) + plot_post_preds(model = gp_model)

# check posterior_smooths()?
# https://www.rdocumentation.org/packages/brms/versions/2.22.0/topics/posterior_smooths.brmsfit
# defining a function to plot posterior slope
plot_post_slope <- function (model, data = NULL, eps = 0.1) {
    
    # if no data is specified, use model$data
    if (is.null(data) ) data <- model$data
    
    # defining a sequence of time values to make predictions
    time_seq <- crossing(
        time = seq(min(data$time), max(data$time), length.out = 100),
        condition = c("cond1", "cond2")
    ) %>%
        arrange(condition)
    
    # retrieving posterior samples
    posterior_samples <- posterior_epred(object = model, newdata = time_seq)
    
    # computing the difference between cond2 and cond1
    posterior_samples <- posterior_samples[, 101:200] - posterior_samples[, 1:100] 
    
    # computing the probability that y is above 0
    prob_y_above_0 <- data.frame(
        time = seq(min(data$time), max(data$time), length.out = 100)
    ) %>%
        mutate(
            m = colMeans(posterior_samples),
            lower = apply(
                X = posterior_samples,
                MARGIN = 2, FUN = quantile, probs = 0.025
            ),
            upper = apply(
                X = posterior_samples,
                MARGIN = 2, FUN = quantile, probs = 0.975
            )
        )
    
    # plotting it
    prob_y_above_0 %>%
        ggplot(aes(x = time, y = m) ) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_vline(xintercept = true_onset, linetype = 2) +
        geom_line() +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        labs(x = "Time (s)", y = expression(Pr(beta~"|"~data, model) ) )
    
}

# using the function with model
plot_post_slope(model = gam, eps = 0.1) +
    plot_post_slope(model = gp_model, eps = 0.1)

# defining a function to plot posterior prob of slope above 0
plot_post_test <- function (model, data = NULL, eps = 0.1) {
    
    # if no data is specified, use model$data
    if (is.null(data) ) data <- model$data
    
    # defining a sequence of time values to make predictions
    time_seq <- crossing(
        time = seq(min(data$time), max(data$time), length.out = 100),
        condition = c("cond1", "cond2")
    ) %>%
        arrange(condition)
    
    # retrieving posterior samples
    posterior_samples <- posterior_epred(object = model, newdata = time_seq)
    
    # computing the difference between cond2 and cond1
    posterior_samples <- posterior_samples[, 101:200] - posterior_samples[, 1:100]
    
    # computing the probability that y is above 0
    prob_y_above_0 <- data.frame(
        time = time_seq$time,
        m = colMeans(posterior_samples > (0 + eps) ),
        lower = apply(
            X = posterior_samples > (0 + eps),
            MARGIN = 2, FUN = quantile, probs = 0.025
        ),
        upper = apply(
            X = posterior_samples > (0 + eps),
            MARGIN = 2, FUN = quantile, probs = 0.975
        )
    )
    
    # plotting it
    prob_y_above_0 %>%
        ggplot(aes(x = time, y = m) ) +
        geom_hline(yintercept = 0.5, linetype = 2) +
        geom_vline(xintercept = true_onset, linetype = 2) +
        geom_vline(xintercept = true_offset, linetype = 2) +
        geom_line() +
        labs(x = "Time (s)", y = expression(Pr(beta>0~"|"~data, model) ) )
    
}

# using the function with the two models
plot_post_test(model = gam, eps = 0.1) +
    plot_post_test(model = gp_model, eps = 0.1)

# defining a sequence of x values to make predictions
n_approx <- 1000
time_seq <- crossing(
    time = seq(min(ppt_df$time), max(ppt_df$time), length.out = n_approx),
    condition = c("cond1", "cond2")
) %>%
    arrange(condition)

# retrieving posterior samples
posterior_samples <- posterior_epred(object = gam, newdata = time_seq)

# computing the difference between cond2 and cond1
posterior_samples <- posterior_samples[, (ncol(posterior_samples)/2+1):ncol(posterior_samples)] -
    posterior_samples[, 1:(ncol(posterior_samples)/2)]

# computing the probability that y is above 0
eps <- 0.05
prob_y_above_0 <- data.frame(
    time = seq(min(ppt_df$time), max(ppt_df$time), length.out = n_approx)
) %>%
    mutate(post_prob = colMeans(posterior_samples) ) %>%
    mutate(m = colMeans(posterior_samples > eps) ) %>%
    mutate(prob_ratio = m / (1 - m) )

# printing the identified clusters
# the signal was a truncated Gaussian defining an objective onset at 160 ms,
# a maximum at 250 ms, and an offset at 342 ms.
threshold <- 10
exceeding_times_gam <- prob_y_above_0 %>%
    dplyr::filter(prob_ratio > threshold) %>%
    summarise(cluster_onset = min(time), cluster_offset = max(time) ) %>%
    mutate(
        cluster_peak = prob_y_above_0 %>%
            dplyr::filter(post_prob == max(post_prob) ) %>%
            pull(time)
    ) %>%
    mutate(true_onset = 0.160, true_peak = 0.250, true_offset = 0.342) %>%
    mutate(model = "GAM")

# plotting it
gam_ratio_plot <- prob_y_above_0 %>%
    mutate(above_thres = ifelse(prob_ratio > threshold, 1, NA) ) %>%
    ggplot(aes(x = time, y = prob_ratio) ) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 3, linetype = "dashed", color = "darkred") +
    geom_hline(yintercept = 1/3, linetype = "dashed", color = "darkred") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 1/10, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 100, linetype = "dashed", color = "orangered") +
    geom_hline(yintercept = 1/100, linetype = "dashed", color = "orangered") +
    geom_vline(xintercept = true_onset, linetype = 2) +
    geom_vline(xintercept = true_offset, linetype = 2) +
    geom_line() +
    geom_point(
        data = . %>% dplyr::filter(!is.na(above_thres) ),
        aes(y = threshold),
        colour = "darkgreen",
        shape = 15,
        size = 3,
        na.rm = TRUE
    ) +
    scale_y_log10(labels = label_log(digits = 2), limits = c(NA, 1e4) ) +
    labs(
        x = "Time (s)",
        y = expression(Pr(beta>0~"|"~data) / (1 - Pr(beta>0~"|"~data) ) )
    )

# retrieving posterior samples
posterior_samples <- posterior_epred(object = gp_model, newdata = time_seq)

# computing the difference between cond2 and cond1
posterior_samples <- posterior_samples[, (ncol(posterior_samples)/2+1):ncol(posterior_samples)] -
    posterior_samples[, 1:(ncol(posterior_samples)/2)]

# computing the probability that y is above 0
prob_y_above_0 <- data.frame(
    time = seq(min(ppt_df$time), max(ppt_df$time), length.out = n_approx)
) %>%
    mutate(post_prob = colMeans(posterior_samples) ) %>%
    mutate(m = colMeans(posterior_samples > eps) ) %>%
    mutate(prob_ratio = m / (1 - m) )

# printing the identified clusters
# the signal was a truncated Gaussian defining an objective onset at 160 ms,
# a maximum at 250 ms, and an offset at 342 ms
exceeding_times_gp <- prob_y_above_0 %>%
    dplyr::filter(prob_ratio > threshold) %>%
    summarise(cluster_onset = min(time), cluster_offset = max(time) ) %>%
    mutate(
        cluster_peak = prob_y_above_0 %>%
            dplyr::filter(post_prob == max(post_prob) ) %>%
            pull(time)
    ) %>%
    mutate(true_onset = 0.160, true_peak = 0.250, true_offset = 0.342) %>%
    mutate(model = "GP")

# plotting it
gp_ratio_plot <- prob_y_above_0 %>%
    mutate(above_thres = ifelse(prob_ratio > threshold, 1, NA) ) %>%
    ggplot(aes(x = time, y = prob_ratio) ) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 3, linetype = "dashed", color = "darkred") +
    geom_hline(yintercept = 1/3, linetype = "dashed", color = "darkred") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 1/10, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 100, linetype = "dashed", color = "orangered") +
    geom_hline(yintercept = 1/100, linetype = "dashed", color = "orangered") +
    geom_vline(xintercept = true_onset, linetype = 2) +
    geom_vline(xintercept = true_offset, linetype = 2) +
    geom_line() +
    geom_point(
        data = . %>% dplyr::filter(!is.na(above_thres) ),
        aes(y = threshold),
        colour = "darkgreen",
        shape = 15,
        size = 3,
        na.rm = TRUE
    ) +
    scale_y_log10(labels = label_log(digits = 2), limits = c(NA, 1e4) ) +
    labs(
        x = "Time (s)",
        y = expression(Pr(beta>0~"|"~data) / (1 - Pr(beta>0~"|"~data) ) )
    )

# plotting results for the two models
gam_ratio_plot + gp_ratio_plot

# printing the results
# bind_rows(exceeding_times_gam, exceeding_times_gp) %>%
#     mutate(across(.cols = where(is.numeric), .fns = ~round(.x, 2) ) ) %>%
#     relocate(model, .before = cluster_onset)

# defining values of eps and threshold to test
eps_values <- seq(from = 0, to = 0.2, by = 0.01)
threshold_values <- seq(from = 1, to = 50, by = 1)

# using as many values of time as in the original data
n_approx <- n_distinct(ppt_df$time)

# initialising results dataframe
results <- crossing(eps = eps_values, threshold = threshold_values) %>%
    mutate(
        estimated_onset = NA, estimated_peak = NA, estimated_offset = NA,
        error_onset = NA, error_peak = NA, error_offset = NA
    )

# generating a grid of time values to compute predictions for
time_seq <- crossing(
    # participant = unique(ppt_df$participant),
    condition = c("cond1", "cond2"),
    time = seq(min(ppt_df$time), max(ppt_df$time), length.out = n_approx*2-1),
) %>%
    arrange(condition)

# retrieving posterior samples
# posterior_samples <- posterior_epred(object = gp_model, newdata = time_seq)
posterior_samples <- posterior_epred(object = gam, newdata = time_seq)

# computing the difference between cond2 and cond1
posterior_samples <- posterior_samples[, (ncol(posterior_samples)/2+1):ncol(posterior_samples)] -
    posterior_samples[, 1:(ncol(posterior_samples)/2)]

# looping over different values of eps and threshold
for (i in seq_len(nrow(results) ) ) {
    
    # retrieving current eps and threshold values
    eps <- results$eps[i]
    threshold <- results$threshold[i]
    
    # computing probability metrics with the current eps value
    prob_y_above_0 <- data.frame(
        time = seq(min(ppt_df$time), max(ppt_df$time), length.out = n_approx*2-1)
    ) %>%
        mutate(post_prob = colMeans(posterior_samples) ) %>%
        mutate(m = colMeans(posterior_samples > eps) ) %>%
        mutate(prob_ratio = m / (1 - m) )
    
    # finding onset, offset, and peak for the current threshold
    exceeding_times <- prob_y_above_0 %>%
        dplyr::filter(prob_ratio > threshold) %>%
        summarise(
            cluster_onset = min(time, na.rm = TRUE),
            cluster_offset = max(time, na.rm = TRUE)
        ) %>%
        mutate(
            cluster_peak = prob_y_above_0 %>%
                dplyr::filter(post_prob == max(post_prob, na.rm = TRUE) ) %>%
                pull(time)
        ) %>%
        mutate(true_onset = 0.160, true_peak = 0.250, true_offset = 0.342)
    
    # storing the results in the dataframe
    if (nrow(exceeding_times) > 0) {
        
        results$estimated_onset[i] <- exceeding_times$cluster_onset
        results$estimated_peak[i] <- exceeding_times$cluster_peak
        results$estimated_offset[i] <- exceeding_times$cluster_offset
        
        # computing errors
        results$error_onset[i] <- abs(
            exceeding_times$cluster_onset - exceeding_times$true_onset
        )
        results$error_peak[i] <- abs(
            exceeding_times$cluster_peak - exceeding_times$true_peak
        )
        results$error_offset[i] <- abs(
            exceeding_times$cluster_offset - exceeding_times$true_offset
        )
        
    }
    
}

# plotting the results
onset_plot <- results %>%
    ggplot(aes(x = eps, y = threshold, fill = error_onset) ) +
    geom_tile(show.legend = FALSE) +
    geom_point(
        data = results %>% dplyr::filter(error_onset == min(error_onset, na.rm = TRUE) ),
        aes(x = eps, y = threshold),
        color = "orangered", size = 2, shape = 4,
        show.legend = FALSE
    ) +
    scale_x_continuous(expand = c(0, 0) ) +
    scale_y_continuous(expand = c(0, 0) ) +
    scale_fill_gradientn(colors = rev(met.brewer("Hokusai1") ) ) +
    labs(
        title = "Onset error",
        x = "eps",
        y = "Threshold",
        fill = "Error"
    )

offset_plot <- results %>%
    ggplot(aes(x = eps, y = threshold, fill = error_offset) ) +
    geom_tile(show.legend = FALSE) +
    geom_point(
        data = results %>% dplyr::filter(error_offset == min(error_offset, na.rm = TRUE) ),
        aes(x = eps, y = threshold),
        color = "orangered", size = 2, shape = 4,
        show.legend = FALSE
    ) +
    scale_x_continuous(expand = c(0, 0) ) +
    scale_y_continuous(expand = c(0, 0) ) +
    scale_fill_gradientn(colors = rev(met.brewer("Hokusai1") ) ) +
    labs(
        title = "Offset error",
        x = "eps",
        y = "Threshold",
        fill = "Error"
    )

# combining the plots
onset_plot + offset_plot

# displaying best parameters values (according to onset error)
results %>%
    dplyr::select(-estimated_peak, -error_peak) %>%
    arrange(error_onset) %>%
    data.frame() %>%
    head()

# displaying best parameters values (according to offset error)
results %>%
    dplyr::select(-estimated_peak, -error_peak) %>%
    arrange(error_offset) %>%
    data.frame() %>%
    head()

# displaying best parameters values (according to both onset/offset error)
results %>%
    filter(error_onset == min(error_onset) & error_offset == min(error_offset) ) %>%
    data.frame()

# averaging across participants
summary_df <- raw_df %>%
    summarise(
        eeg_mean = mean(eeg),
        eeg_sd = sd(eeg),
        .by = c(participant, condition, time)
    )

# defining a contrast for condition
contrasts(summary_df$condition) <- c(-0.5, 0.5)

# fitting the GAM
meta_gam <- brm(
    # using by-participant SD of ERPs across trials
    eeg_mean | se(eeg_sd) ~
        condition + s(time, bs = "cr", k = 10, by = condition) +
        (1 | participant),
    data = summary_df,
    family = gaussian(),
    iter = 5000,
    chains = 4,
    cores = 4,
    file = "models/meta_gam.rds"
)

# fitting the GP
meta_gp <- brm(
    # using by-participant SD of ERPs across trials
    eeg_mean | se(eeg_sd) ~
        condition + gp(time, k = 20, by = condition) +
        (1 | participant),
    data = summary_df,
    family = gaussian(),
    control = list(adapt_delta = 0.99),
    iter = 5000,
    chains = 4,
    cores = 4,
    file = "models/meta_gp.rds"
)

# plotting the posterior predictions
meta_gam_preds <- plot(
    conditional_effects(x = meta_gam, effect = "time:condition"),
    points = FALSE, theme = theme_light(), plot = FALSE
)[[1]] +
    labs(x = "Time (s)", y = "EEG signal (a.u.)")
meta_gp_preds <- plot(
    conditional_effects(x = meta_gp, effect = "time:condition"),
    points = FALSE, theme = theme_light(), plot = FALSE
)[[1]] +
    labs(x = "Time (s)", y = "EEG signal (a.u.)")

# combining the two plots
meta_gam_preds + meta_gp_preds

# simulating decoding accuracies (ROC AUC)
n_participants <- 20
n_timepoints <- 240
timepoints <- seq(from = 0.1, to = 1.3, length.out = n_timepoints) 

# defining a function to generate AUC using a Beta distribution
simulate_auc <- function (
        timepoints, onset = 0.1,
        peak_time = 0.4, sigma = 0.5, amplitude = 30
) {
    
    # ensuring timepoints - onset > 0
    adjusted_time <- timepoints - onset + 0.001
    
    # log-normal peak function with right skew
    log_peak <- log(peak_time - onset + 0.001)
    log_t <- log(adjusted_time)
    
    # computing right-skewed alpha with log-normal peak
    alpha_t <- 10 + amplitude * exp(-((log_t - log_peak)^2) / (2 * sigma^2) )
    
    # keeping baseline at chance level (0.5)
    beta_t <- 10
    
    # simulate Beta-distributed AUC values
    auc <- rbeta(length(timepoints), shape1 = alpha_t, shape2 = beta_t)
    
    # returning the AUC
    return (auc)
    
}

# simulating data for all participants
decoding_data <- data.frame(
    participant = rep(1:n_participants, each = n_timepoints),
    time = rep(timepoints, times = n_participants),
    auc = unlist(lapply(1:n_participants, function (i) {
        simulate_auc(
            timepoints,
            onset = 0.1,
            # variable peak time
            peak_time = runif(1, 0.6, 0.7),
            # variable spread
            sigma = runif(1, 0.1, 0.3),
            # variable peak amplitude
            amplitude = runif(1, 30, 40)
        )})
    )
) %>%
    mutate(time = time - 0.3)

# computing group mean and 95% quantile interval
group_data <- decoding_data %>%
    group_by(time) %>%
    summarise(
        mean_auc = mean(auc),
        lower = quantile(auc, 0.025),
        upper = quantile(auc, 0.975)
    )

# plotting it
# group_data %>%
#     ggplot(aes(x = time, y = mean_auc) ) +
#     geom_hline(yintercept = 0.5, linetype = "dashed") +
#     geom_vline(xintercept = 0.0, linetype = "dashed") +
#     geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
#     geom_line(color = "steelblue", size = 1) +
#     ylim(0, 1) +
#     labs(x = "Time (s)", y = "Decoding accuracy (ROC AUC)")

## # loading the reticulate package
## library(reticulate)
## 
## # importing the numpy and mne python modules
## np <- import("numpy")
## mne <- import("mne")
## sklearn <- import("sklearn")
## 
## # defining the function in R (it will be executed in Python)
## mne_decoding <- function (X, labels, ncores = 8) {
## 
##     # converting R dataframe to NumPy array (reshaping if needed)
##     # X should be a matrix before conversion
##     X_np <- np$array(X)
## 
##     if (length(dim(X_np) ) == 2) {
## 
##         # adding a second dimension (channels) if missing
##         X_np <- np$expand_dims(X_np, axis = as.integer(1) )
## 
##     }
## 
##     # defining the classifier
##     clf <- sklearn$linear_model$LogisticRegression(solver = "liblinear")
## 
##     # sliding the estimator on all time frames
##     time_decod <- mne$decoding$SlidingEstimator(
##         clf, n_jobs = as.integer(ncores),
##         scoring = "roc_auc", verbose = TRUE
##         )
## 
##     # or using N-fold cross-validation
##     scores <- mne$decoding$cross_val_multiscore(
##         time_decod,
##         X_np,
##         labels,
##         cv = as.integer(4),
##         n_jobs = as.integer(ncores),
##         verbose = TRUE
##         )
## 
##     # returning the scores (averaged over CV folds)
##     return (scores)
## 
## }
## 
## # listing all participants
## participants <- unique(raw_df$participant)
## 
## # initialising empty decoding results
## group_decoding_scores <- data.frame()
## 
## # running decoding for each participant
## for (ppt in participants) {
## 
##     # printing progress
##     print(ppt)
## 
##     # retrieve data from one participant
##     ppt_data <- raw_df %>%
##         filter(participant == ppt) %>%
##         select(-participant) %>%
##         pivot_wider(names_from = time, values_from = eeg) %>%
##         select(-condition, -trial)
## 
##     # extracting the labels
##     labels <- raw_df %>%
##         filter(participant == ppt) %>%
##         select(-participant) %>%
##         pivot_wider(names_from = time, values_from = eeg) %>%
##         pull(condition) %>%
##         as.numeric()
## 
##     # extracting the timesteps
##     timesteps <- raw_df %>%
##         filter(participant == ppt) %>%
##         select(-participant) %>%
##         pull(time) %>%
##         unique()
## 
##     # running the decoding
##     decoding_scores <- data.frame(
##         mne_decoding(X = ppt_data, labels = labels-1)
##         ) %>%
##         # computing the average over CV folds
##         summarise(across(where(is.numeric), mean) )
## 
##     # appending to previous results
##     group_decoding_scores <- bind_rows(group_decoding_scores, decoding_scores)
## 
## }
## 
## # saving the scores
## saveRDS(object = group_decoding_scores, file = "results/decoding_scores.rds")

# importing the decoding scores
group_decoding_scores <- readRDS(file = "results/decoding_scores.rds")

# extracting the timesteps
timesteps <- raw_df %>%
    pull(time) %>%
    unique()

# plotting it
decoding_summary <- group_decoding_scores %>%
    t() %>%
    data.frame() %>%
    mutate(
        mean = rowMeans(across(1:20) ),
        se = apply(across(1:20), 1, function (x) sd(x) / sqrt(length(x) ) )
    ) %>%
    mutate(lower = mean - se, upper = mean + se) %>%
    mutate(time = timesteps) %>%
    mutate(
        ymin = pmin(mean, 0.5),
        ymax = pmax(mean, 0.5)
    )

# reshaping the decoding data
decoding_data <- group_decoding_scores %>%
    t() %>%
    data.frame() %>%
    mutate(time = timesteps) %>%
    pivot_longer(cols = 1:20) %>%
    mutate(
        participant = rep(
            unique(raw_df$participant),
            n_distinct(raw_df$time)
        )
    ) %>%
    select(participant, time, auc = value)

# plotting it
# group_decoding_scores %>%
#     t() %>%
#     data.frame() %>%
#     mutate(time = timesteps) %>%
#     pivot_longer(cols = 1:20) %>%
#     mutate(
#         participant = rep(
#             unique(raw_df$participant),
#             n_distinct(raw_df$time)
#             )
#         ) %>%
decoding_data %>%
    ggplot(aes(x = time) ) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_vline(xintercept = true_onset, linetype = "dashed") +
    geom_line(aes(y = auc), linewidth = 0.5) +
    facet_wrap(~participant) +
    geom_ribbon(
        data = decoding_summary,
        aes(x = time, ymin = ymin, ymax = ymax),
        alpha = 0.5
    ) +
    labs(x = "Time (s)", y = "Decoding accuracy (ROC AUC)")

# fitting the GAM
decoding_gam <- brm(
    auc ~ s(time, bs = "cr", k = 10),
    data = decoding_data,
    family = Beta(),
    iter = 5000,
    chains = 4,
    cores = 4,
    file = "models/decoding_gam.rds"
)

# fitting the GP
decoding_gp <- brm(
    auc ~ gp(time, k = 20),
    data = decoding_data,
    family = Beta(),
    control = list(adapt_delta = 0.99),
    iter = 5000,
    chains = 4,
    cores = 4,
    file = "models/decoding_gp.rds"
)

# plotting the posterior predictions
gam_decoding_preds <- plot(
    conditional_effects(x = decoding_gam, prob = 0.99),
    points = TRUE,
    point_args = list(size = 1, pch = 21, colour = "grey"),
    line_args = list(colour = "black"),
    theme = theme_light(), plot = FALSE
)[[1]] +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_vline(xintercept = true_onset, linetype = 2) +
    ylim(0, 1) +
    labs(x = "Time (s)", y = "Decoding accuracy (ROC AUC)")

gp_decoding_preds <- plot(
    conditional_effects(x = decoding_gp, prob = 0.99),
    points = TRUE,
    point_args = list(size = 1, pch = 21, colour = "grey"),
    line_args = list(colour = "black"),
    theme = theme_light(), plot = FALSE
)[[1]] +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_vline(xintercept = true_onset, linetype = 2) +
    ylim(0, 1) +
    labs(x = "Time (s)", y = "Decoding accuracy (ROC AUC)")

# combining the two plots
gam_decoding_preds + gp_decoding_preds

# defining a function to plot posterior prob of slope above 0
plot_post_decoding <- function (model, data, chance_level = 0.5, eps = 0.1) {
    
    # defining a sequence of time values to make predictions
    time_seq <- crossing(time = seq(min(data$time), max(data$time), length.out = 100) )
    
    # retrieving posterior samples
    posterior_samples <- posterior_epred(object = model, newdata = time_seq)
    
    # computing the probability that y is above 0
    prob_y_above_0 <- data.frame(
        time = time_seq$time,
        m = colMeans(posterior_samples > (chance_level + eps) ),
        lower = apply(
            X = posterior_samples > (chance_level + eps),
            MARGIN = 2, FUN = quantile, probs = 0.025
        ),
        upper = apply(
            X = posterior_samples > (chance_level + eps),
            MARGIN = 2, FUN = quantile, probs = 0.975
        )
    )
    
    # plotting it
    prob_y_above_0 %>%
        ggplot(aes(x = time, y = m) ) +
        geom_hline(yintercept = chance_level, linetype = 2) +
        geom_vline(xintercept = true_onset, linetype = 2) +
        geom_line() +
        labs(x = "Time (s)", y = expression(Pr(AUC>chance~"|"~data, model) ) )
    
}

# using the function with the two models
plot_post_decoding(
    model = decoding_gam, data = decoding_data,
    chance_level = 0.5, eps = 0.1
) +
    plot_post_decoding(
        model = decoding_gp, data = decoding_data,
        chance_level = 0.5, eps = 0.1
    ) 

# defining a function to plot posterior prob of slope above 0
plot_post_decoding_ratio <- function (model, data, chance_level = 0.5, eps = 0.1) {
    
    # defining a sequence of time values to make predictions
    time_seq <- data.frame(time = data$time)
    
    # retrieving posterior samples
    posterior_samples <- posterior_epred(object = model, newdata = time_seq)
    
    # computing the probability that y is above 0
    prob_y_above_0 <- data.frame(time = data$time) %>%
        mutate(m = colMeans(posterior_samples > (chance_level + eps) ) ) %>%
        mutate(prob_ratio = m / (1 - m) )
    
    # printing the identified clusters
    threshold <- 10
    exceeding_times <- prob_y_above_0 %>%
        dplyr::filter(prob_ratio > threshold) %>%
        summarise(cluster_onset = min(time), cluster_offset = max(time) )
    
    # printing the identified clusters
    # print(exceeding_times)
    
    # plotting it
    prob_y_above_0 %>%
        mutate(above_thres = ifelse(prob_ratio > threshold, 1, NA) ) %>%
        ggplot(aes(x = time, y = prob_ratio) ) +
        geom_hline(yintercept = 1, linetype = "dashed") +
        geom_hline(yintercept = 3, linetype = "dashed", color = "darkred") +
        geom_hline(yintercept = 1/3, linetype = "dashed", color = "darkred") +
        geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
        geom_hline(yintercept = 1/10, linetype = "dashed", color = "red") +
        geom_hline(yintercept = 100, linetype = "dashed", color = "orangered") +
        geom_hline(yintercept = 1/100, linetype = "dashed", color = "orangered") +
        geom_vline(xintercept = true_onset, linetype = 2) +
        geom_line() +
        geom_point(
            data = . %>% dplyr::filter(!is.na(above_thres) ),
            aes(y = threshold),
            colour = "darkgreen",
            shape = 15,
            size = 3,
            na.rm = TRUE
        ) +
        scale_y_log10(labels = label_log(digits = 2), limits = c(NA, 1e4) ) +
        labs(
            x = "Time (s)",
            y = expression(
                Pr(AUC>chance~"|"~data) / (1 - Pr(AUC>chance~"|"~data) )
            )
        )
    
}

# using the function with the two models
plot_post_decoding_ratio(
    model = decoding_gam, data = decoding_data,
    chance_level = 0.5, eps = 0.05
) +
    plot_post_decoding_ratio(
        model = decoding_gp, data = decoding_data,
        chance_level = 0.5, eps = 0.05
    )

## # centering AUC relative to chance level (0.5)
## decoding_data_perm <- decoding_data %>%
##     mutate(auc_centered = auc - 0.5) %>%
##     dplyr::select(-auc)
## 
## # formatting the AUC data to be used with permuco
## dv <- decoding_data_perm %>%
##     pivot_wider(names_from = time, values_from = auc_centered) %>%
##     dplyr::select(-participant) %>%
##     data.frame()
## 
## # formatting the intercept to be used with permuco
## iv <- decoding_data %>%
##     dplyr::select(participant) %>%
##     distinct() %>%
##     mutate(participant = as.factor(participant) ) %>%
##     data.frame() %>%
##     mutate(intercept = 1)
## 
## # performing cluster-based permutation test
## cluster_results <- clusterlm(
##     formula = dv ~ intercept,
##     data = iv,
##     # data = decoding_data_perm,
##     # formula = auc_centered ~ participant,
##     # number of permutations (increase for accuracy)
##     # np = 2000,
##     # type of transformation
##     type = "signflip"
##     # using a t-test
##     # test = "t"
##     # multiple comparison correction
##     # multcomp = "clustermass"
##     # multcomp = "clusterdepth_head",
##     # multcomp = "tfce"
##     )
## 
## # printing results
## print(cluster_results)
## 
## # plotting significant clusters
## plot(cluster_results)
## 
## # retrieving t and p-values
## f_values <- cluster_results$multiple_comparison$intercept$clustermass$main[, 1]
## p_values <- cluster_results$multiple_comparison$intercept$clustermass$main[, 2]

# https://github.com/GRousselet/onsetsim/blob/main/docs/weakstrongfwer_demo.md
# BH95: p.adjust(pvals, method = "BH")
# BY01: p.adjust(perm.pvals, method = "BY")

# defining an arbitrary alpha threshold
aath <- 0.05

# summarising raw_data per participant
ppt_data <- raw_df %>%
    summarise(eeg = mean(eeg), .by = c(participant, condition, time) ) %>%
    pivot_wider(names_from = condition, values_from = eeg) %>%
    mutate(eeg_diff = cond2 - cond1)

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
        # pval_bonf = p.adjust(p = pval, method = "bonferroni"),
        pval_holm = p.adjust(p = pval, method = "holm")
    ) %>%
    ungroup()

# finding onsets
from_timestep <- 1
onset_p <- find_onset(
    tests_results$pval[from_timestep:nrow(tests_results)] <= aath,
    tests_results$time[from_timestep:nrow(tests_results)]
)
onset_bh <- find_onset(
    tests_results$pval_bh[from_timestep:nrow(tests_results)] <= aath,
    tests_results$time[from_timestep:nrow(tests_results)]
)
onset_by <- find_onset(
    tests_results$pval_by[from_timestep:nrow(tests_results)] <= aath,
    tests_results$time[from_timestep:nrow(tests_results)]
)
# onset_bonf <- find_onset(
#     tests_results$pval_bonf[from_timestep:nrow(tests_results)] <= aath,
#     tests_results$time[from_timestep:nrow(tests_results)]
#     )
onset_holm <- find_onset(
    tests_results$pval_holm[from_timestep:nrow(tests_results)] <= aath,
    tests_results$time[from_timestep:nrow(tests_results)]
)

# plotting the results
default_colours <- hue_pal()(10)
tvals_plot <- tests_results %>%
    ggplot(aes(x = time, y = tval) ) +
    geom_area(position = "identity") +
    geom_vline(xintercept = true_onset, linetype = 2) +
    geom_vline(xintercept = onset_p, linetype = 2, color = default_colours[1]) +
    geom_vline(xintercept = onset_bh, linetype = 2, color = default_colours[2]) +
    geom_vline(xintercept = onset_by, linetype = 2, color = default_colours[4]) +
    # geom_vline(xintercept = onset_bonf, linetype = 2, color = default_colours[3]) +
    geom_vline(xintercept = onset_holm, linetype = 2, color = default_colours[5]) +
    labs(x = "Time (s)", y = bquote(t^2) )

plot_title <- paste0(
    "True onset: ", true_onset, "s, ",
    "Uncorrected onset: ", round(onset_p, 3), "s, ",
    "BH onset: ", round(onset_bh, 3), "s, ",
    "BY onset: ", round(onset_by, 3), "s, ",
    # "Bonf onset: ", round(onset_bonf, 3), "s, ",
    "Holm onset: ", round(onset_holm, 3),  "s"
)

# plotting -log2(p) (surprisal under H0) with thresholds (on p-values)
# https://lesslikely.com/statistics/s-values/
# a s-value of 4 (bits of information) is no more surprising than
# getting all heads on 4 fair coin tosses
pvals_plot <- tests_results %>%
    pivot_longer(cols = pval:pval_holm) %>%
    mutate(sval = -log2(value) ) %>%
    ggplot(aes(x = time, y = sval, colour = name) ) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = -log2(aath), linetype = 2) +
    geom_vline(xintercept = true_onset, linetype = 2) +
    geom_vline(xintercept = onset_p, linetype = 2, color = default_colours[1]) +
    geom_vline(xintercept = onset_bh, linetype = 2, color = default_colours[2]) +
    geom_vline(xintercept = onset_by, linetype = 2, color = default_colours[4]) +
    # geom_vline(xintercept = onset_bonf, linetype = 2, color = default_colours[3]) +
    geom_vline(xintercept = onset_holm, linetype = 2, color = default_colours[5]) +
    labs(
        x = "Time (s)",
        y = "-log2(p-value)",
        colour = "correction"
    )

# combining the two plots
tvals_plot + pvals_plot +
    plot_annotation(title = plot_title) & theme(plot.title = element_text(hjust = 0.5) )

# using the changepoint package to identify onsets and offsets
res <- cpt.meanvar(data = tests_results$tval, method = "BinSeg", Q = 2)

# plotting the results
tests_results %>%
    ggplot(aes(x = time, y = tval) ) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = true_onset, linetype = 2) +
    geom_vline(xintercept = true_offset, linetype = 2) +
    geom_vline(
        xintercept = tests_results$time[res@cpts[1]],
        linetype = 3, color = "darkgreen"
    ) +
    geom_vline(
        xintercept = tests_results$time[res@cpts[2]],
        linetype = 3, color = "darkgreen"
    ) +
    labs(x = "Time (s)", y = bquote(t^2) ) +
    ggtitle(
        paste0(
            "True onset: ", true_onset, "s, ", "Change point onset: ",
            round(tests_results$time[res@cpts[1]], 3), "s, ",
            "True offset: ", true_offset, "s, ", "Change point offset: ",
            round(tests_results$time[res@cpts[2]], 3), "s"
        )
    )

## # loading the reticulate package
## library(reticulate)
## 
## # importing the numpy and mne python modules
## np <- import("numpy")
## mne <- import("mne")
## 
## # defining the function in R (it will be executed in Python)
## freq_stats_gat_matrix <- function (X) {
## 
##     # converting R matrix to NumPy array
##     X_np <- np$array(X, dtype = "float64")
## 
##     # running the statistical test
##     results <- mne$stats$spatio_temporal_cluster_1samp_test(
##         X_np,
##         out_type = "mask",
##         n_permutations = as.integer(1000),
##         n_jobs = as.integer(4),
##         verbose = TRUE
##         )
## 
##     # extracting results
##     T_obs_ <- results[[1]]
##     clusters <- results[[2]]
##     p_values <- results[[3]]
## 
##     # using the first slice of the 3D array
##     p_values_ <- np$transpose(np$ones_like(X_np[1, , drop = FALSE]) )
## 
##     # assigning p-values to clusters
##     for (i in seq_along(clusters) ) {
## 
##         # retrieving the current cluster
##         cluster_mask <- clusters[[i]][[1]]
##         idx_start <- cluster_mask$start+1
##         idx_stop <- cluster_mask$stop
## 
##         # retrieving the p-value for this cluster
##         pval <- p_values[[i]]
## 
##         # assigning this p-value to timesteps belonging
##         # to the current cluster
##         p_values_[idx_start:idx_stop] <- pval
## 
##     }
## 
##     # converting result back to R format
##     return (as.matrix(np$squeeze(np$transpose(p_values_) ) ) )
## 
## }
## 
## # converting decoding_data to a matrix
## data_matrix <- matrix(
##     decoding_data$auc,
##     ncol = length(unique(decoding_data$time) ),
##     byrow = TRUE
##     )
## 
## # running the function
## chance_level <- 0.5
## p_values_results <- freq_stats_gat_matrix(data_matrix-chance_level)
## 
## # converting back to a dataframe if needed
## p_values_df <- data.frame(
##     pval = p_values_results,
##     time = decoding_data$time
##     )
## 
## # saving the scores
## saveRDS(object = p_values_df, file = "results/mne_permutation_decoding_scores.rds")

# importing the decoding scores
p_values_df <- readRDS(file = "results/mne_permutation_decoding_scores.rds")

# plotting the s-values (-log2(p-values))
p_values_df %>%
    mutate(sval = -log2(pval) ) %>%
    ggplot(aes(x = time, y = sval) ) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = -log2(0.05), linetype = 2) +
    labs(
        x = "Time (s)",
        y = "-log2(p-value)"
    )
