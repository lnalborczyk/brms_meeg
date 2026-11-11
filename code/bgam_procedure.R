###########################################################
# Illustrate the BGAMM procedure                          #
# ------------------------------------------------------- #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@cnrs.fr                    #
# Last updated on November 4, 2026                        #
###########################################################

library(tidybayes)
library(patchwork)
library(scales)

###################################################################
# 0. Plotting functions
###################################################################

plot_post_preds <- function (model, data = NULL, re_formula = NA) {
    
    if (is.null(data) ) data <- model$data
    
    # compute posterior estimate of the condition difference over time
    post_preds <- data %>%
        add_epred_draws(
            object = model,
            re_formula = re_formula
            ) %>%
        group_by(condition, time) %>%
        summarise(
            eeg_mean = median(.epred, na.rm = TRUE),
            lower = quantile(.epred, probs = 0.025, na.rm = TRUE),
            upper = quantile(.epred, probs = 0.975, na.rm = TRUE),
            .groups = "drop"
            )
    
    post_preds %>%
        ggplot(
            aes(
                x = time,
                y = eeg_mean,
                colour = condition,
                fill = condition
                )
            ) +
        geom_vline(xintercept = true_onset, linewidth = 0.5, linetype = 2) +
        geom_vline(xintercept = true_offset, linewidth = 0.5, linetype = 2) +
        geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2) +
        geom_point(
            data = model$data %>% summarise(
                eeg_mean = mean(eeg_mean),
                .by = c(time, condition)
                ),
            aes(fill = NULL),
            size = 1,
            pch = 21, show.legend = FALSE
            ) +
        geom_line(aes(y = eeg_mean), show.legend = FALSE) +
        geom_ribbon(
            aes(ymin = lower, ymax = upper, colour = NULL),
            alpha = 0.25, show.legend = FALSE
            ) +
        scale_colour_met_d(name = "Johnson") +
        scale_fill_met_d(name = "Johnson") +
        labs(x = "Time (s)", y = "EEG signal (a.u.)")
    
}

# define a function to plot the posterior estimate of the slope
plot_post_slope <- function (model, data = NULL, sesoi = NULL, re_formula = NA) {
    
    # if no data is specified, use model$data
    if (is.null(data) ) data <- model$data
    
    # compute posterior estimate of the condition difference over time
    post_prob_slope <- data %>%
        add_epred_draws(
            object = model,
            re_formula = re_formula
            ) %>%
        data.frame() %>%
        dplyr::select(participant, time, condition, .epred, .draw) %>%
        pivot_wider(names_from = condition, values_from = .epred) %>%
        mutate(epred_diff = cond2 - cond1) %>%
        group_by(time) %>%
        summarise(
            post_prob = median(epred_diff, na.rm = TRUE),
            lower = quantile(epred_diff, probs = 0.025, na.rm = TRUE),
            upper = quantile(epred_diff, probs = 0.975, na.rm = TRUE),
            .groups = "drop"
            )
    
    # define the ROPE / SESOI ribbon
    rope_df <- post_prob_slope %>%
        transmute(
            time = time,
            rope_low = -sesoi,
            rope_high = sesoi
            )
    
    # plot
    post_prob_slope %>%
        ggplot(aes(x = time, y = post_prob) ) +
        geom_vline(xintercept = true_onset, linewidth = 0.5, linetype = 2) +
        geom_vline(xintercept = true_offset, linewidth = 0.5, linetype = 2) +
        geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2) +
        
        # ROPE / SESOI ribbon
        geom_ribbon(
            data = rope_df,
            aes(x = time, ymin = rope_low, ymax = rope_high),
            inherit.aes = FALSE,
            fill = "orangered",
            alpha = 0.2
            ) +
        
        # 95% posterior credible interval
        geom_ribbon(
            aes(ymin = lower, ymax = upper),
            alpha = 0.2
            ) +
        
        # posterior median
        geom_line(linewidth = 0.8) +
        labs(
            x = "Time (s)",
            y = expression("Posterior slope estimate, " * beta(t) )
            )
    
}

###################################################################
# 1. Compute posterior odds clusters and boundary distributions
###################################################################

compute_boundary_posteriors <- function (
    model,
    sesoi = 0,
    threshold = 10,
    re_formula = NA
    ) {

    ###################################################################
    # Extract model information
    ###################################################################

    time_vec <- sort(unique(model$data$time) )
    condition_levels <- levels(model$data$condition)
    
    stopifnot(length(condition_levels) == 2)
    stopifnot(length(time_vec) >= 2)
    stopifnot(threshold > 1)
    stopifnot(is.finite(sesoi), sesoi >= 0)

    ###################################################################
    # Build group-level prediction grid
    ###################################################################

    newdata <- tidyr::crossing(
        time = time_vec,
        condition = condition_levels
        ) %>%
        mutate(
            condition = factor(
                condition,
                levels = condition_levels
                ),
            participant = NA
            )

    # Include SEM if required by the model
    if ("eeg_sem" %in% names(model$data) ) {
        
        newdata$eeg_sem <- 1
        
    }

    # Preserve original contrast coding
    contrasts(newdata$condition) <- contrasts(model$data$condition)

    ###################################################################
    # Extract joint posterior expected-value draws
    ###################################################################

    post_preds <- newdata %>%
        add_epred_draws(
            object = model,
            re_formula = re_formula
            ) %>%
        as_tibble() %>%
        select(.draw, time, condition, .epred) %>%
        mutate(condition = as.character(condition) ) %>%
        pivot_wider(
            names_from = condition,
            values_from = .epred
            ) %>%
        mutate(
            epred_diff =
                .data[[condition_levels[2]]] -
                .data[[condition_levels[1]]]
            ) %>%
        select(.draw, time, epred_diff) %>%
        arrange(.draw, time)

    n_draws <- n_distinct(post_preds$.draw)

    ###################################################################
    # METHOD 1: Pointwise posterior odds clusters
    ###################################################################

    prob_data <- post_preds %>%
        group_by(time) %>%
        summarise(
            m = mean(epred_diff > sesoi),
            .groups = "drop"
            ) %>%
        mutate(
            prob_ratio = m / (1 - m),
            # Avoid zero and infinite values on log scale
            prob_ratio = pmin(
                pmax(prob_ratio, 1 / n_draws),
                n_draws
                ),
            above_threshold = prob_ratio >= threshold
            ) %>%
        arrange(time) %>%
        mutate(
            cluster_id = cumsum(
                above_threshold !=
                    lag(above_threshold, default = FALSE)
                )
            )

    # Identify contiguous posterior-odds clusters
    cluster_data <- prob_data %>%
        filter(above_threshold) %>%
        group_by(cluster_id) %>%
        summarise(
            onset = min(time),
            offset = max(time),
            .groups = "drop"
            )

    ###################################################################
    # METHOD 2: Posterior draws onset and offset
    ###################################################################
    #
    # Completely independent of posterior odds clusters.
    #
    # For each draw:
    #   onset  = first time epred_diff > sesoi
    #   offset = last time epred_diff > sesoi
    #
    # Multiple excursions are allowed.
    # Draws with no exceedance receive NA boundaries.
    ###################################################################

    boundary_draws <- post_preds %>%
        group_by(.draw) %>%
        summarise(
            onset = if (any(epred_diff > sesoi) ) {
                min(time[epred_diff > sesoi])
            } else {
                NA_real_
            },

            offset = if (any(epred_diff > sesoi) ) {
                max(time[epred_diff > sesoi])
            } else {
                NA_real_
            },

            .groups = "drop"
        ) %>%
        mutate(
            effect_present = !is.na(onset)
        )

    ###################################################################
    # Summarise posterior boundary distributions
    ###################################################################

    summarise_boundary <- function (x) {

        x <- x[is.finite(x)]

        if (length(x) == 0) {
            return(c(
                median = NA_real_,
                q025 = NA_real_,
                q975 = NA_real_
            ))
        }

        c(
            median = median(x),
            q025 = unname(quantile(x, 0.025)),
            q975 = unname(quantile(x, 0.975))
        )
    }

    onset_summary <- summarise_boundary(boundary_draws$onset)
    offset_summary <- summarise_boundary(boundary_draws$offset)

    boundary_summary <- tibble(
        n_draws = nrow(boundary_draws),
        n_detected = sum(boundary_draws$effect_present),
        prob_present = mean(boundary_draws$effect_present),
        onset_median = unname(onset_summary["median"]),
        onset_q025 = unname(onset_summary["q025"]),
        onset_q975 = unname(onset_summary["q975"]),
        offset_median = unname(offset_summary["median"]),
        offset_q025 = unname(offset_summary["q025"]),
        offset_q975 = unname(offset_summary["q975"])
        )

    ###################################################################
    # Empirical posterior onset and offset densities
    ###################################################################

    boundary_long <- boundary_draws %>%
        pivot_longer(
            cols = c(onset, offset),
            names_to = "boundary",
            values_to = "time"
            ) %>%
        filter(!is.na(time) )

    # Kernel density estimates for plotting
    boundary_densities <- boundary_long %>%
        group_by(boundary) %>%
        group_split() %>%
        purrr::map_dfr(function (df) {

            if (n_distinct(df$time) < 2) {
                
                return (tibble() )
                
            }

            dens <- stats::density(
                df$time,
                n = 512,
                from = min(time_vec),
                to = max(time_vec)
                )

            tibble(
                boundary = df$boundary[1],
                time = dens$x,
                density = dens$y / max(dens$y)
                )
            
        })

    # --------------------------------------------------------
    # Return results
    # --------------------------------------------------------

    return (
        list(
            post_preds = post_preds,
            prob_data = prob_data,
            cluster_data = cluster_data,
            boundary_draws = boundary_draws,
            boundary_summary = boundary_summary,
            boundary_long = boundary_long,
            boundary_densities = boundary_densities,
            n_draws = n_draws,
            sesoi = sesoi,
            threshold = threshold
            )
        )
    
}

########################################################################
# 2. Plot posterior odds and empirical boundary distributions
########################################################################

plot_boundary_posteriors <- function (
    results,
    true_onset = NULL,
    true_offset = NULL,
    density_height = 1.5
    ) {

    n_draws <- results$n_draws
    threshold <- results$threshold
    density_base <- 1 / n_draws

    # Position densities near bottom of log10 odds scale
    # The vertical scaling is purely for visualisation.
    density_data <- results$boundary_densities %>%
        mutate(y_density = density_base * 10^(density_height * density) )

    p <- ggplot(
        results$prob_data,
        aes(x = time, y = prob_ratio)
        ) +

        ###################################################################
        # Empirical posterior distributions (independent)
        ###################################################################

        geom_ribbon(
            data = density_data,
            aes(
                x = time,
                ymin = density_base,
                ymax = y_density,
                group = boundary
                ),
            inherit.aes = FALSE,
            fill = "darkgreen",
            alpha = 0.2
            ) +

        geom_line(
            data = density_data,
            aes(
                x = time,
                y = y_density,
                group = boundary
                ),
            inherit.aes = FALSE,
            colour = "darkgreen",
            linewidth = 0.8
            ) +

        ###################################################################
        # Posterior odds
        ###################################################################

        geom_hline(
            yintercept = threshold,
            linewidth = 0.5,
            linetype = "dashed"
            ) +

        geom_line(linewidth = 0.8) +

        ###################################################################
        # Posterior odds clusters (separate method)
        ###################################################################

        geom_segment(
            data = results$cluster_data,
            aes(
                x = onset,
                xend = offset,
                y = threshold,
                yend = threshold
                ),
            inherit.aes = FALSE,
            colour = "darkgreen",
            linewidth = 2,
            lineend = "round"
            ) +

        ###################################################################
        # Axes and labels
        ###################################################################

        scale_y_log10(
            labels = label_log(digits = 2),
            limits = c(1 / n_draws, n_draws),
            breaks = 10^seq(
                floor(log10(1 / n_draws) ),
                ceiling(log10(n_draws) ),
                by = 1
                )
            ) +

        labs(
            x = "Time (s)",
            y = expression("Posterior odds for " * beta(t) > delta)
            )

    ###################################################################
    # Optional true onset/offset
    ###################################################################

    if (!is.null(true_onset) ) {
        
        p <- p + geom_vline(
            xintercept = true_onset,
            linewidth = 0.5,
            linetype = 2
            )
        
    }

    if (!is.null(true_offset) ) {
        
        p <- p + geom_vline(
            linewidth = 0.5,
            xintercept = true_offset,
            linetype = 2
            )
        
    }

    p
    
}

#############################################################
# 3. Settings
#############################################################

# Compute the SESOI as the baseline upper 95% quantile (model predictions)
# define baseline window
baseline_window <- c(0, 0.1)

# Build group-level newdata
time_vec <- sort(unique(meta_gam$data$time) )
condition_levels <- levels(meta_gam$data$condition)
reference_condition <- condition_levels[1]
target_condition <- condition_levels[2]

newdata_contrast <- tidyr::crossing(
    time = time_vec,
    condition = condition_levels
    ) %>%
    mutate(
        condition = factor(condition, levels = condition_levels),
        # dummy eeg_sem and participant variables to satisfy validate_data()
        # and exclude participant-speciic smooth when retrieving population
        # (group) level predictions
        eeg_sem = 1,
        participant = NA
        )

# Preserve contrast coding
contrasts(newdata_contrast$condition) <- contrasts(meta_gam$data$condition)

# Retrieve group-level posterior predictions
post_preds <- newdata_contrast %>%
    add_epred_draws(
        object = meta_gam,
        re_formula = NA
        ) %>%
    data.frame() %>%
    mutate(condition = as.character(condition) ) %>%
    select(.draw, time, condition, .epred) %>%
    pivot_wider(
        names_from = condition,
        values_from = .epred
        ) %>%
    mutate(epred_diff = .data[[target_condition]] - .data[[reference_condition]])

# Summarise posterior contrast during baseline
baseline_post_preds <- post_preds %>%
    filter(
        time >= baseline_window[1],
        time < baseline_window[2]
        ) %>%
    group_by(time) %>%
    summarise(
        epred_median = median(epred_diff, na.rm = TRUE),
        epred_q025 = quantile(epred_diff, 0.025, na.rm = TRUE),
        epred_q975 = quantile(epred_diff, 0.975, na.rm = TRUE),
        epred_abs_median = median(abs(epred_diff), na.rm = TRUE),
        .groups = "drop"
        )

###################################################################
# SESOI definitions
###################################################################

# Max positive group-level posterior median during baseline
meta_gam_sesoi_positive <- max(
    baseline_post_preds$epred_median,
    na.rm = TRUE
    )

# Max absolute group-level posterior median during baseline
meta_gam_sesoi_absolute <- max(
    abs(baseline_post_preds$epred_median),
    na.rm = TRUE
    )

# Upper 95% posterior quantile during baseline
meta_gam_sesoi_upper95 <- max(
    baseline_post_preds$epred_q975,
    na.rm = TRUE
    )

# Symmetric 95% baseline bound
meta_gam_sesoi_symmetric95 <- max(
    abs(baseline_post_preds$epred_q025),
    abs(baseline_post_preds$epred_q975),
    na.rm = TRUE
    )

# Simultaneous baseline-derived SESOI
# 95% percentile of the maximum absolute posterior contrast across baseline
# computed for each posterior draw
meta_gam_sesoi_simultaneous95 <- post_preds |>
    filter(
        time >= baseline_window[1],
        time < baseline_window[2]
        ) |>
    group_by(.draw) |>
    summarise(
        max_abs_baseline = max(abs(epred_diff) ),
        .groups = "drop"
        ) |>
    summarise(sesoi = quantile(max_abs_baseline, 0.95) ) |>
    pull(sesoi)

# Compare SESOI definitions
# c(
#     meta_gam_sesoi_positive, meta_gam_sesoi_absolute,
#     meta_gam_sesoi_upper95, meta_gam_sesoi_symmetric95, meta_gam_sesoi_simultaneous95
#     )

# SESOI based on baseline posterior median
# sesoi <- meta_gam_sesoi_absolute

# Alternative SESOI based on the point-wise 95% CrI
# sesoi <- meta_gam_sesoi_upper95
sesoi <- meta_gam_sesoi_symmetric95

# Alternative SESOI based on the simultaneous 95% CrI
# sesoi <- meta_gam_sesoi_simultaneous95

# Posterior odds threshold
threshold <- 3

###################################################################
# 4. Compute posterior results
###################################################################

results <- compute_boundary_posteriors(
    model = meta_gam,
    sesoi = sesoi,
    threshold = threshold,
    re_formula = NA
    )

# Posterior onset/offset summaries
# print(results$boundary_summary)

# Detected posterior odds clusters
# print(results$cluster_data)

###################################################################
# 5. Plot posterior predictions
###################################################################

p1 <- plot_post_preds(
    model = meta_gam,
    data = newdata_contrast
    )

###################################################################
# 6. Plot posterior slope
###################################################################

p2 <- plot_post_slope(
    model = meta_gam,
    sesoi = sesoi,
    data = newdata_contrast
    )

###################################################################
# 7. Plot posterior odds and boundary distributions
###################################################################

p3 <- plot_boundary_posteriors(
    results = results,
    true_onset = true_onset,
    true_offset = true_offset,
    density_height = 2
    )

###################################################################
# 8. Combine panels
###################################################################

combined_plot <- (
    free(p1, type = "label", side = "l") /
    free(p2, type = "label", side = "l") /
    free(p3, type = "label", side = "l")
    ) +
    plot_layout(axes = "collect_x") +
    plot_annotation(tag_levels = "A") # & theme(plot.tag = element_text(face = "bold") )

combined_plot

###################################################################
# 9. Save figure
###################################################################

ggsave(
    filename = "figures/bgam_procedure.png",
    plot = combined_plot,
    width = 9,
    height = 12,
    dpi = 300,
    # bg = "transparent"
    bg = "white"
    )
