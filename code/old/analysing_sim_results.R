###########################################################
# Monte-Carlo simulation of onset/offset error properties #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on February 13, 2025                       #
###########################################################

library(patchwork)
library(tidyverse)
library(ggrepel)

# loading simulation results
sim_results <- readRDS(file = "results/sim_results.rds")

# identifying the best brms results
# eps=0 and threshold = 10 or 20 seems pretty good...
sim_results %>%
    select(simulation_id:offset_brms) %>%
    summarise(
        MAE_onset = median(abs(onset_brms - true_onset) ),
        variance_onset = var(onset_brms),
        MAE_offset = median(abs(offset_brms - true_offset) ),
        variance_offset = var(offset_brms),
        .by = c(eps, threshold)
        ) %>%
    arrange(MAE_onset, MAE_offset) %>%
    head(10)

# plotting some results
sim_results %>%
    summarise(error_brms = mean(abs(onset_brms-true_onset) ), .by = c(eps, threshold) ) %>%
    ggplot(aes(x = eps, y = threshold, fill = error_brms) ) +
    geom_tile(show.legend = FALSE) +
    geom_point(
        data = sim_results %>%
            summarise(error_brms = mean(abs(onset_brms-true_onset) ), .by = c(eps, threshold) ) %>%
            dplyr::filter(error_brms == min(error_brms, na.rm = TRUE) ),
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

sim_results %>%
    summarise(error_brms = mean(abs(offset_brms-true_offset) ), .by = c(eps, threshold) ) %>%
    # arrange(error_brms) %>% head()
    ggplot(aes(x = eps, y = threshold, fill = error_brms) ) +
    geom_tile(show.legend = FALSE) +
    geom_point(
        data = sim_results %>%
            summarise(error_brms = mean(abs(offset_brms-true_offset) ), .by = c(eps, threshold) ) %>%
            dplyr::filter(error_brms == min(error_brms, na.rm = TRUE) ),
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

# computing and plotting bias, MAE, and variance
# bias := true_onset - median(sim_onset_distrib)
# MAE :=  mean(abs(sim_onset-true_onset))
# variance := variance(sim_onset_distrib)
sim_results %>%
    # keeping only the best brms results
    filter(eps == 0 & threshold == 20) %>%
    select(-eps,-threshold) %>%
    pivot_longer(cols = -simulation_id, names_to = "measure", values_to = "value") %>%
    data.frame() %>%
    mutate(
        error = case_when(
            grepl("onset", measure) ~ abs(value - true_onset),
            grepl("offset", measure) ~ abs(value - true_offset)
            )
        ) %>%
    separate(measure, into = c("type", "method"), sep = "_", extra = "merge") %>%
    group_by(type, method) %>%
    summarise(
        MAE = median(error, na.rm = TRUE),
        variance = var(error, na.rm = TRUE),
        .groups = "drop"
        ) %>%
    ungroup() %>%
    mutate(type = factor(x = type, levels = c("onset", "offset") ) ) %>%
    mutate(
        method = factor(
            x = method,
            levels = c("p", "bh", "by", "holm", "cluster_mass", "cluster_tfce", "cpt", "brms"),
            labels = c("p-value", "FDR BH95", "FDR BY01", "Holm", "Cluster-mass", "TFCE", "Change point", "brms"),
            )
        ) %>%
    ggplot(
        aes(
            x = MAE,
            y = variance,
            colour = method,
            fill = method
            )
        ) +
    geom_point(
        size = 3, pch = 21,
        colour = "white",
        show.legend = FALSE
        ) +
    geom_label_repel(
        aes(label = method),
        colour = "white",
        show.legend = FALSE
        ) +
    scale_y_log10() +
    facet_wrap(~type) +
    theme_light(base_size = 12, base_family = "Open Sans") +
    scale_fill_manual(values = met.brewer(name = "Johnson", n = 8) ) +
    scale_colour_manual(values = met.brewer(name = "Johnson", n = 8) ) +
    labs(x = "Median absolute error (s)", y = "Variance (log-scale)")

# saving the plot
ggsave(
    filename = "figures/simulation_results_mae_variance.png",
    width = 12, height = 6, dpi = 300,
    device = "png"
    )

# plotting the distributions of onsets and offsets
sim_results %>%
    # keeping only the best brms results
    filter(eps == 0 & threshold == 20) %>%
    select(-eps,-threshold) %>%
    pivot_longer(
        cols = -simulation_id, 
        names_to = "measure",
        values_to = "value"
        ) %>%
    separate(measure, into = c("type", "method"), sep = "_", extra = "merge") %>%
    data.frame() %>%
    mutate(type = factor(x = type, levels = c("onset", "offset") ) ) %>%
    mutate(
        method = factor(
            x = method,
            levels = c("p", "bh", "by", "holm", "cluster_mass", "cluster_tfce", "cpt", "brms"),
            labels = c("p-value", "FDR BH95", "FDR BY01", "Holm", "Cluster-mass", "TFCE", "Change point", "brms"),
            )
        ) %>%
    mutate(plot_row = ifelse(method %in% c("p-value", "FDR BH95", "FDR BY01", "Holm"), 0, 1) ) %>%
    ggplot(
        aes(
            x = value,
            colour = method,
            fill = method
            )
        ) +
    geom_vline(xintercept = true_onset) +
    geom_vline(xintercept = true_offset) +
    geom_density(aes(fill = NULL) ) +
    # facet_wrap(~type, scales = "free_x") +
    # facet_wrap(~type) +
    facet_grid(plot_row~type, scales = "free_x") +
    theme_light(base_size = 12, base_family = "Open Sans") +
    scale_fill_manual(values = met.brewer(name = "Johnson", n = 8) ) +
    scale_colour_manual(values = met.brewer(name = "Johnson", n = 8) ) +
    labs(x = "Onset/Offset (s)", y = "Density")

# saving the plot
ggsave(
    filename = "figures/simulation_results_distributions.png",
    width = 12, height = 6, dpi = 300,
    device = "png"
    )
