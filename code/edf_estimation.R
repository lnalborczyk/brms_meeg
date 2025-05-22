###########################################################
# Comparing various EDF estimates                         #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on May 22, 2025                            #
###########################################################

library(tidyverse)
library(bamlss)
library(brms)

# importing the reshaped MEG decoding results
decoding_df <- read.csv(file = "data/decoding_results_reshaped.csv") %>%
    # removing some participants
    mutate(participant_id = cur_group_id(), .by = participant) %>%
    filter(participant_id < 33) %>%
    select(-participant_id)

# averaging across participants
decoding_summary_df <- decoding_df %>%
    summarise(
        auc_mean = mean(auc),
        auc_sd = sd(auc),
        .by = c(participant, time)
        )

# setting possible k values
k_values <- seq(from = 10, to = 50, by = 5)

# initialising en empty dataframe to store the results
results <- data.frame()

# fitting the model for each k
for (k_value in k_values) {
    
    # printing progress
    cat("\nProcessing k =", k_value, "\n")
    
    # fitting the model with bamlss
    b_model <- bamlss(
        formula = auc_mean ~ s(time, bs = "tp", k = k_value),
        family = beta_bamlss,
        data = decoding_summary_df
        )
    
    # computing the EDF
    ss <- samplestats(samples = b_model)
    
    # constructing the smooth term dynamically
    smooth_term <- glue::glue("s(time, bs = 'tp', k = {k_value})")
    
    # full formula
    formula_str <- glue::glue("auc_mean ~ {smooth_term}")
    
    # convert to formula
    formula_obj <- brms::bf(formula_str)
    
    # fitting the Beta GAM
    meg_decoding_gam <- brm(
        formula = formula_obj,
        data = decoding_summary_df,
        family = brms::Beta(),
        warmup = 1000,
        iter = 4000,
        chains = 4,
        cores = 4
        )
    
    # computing the p_WAIC
    pwaic <- waic(meg_decoding_gam)
    looic <- loo(meg_decoding_gam)
    
    # combining the estimates in a dataframe
    temp_results <- data.frame(
        k = k_value,
        pd = ss$pd,
        pwaic = pwaic$estimates[2, 1],
        ploo = looic$estimates[2, 1]
        )
    
    # appending current results to previous results
    results <- bind_rows(results, temp_results)
    
}

# plotting the results
results %>%
    pivot_longer(cols = pd:ploo) %>%
    ggplot(aes(x = k, y = value, colour = name) ) +
    geom_abline(slope = 1, linetype = 2) +
    geom_line() +
    geom_point() +
    theme_bw() +
    labs(x = "User-defined k value", y = "EDF estimate", colour = "")

# plotting the correlation between PD and pWAIC
results %>%
    # pivot_longer(cols = pd:pwaic) %>%
    ggplot(aes(x = pd, y = pwaic) ) +
    # geom_abline(slope = 1) +
    # geom_line() +
    geom_point() +
    geom_smooth() +
    theme_bw()

# saving the plot
ggsave(
    # filename = "figures/edf.png",
    filename = "figures/meta_gam_ppc.png",
    width = 12, height = 9, dpi = 300,
    bg = "white",
    device = "png"
    )
