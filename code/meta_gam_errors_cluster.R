###########################################################
# Monte-Carlo simulation of onset/offset error properties #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on April 7, 2025                           #
###########################################################

# importing R packages
library(tidybayes)
library(dplyr)
library(tidyr)
library(brms)

# setting up parallel workers
library(purrr)
library(furrr)

# total number of total CPU cores available on the HPC node
total_cores <- 12

# number of cores per brms model
cores_per_model <- 4

# number of parallel simulations to run at once
n_parallel <- total_cores / cores_per_model

# setting up parallel workers
future::plan(future::multisession(workers = n_parallel) )

# number of simulations to run per combination
nsims <- 100

# importing the data generation functions 
source(file = "code/generate_eeg_data.R")
source(file = "code/utils.R")

# defining a function to run a single iteration of the simulation
single_iteration <- function (i) {
    
    # sanity check
    system(
        sprintf(
            "echo 'Running simulation %d on PID %d with %s cores' >> ./sim_log.txt",
            i, Sys.getpid(), cores_per_model
            )
        )
    
    # initialising an empty dataframe
    temp_results <- data.frame()
    
    # number of posterior samples to use to compute the posterior probability ratio
    n_post_samples <- 1e4
    
    # defining the parameter grid
    param_grid <- crossing(
        kvalue = c(5, 10, 15, 20, 25, 30, 35, 40),
        threshold = seq(from = 1, to = 50, by = 1)
        )
    
    for (j in seq_len(nrow(param_grid) ) ) {
        
        # retrieving the current kvalue and threshold
        kvalue <- param_grid$kvalue[j]
        threshold <- param_grid$threshold[j]
        
        # printing progress
        message_parallel(sprintf("Running sim %d, k = %d, threshold = %d", i, kvalue, threshold) )
        
        # when we assess a new sim_id, we generate new data
        if (kvalue == 5 & threshold == 1) {
    
            # generate some data
            message_parallel("\nSimulating new data...\n")
            raw_df <- generate_data(
                n_trials = n_trials, n_ppt = n_ppt, outvar = outvar,
                srate = srate, ronset = ronset
                )
    
            # averaging across participants
            summary_df <- raw_df %>%
                summarise(
                    eeg_mean = mean(eeg),
                    eeg_sd = sd(eeg),
                    .by = c(participant, condition, time)
                    )
    
            # defining a contrast for condition
            contrasts(summary_df$condition) <- c(-0.5, 0.5)
            
        }
        
        # for each new k value, we fit a new GAM
        if (threshold == 1) {
            
            # construct the smooth term dynamically
            smooth_term <- glue::glue("s(time, bs = 'tp', k = {kvalue}, by = condition)")
            
            # full formula
            formula_str <- glue::glue("eeg_mean | se(eeg_sd) ~ condition + {smooth_term} + (1 | participant)")
            
            # convert to formula
            formula_obj <- brms::bf(formula_str)
            
            # fitting the GAM
            message_parallel("\nFitting the GAM...\n")
            error_model <- brm(
                formula = formula_obj,
                data = summary_df,
                family = gaussian(),
                warmup = 2000,
                iter = 5000,
                chains = cores_per_model,
                cores = cores_per_model
                # backend = "cmdstanr",
                # stan_model_args = list(stanc_options = list("O1") )
                )
            
            prob_y_above_0 <- error_model$data %>%
                add_epred_draws(object = error_model, ndraws = n_post_samples) %>%
                # add_epred_draws(object = error_model) %>%
                data.frame() %>%
                dplyr::select(participant, time, condition, .epred, .draw) %>%
                pivot_wider(names_from = condition, values_from = .epred) %>%
                mutate(epred_diff = cond2 - cond1) %>%
                # computing mean posterior probability at the group level
                group_by(time) %>%
                summarise(m = mean(epred_diff > 0) ) %>%
                mutate(prob_ratio = m / (1 - m) ) %>%
                ungroup() %>%
                # ensuring there is no 0 or +Inf values
                mutate(
                    prob_ratio = ifelse(
                        is.infinite(prob_ratio),
                        # ndraws(error_model),
                        n_post_samples,
                        prob_ratio
                        )
                    ) %>%
                mutate(
                    prob_ratio = ifelse(
                        prob_ratio == 0,
                        # 1 / ndraws(error_model),
                        1 / n_post_samples,
                        prob_ratio
                        )
                    )
        
        }
        
        # finding clusters
        clusters <- find_clusters(
            data = prob_y_above_0 |> select(time, value = prob_ratio),
            threshold = threshold
            ) |>
            mutate(
                kvalue = kvalue,
                threshold = threshold,
                sim_id = formatC(x = i, width = 3, flag = 0),
                ) |>
            select(kvalue, threshold, sim_id, cluster_id, cluster_onset, cluster_offset)
        
        # binding to previous results
        temp_results <- bind_rows(temp_results, clusters)
        
    }
    
    # returning the clusters
    return (temp_results)

}

# running the simulation n_sims times (in parallel)
results <- future_map_dfr(
    .x = 1:nsims,
    .f = single_iteration,
    .progress = TRUE
    )

# saving the simulation results
# saveRDS(object = results, file = "./results/meta_gam_error_properties_cluster.rds")
saveRDS(object = results, file = "./results/meta_gam_error_properties_100sims.rds")

# resetting back to sequential and releasing parallel workers
plan(sequential)
