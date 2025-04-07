###########################################################
# Monte-Carlo simulation of onset/offset error properties #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on April 7, 2025                           #
###########################################################

# importing R packages
library(changepoint)
library(reticulate)
library(tidybayes)
library(dplyr)
library(tidyr)
library(brms)

# for setting up the cluster
library(doParallel)
library(foreach)

# total number of total CPU cores available on the HPC node
total_cores <- 12

# number of cores per brms model
cores_per_model <- 4

# number of parallel simulations to run at once
n_parallel <- total_cores / cores_per_model

# registering parallel backend
cl <- makeCluster(n_parallel)
registerDoParallel(cl)

# importing R version of Matlab code from Yeung et al. (2004)
source("code/eeg_noise.R")

# importing the ERP template with true onset = 160 ms, F = 81, and max at F = 126
source("code/erp_template.R")

# importing helper functions from Guillaume Rouselet
# source("code/functions.R")

# importing home-made helper functions
source("code/utils.R")

# to use with the eeg_noise function
meanpower <- unlist(read.table("code/meanpower.txt") )

# defining the total number of simulations to be performed
nsims <- 100

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

# defining possible (posterior probability ratio) threshold values
threshold_values <- c(1, 3, 6, 10, 20, 50, 100)

# initialising empty simulation results
sim_results <- data.frame()

# for each simulation
sim_results_list <- foreach(
    i = 1:nsims,
    .packages = c("brms", "dplyr", "tidyr", "tidybayes", "changepoint", "reticulate")
    ) %dopar% {
        
        # importing Python modules locally
        message_parallel("\nImporting Python modules...")
        use_condaenv("r-reticulate3", conda = "~/miniforge3/bin/conda", required = TRUE)
        
        # importing Python modules on the HPC cluster
        # py_require(packages = c("numpy", "mne"), python_version = "3.12.0")
        np <- import("numpy")
        mne <- import("mne")
        # mne$set_log_level("WARNING") # or "ERROR"
        
        # printing progress
        message_parallel(sprintf("\nRunning simulation %d out of %d...", i, nsims) )
        
        # simulating some EEG data
        raw_df <- generate_data(
            n_trials = n_trials, n_ppt = n_ppt, outvar = outvar,
            srate = srate, ronset = ronset
            )
        
        # summarising raw_data per participant
        ppt_data <- raw_df |>
            summarise(eeg = mean(eeg), .by = c(participant, condition, time) ) |>
            pivot_wider(names_from = condition, values_from = eeg) |>
            mutate(eeg_diff = cond2 - cond1)
        
        # summarising raw_data per participant (mean + SD)
        summary_df <- raw_df |>
            summarise(
                eeg_mean = mean(eeg),
                eeg_sd = sd(eeg),
                .by = c(participant, condition, time)
                )
        
        # defining a contrast for condition
        contrasts(summary_df$condition) <- c(-0.5, 0.5)
        
        # fitting the GAMM
        meta_gam <- brm(
            # using by-participant SD of ERPs across trials
            eeg_mean | se(eeg_sd) ~
                condition + s(time, bs = "cr", k = 20, by = condition) +
                (1 | participant),
            data = summary_df,
            family = gaussian(),
            warmup = 2000,
            iter = 5000,
            chains = cores_per_model,
            cores = cores_per_model
            # backend = "cmdstanr",
            # stan_model_args = list(stanc_options = list("O1") )
            )
        
        # computing probability metrics with the current threshold value
        message_parallel("\nModel fitted...")
        prob_y_above_0 <- meta_gam$data |>
            add_epred_draws(object = meta_gam) |>
            data.frame() |>
            dplyr::select(participant, time, condition, .epred, .draw) |>
            pivot_wider(names_from = condition, values_from = .epred) |>
            mutate(epred_diff = cond2 - cond1) |>
            # computing mean posterior probability at the group level
            group_by(time) |>
            summarise(m = mean(epred_diff > 0) ) |>
            mutate(prob_ratio = m / (1 - m) ) |>
            ungroup() |>
            # ensuring there is no 0 or +Inf values
            mutate(prob_ratio = ifelse(is.infinite(prob_ratio), ndraws(meta_gam), prob_ratio) ) |>
            mutate(prob_ratio = ifelse(prob_ratio == 0, 1 / ndraws(meta_gam), prob_ratio) )
        
        # initialising empty dataframe
        temp_brms_results <- data.frame()
        
        # looping over different threshold values
        for (j in seq_len(length(threshold_values) ) ) {
            
            # retrieving current threshold value
            threshold <- threshold_values[j]
            
            # finding clusters
            onset_offset_brms <- find_clusters(
                data = prob_y_above_0 |> select(time, value = prob_ratio),
                threshold = threshold
                ) |>
                mutate(simulation_id = formatC(x = i, width = 3, flag = 0) ) |>
                select(simulation_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
                mutate(method = "brms") |>
                pivot_longer(cols = onset:offset) |>
                select(simulation_id, method, cluster_id, onset_offset = name, value) |>
                data.frame() |>
                mutate(threshold = threshold) |>
                select(simulation_id, threshold, method, cluster_id, onset_offset, value)
            
            # binding to previous results
            temp_brms_results <- bind_rows(temp_brms_results, onset_offset_brms)
            
        }
        
        # converting data to a matrix (for later use in MNE functions)
        data_matrix <- matrix(
            ppt_data$eeg_diff,
            ncol = length(unique(ppt_data$time) ),
            byrow = TRUE
            )
        
        # massive univariate t-tests
        tests_results <- ppt_data |>
            group_by(time) |>
            summarise(
                tval = t.test(x = eeg_diff, mu = 0)$statistic^2,
                pval = t.test(x = eeg_diff, mu = 0)$p.value
                ) |>
            mutate(
                pval_bh = p.adjust(p = pval, method = "BH"),
                pval_by = p.adjust(p = pval, method = "BY"),
                pval_holm = p.adjust(p = pval, method = "holm")
                ) |>
            ungroup()
        
        # initialising empty dataframe
        temp_sim_results <- data.frame()
        
        # for each significance (alpha) threshold
        for (alpha_level in significance_thresholds) {
            
            # using the binary segmentation method to identify onsets and offsets
            res <- cpt.meanvar(data = tests_results$tval, method = "BinSeg", Q = 2)
            
            onset_offset_cpt <- data.frame(
                simulation_id = formatC(x = i, width = 3, flag = 0),
                threshold = alpha_level,
                cluster_id = 1,
                onset = tests_results$time[res@cpts[1]],
                offset = tests_results$time[res@cpts[2]]
                ) |>
                mutate(method = "cpt") |>
                pivot_longer(cols = onset:offset) |>
                select(simulation_id, threshold, method, cluster_id, onset_offset = name, value) |>
                data.frame()
            
            # finding clusters
            onset_offset_raw_p <- find_clusters(
                data = tests_results |> mutate(pval = pval * (-1) ) |> select(time, value = pval),
                threshold = -alpha_level
                ) |>
                mutate(
                    simulation_id = formatC(x = i, width = 3, flag = 0),
                    threshold = alpha_level
                    ) |>
                select(simulation_id, threshold, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
                mutate(method = "raw_p") |>
                pivot_longer(cols = onset:offset) |>
                select(simulation_id, threshold, method, cluster_id, onset_offset = name, value) |>
                data.frame()
            
            onset_offset_pval_bh <- find_clusters(
                data = tests_results |> mutate(pval = pval_bh * (-1) ) |> select(time, value = pval),
                threshold = -alpha_level
                ) |>
                mutate(
                    simulation_id = formatC(x = i, width = 3, flag = 0),
                    threshold = alpha_level
                    ) |>
                select(simulation_id, threshold, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
                mutate(method = "pval_bh") |>
                pivot_longer(cols = onset:offset) |>
                select(simulation_id, threshold, method, cluster_id, onset_offset = name, value) |>
                data.frame()
            
            onset_offset_pval_by <- find_clusters(
                data = tests_results |> mutate(pval = pval_by * (-1) ) |> select(time, value = pval),
                threshold = -alpha_level
                ) |>
                mutate(
                    simulation_id = formatC(x = i, width = 3, flag = 0),
                    threshold = alpha_level
                    ) |>
                select(simulation_id, threshold, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
                mutate(method = "pval_by") |>
                pivot_longer(cols = onset:offset) |>
                select(simulation_id, threshold, method, cluster_id, onset_offset = name, value) |>
                data.frame()
            
            onset_offset_pval_holm <- find_clusters(
                data = tests_results |> mutate(pval = pval_holm * (-1) ) |> select(time, value = pval),
                threshold = -alpha_level
                ) |>
                mutate(
                    simulation_id = formatC(x = i, width = 3, flag = 0),
                    threshold = alpha_level
                    ) |>
                select(simulation_id, threshold, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
                mutate(method = "pval_holm") |>
                pivot_longer(cols = onset:offset) |>
                select(simulation_id, threshold, method, cluster_id, onset_offset = name, value) |>
                data.frame()
            
            # running the MNE cluster-based permutation
            p_values_cluster_mass <- freq_stats_cluster_matrix(
                X = data_matrix,
                timesteps = unique(ppt_data$time),
                cluster_type = "mass",
                alpha_level = alpha_level
                )
            p_values_tfce <- freq_stats_cluster_matrix(
                X = data_matrix,
                timesteps = unique(ppt_data$time),
                cluster_type = "tfce",
                alpha_level = alpha_level
                )
            
            onset_offset_cluster_mass <- find_clusters(
                data = p_values_cluster_mass |> mutate(pval = pval*(-1) ) |> select(time, value = pval),
                threshold = -alpha_level
                ) |>
                mutate(
                    simulation_id = formatC(x = i, width = 3, flag = 0),
                    threshold = alpha_level
                    ) |>
                select(simulation_id, threshold, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
                mutate(method = "cluster_mass") |>
                pivot_longer(cols = onset:offset) |>
                select(simulation_id, threshold, method, cluster_id, onset_offset = name, value) |>
                data.frame()
            
            onset_offset_cluster_tfce <- find_clusters(
                data = p_values_tfce |> mutate(pval = pval*(-1) ) |> select(time, value = pval),
                threshold = -alpha_level
                ) |>
                mutate(
                    simulation_id = formatC(x = i, width = 3, flag = 0),
                    threshold = alpha_level
                    ) |>
                select(simulation_id, threshold, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
                mutate(method = "cluster_tfce") |>
                pivot_longer(cols = onset:offset) |>
                select(simulation_id, threshold, method, cluster_id, onset_offset = name, value) |>
                data.frame()
            
            # binding all results together
            alpha_temp_sim_results <- bind_rows(
                onset_offset_raw_p, onset_offset_pval_bh,
                onset_offset_pval_by, onset_offset_pval_holm, onset_offset_cpt,
                onset_offset_cluster_mass, onset_offset_cluster_tfce
                )
            
            # appending to previous results
            temp_sim_results <- bind_rows(temp_sim_results, alpha_temp_sim_results)
            
        }
        
        # bind these results with the brms results
        temp_sim_results <- bind_rows(temp_brms_results, temp_sim_results)
        
        # returning the results
        return (temp_sim_results)
    
    }

# binding all simulation results
sim_results <- bind_rows(sim_results_list)

# saving the simulation results
# saveRDS(object = sim_results, file = "./results/errors_results_cluster.rds")
saveRDS(object = sim_results, file = "./results/errors_results_100sims.rds")

# stopping the cluster
stopCluster(cl)
