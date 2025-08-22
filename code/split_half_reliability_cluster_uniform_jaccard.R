###########################################################
# Split-half reliability of onset/offset estimates        #
# Written by Ladislas Nalborczyk                          #
# Contact: ladislas.nalborczyk@gmail.com                  #
# Last updated on August 7, 2025                          #
###########################################################

# logging error messages on the HPC cluster
# options(error = function () {traceback(2); quit("no", status = 1, runLast = FALSE)})

# importing R packages
library(changepoint)
library(reticulate)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(glue)
library(brms)

# for setting up the cluster
library(doParallel)
library(foreach)

# retrieving the array ID
# ensuring this is passed in the SLURM job array script: #SBATCH --array=1-50
args <- commandArgs(trailingOnly = TRUE)
# array_id <- as.numeric(args[1])
array_id_increment <- 200 # id increment for remaining simulations
array_id <- array_id_increment + as.numeric(args[1])
cat("Running job for array_id =", array_id, "\n")

# total number of total CPU cores available on the HPC node
total_cores <- 16

# number of cores per brms model
cores_per_model <- 4

# number of parallel simulations to run at once
n_parallel <- total_cores / cores_per_model

# registering parallel backend
cl <- makeCluster(n_parallel)
registerDoParallel(cl)

# importing home-made helper functions
source("code/utils.R")

# importing the reshaped MEG decoding results
raw_df <- read.csv(file = "data/decoding_results_reshaped.csv") |>
    # removing the last participant (to get a round number of participants)
    mutate(participant_id = cur_group_id(), .by = participant) |>
    dplyr::filter(participant_id < 33) |>
    select(-participant_id)

# averaging across participants
summary_df <- raw_df |>
    summarise(
        auc_mean = mean(auc),
        auc_sd = sd(auc),
        auc_se = sd(auc) / sqrt(n() ),
        .by = c(participant, time)
        )

# defining the decoding chance level
chance_level <- 0.5

# importing the data splits list
# split_df <- read.csv(file = "data/uniform_jaccard_data_splits_100pairs.csv")
# split_df <- read.csv(file = "data/uniform_jaccard_data_splits_1000pairs.csv")
# split_df <- read.csv(file = "data/uniform_jaccard_data_splits_1000pairs_remaining.csv")
split_df <- read.csv(file = "data/uniform_jaccard_data_splits_1000pairs_remaining2.csv")

# total number of unique splits
all_splits <- unique(split_df$split_id)
N_total_splits <- length(all_splits)

# total number of array jobs
N_array_jobs <- 10

# computing the range of splits for this job
splits_per_job <- ceiling(N_total_splits / N_array_jobs)
# split_start <- (array_id - 1) * splits_per_job + 1
# split_end <- min(array_id * splits_per_job, N_total_splits)
split_start <- (array_id - array_id_increment - 1) * splits_per_job + 1
split_end <- min((array_id - array_id_increment) * splits_per_job, N_total_splits)

# retrieving the splits' IDs assigned to this job
assigned_splits <- all_splits[split_start:split_end]

# sanity check
cat("Processing split_numbers:", paste(range(assigned_splits), collapse = " - "), "\n")

# subsetting split_df accordingly
split_df_this_node <- split_df |> dplyr::filter(split_id %in% assigned_splits)

# number of simulations to run on this node (one per split)
n_splits <- length(unique(split_df_this_node$split_id) )

# creating a dataframe with all (split_number, group) combinations to iterate over
split_group_combinations <- crossing(
    split_number = assigned_splits,
    group = c("A", "B")
    ) |>
    # creating a custom simulation ID per group
    mutate(split_id = sprintf("a%03d_s%05d_%s", array_id, split_number, group) )

# defining the threshold on the posterior probability ratio
post_prob_ratio_threshold <- 20

# defining the significance level
alpha_level <- 0.05

# for each simulation
reliability_results <- foreach(
    i = 1:nrow(split_group_combinations),
    .combine = bind_rows,
    .errorhandling = "pass",
    .packages = c("reticulate", "brms", "dplyr", "tidyr", "changepoint", "tidybayes", "ggplot2")
    ) %dopar% {

	tryCatch({
	    
        # importing Python modules locally
        # message_parallel("\nImporting Python modules...")
        # use_condaenv("r-reticulate3", conda = "~/miniforge3/bin/conda", required = TRUE)
        
        # importing Python modules on the HPC cluster
        py_require(packages = c("numpy==2.2.0", "pillow==11.2.1", "mne"), python_version = "3.12.0")
        np <- import("numpy")
        mne <- import("mne")
        mne$set_log_level("WARNING") # or "ERROR"
        
        # extracting split_number and group for this simulation
        current_split <- split_group_combinations$split_number[i]
        current_group <- split_group_combinations$group[i]
        current_split_id <- split_group_combinations$split_id[i]
        
        # printing progress
        message_parallel(sprintf("\nRunning simulation: %s (%d out of %d)...", current_split_id, i, nrow(split_group_combinations) ) )
        
        # selecting participants from the correct group
        participants_split <- split_df_this_node |>
            dplyr::filter(split_id == current_split, group == current_group) |>
            pull(participant)
        
        # subsetting data accordingly
        raw_df_split <- raw_df |> dplyr::filter(participant %in% participants_split)
        summary_df_split <- summary_df |> dplyr::filter(participant %in% participants_split)
        
        # computing the smallest effect size of interest (SESOI)
        # as the SD of decoding performance during the baseline
        # sesoi <- summary_df_split |>
        #     dplyr::filter(time < 0) |>
        #     summarise(auc_mean = mean(auc_mean), .by = time) |>
        #     summarise(sesoi = sd(auc_mean) ) |>
        #     pull(sesoi)
        
        # computing the SESOI as the baseline upper 95% quantile
        sesoi <- summary_df_split |>
            filter(time < 0) |>
            summarise(auc_mean = mean(auc_mean), .by = time) |>
            summarise(quantile(x = auc_mean, probs = 0.90) ) |>
            pull() - chance_level
        
        # fitting the GAM
        gam_split <- brm(
            auc_mean ~ s(time, bs = "cr", k = 50),
            data = summary_df_split,
            # family = Beta(),
            family = gaussian(),
            # warmup = 1000,
            # iter = 5000,
            warmup = 2000,
            iter = 10000,
            chains = cores_per_model,
            cores = cores_per_model
            # backend = "cmdstanr",
            # stan_model_args = list(stanc_options = list("O1") )
            )
        
        # checking model's predictions against raw data
        # plot(
        #     conditional_effects(x = gam_split),
        #     line_args = list(colour = "steelblue", fill = "steelblue", alpha = 0.3),
        #     points = FALSE, plot = FALSE
        #     )[[1]] +
        #     geom_hline(yintercept = 0.5, linetype = 2) +
        #     geom_vline(xintercept = 0.0, linetype = 2) +
        #     annotate(
        #         geom = "rect",
        #         xmin = -Inf, xmax = Inf,
        #         ymin = chance_level - sesoi,
        #         ymax = chance_level + sesoi,
        #         fill ="orangered",
        #         alpha = 0.2
        #         ) +
        #     geom_line(
        #         data = summary_df_split |> summarise(auc_mean = mean(auc_mean), .by = time),
        #         aes(x = time, y = auc_mean), inherit.aes = FALSE
        #         ) +
        #     labs(x = "Time (s)", y = "Decoding accuracy (ROC AUC)")
        
        # saving the plot
        # ggsave(
        #     filename = paste0(
        #         "reliability_plots/gam_split_predictions_split_",
        #         current_split_id, ".png"
        #         ),
        #     width = 12, height = 6, dpi = 300,
        #     device = "png", bg = "white"
        #     )
        
        # computing the posterior probability
        message_parallel("\nModel fitted...")
        prob_y_above_0 <- data.frame(time = unique(gam_split$data$time) ) |>
            # using a subset of posterior samples
            # add_epred_draws(object = gam_split, ndraws = n_post_samples) |>
            # or using all posterior samples
            add_epred_draws(object = gam_split) |>
            # converting to dataframe
            data.frame() |>
            group_by(time) |>
            summarise(m = mean(.epred > (0 + chance_level + sesoi) ) ) |>
            # summarise(m = mean(.epred >= sesoi) ) |>
            mutate(prob_ratio = m / (1 - m) ) |>
            ungroup() |>
            # ensuring there is no 0 or +Inf values
            mutate(prob_ratio = pmin(prob_ratio, ndraws(gam_split) ) ) |>
            mutate(prob_ratio = pmax(prob_ratio, 1 / ndraws(gam_split) ) )
        
        # finding clusters
        onset_offset_brms <- find_clusters(
            data = prob_y_above_0 |> select(time, value = prob_ratio),
            threshold = post_prob_ratio_threshold
            ) |>
            mutate(split_id = current_split_id) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "brms") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        # converting data to a matrix (for later use in MNE functions)
        data_matrix <- matrix(
            summary_df_split$auc_mean,
            ncol = length(unique(summary_df_split$time) ),
            byrow = TRUE
            )
        
        # massive univariate t-tests
        tests_results <- summary_df_split|>
            group_by(time) |>
            summarise(
                tval = t.test(x = auc_mean, mu = chance_level)$statistic^2,
                pval = t.test(x = auc_mean, mu = chance_level)$p.value
                ) |>
            mutate(
                pval_bh = p.adjust(p = pval, method = "BH"),
                pval_by = p.adjust(p = pval, method = "BY"),
                pval_holm = p.adjust(p = pval, method = "holm")
                ) |>
            ungroup()
            
        # using the binary segmentation method to identify onsets and offsets
        res <- cpt.meanvar(data = tests_results$tval, method = "BinSeg", Q = 2)
        onset_offset_cpt <- data.frame(
            split_id = current_split_id,
            cluster_id = 1,
            onset = tests_results$time[res@cpts[1]],
            offset = tests_results$time[res@cpts[2]]
            ) |>
            mutate(method = "cpt") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        # finding clusters
        onset_offset_raw_p <- find_clusters(
            data = tests_results |> mutate(pval = pval * (-1) ) |> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = current_split_id) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "raw_p") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        onset_offset_pval_bh <- find_clusters(
            data = tests_results |> mutate(pval = pval_bh * (-1) ) |> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = current_split_id) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "pval_bh") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        onset_offset_pval_by <- find_clusters(
            data = tests_results |> mutate(pval = pval_by * (-1) ) |> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = current_split_id) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "pval_by") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        onset_offset_pval_holm <- find_clusters(
            data = tests_results |> mutate(pval = pval_holm * (-1) ) |> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = current_split_id) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "pval_holm") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        # running the MNE cluster-based permutation
        p_values_cluster_mass <- freq_stats_cluster_matrix(
            X = data_matrix - chance_level,
            timesteps = unique(summary_df_split$time),
            cluster_type = "mass",
            alpha_level = alpha_level
            )
        p_values_tfce <- freq_stats_cluster_matrix(
            X = data_matrix - chance_level,
            timesteps = unique(summary_df_split$time),
            cluster_type = "tfce",
            alpha_level = alpha_level
            )
        
        onset_offset_cluster_mass <- find_clusters(
            data = p_values_cluster_mass |> mutate(pval = pval*(-1) ) |> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = current_split_id) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "cluster_mass") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        onset_offset_cluster_tfce <- find_clusters(
            data = p_values_tfce |> mutate(pval = pval*(-1) ) |> select(time, value = pval),
            threshold = -alpha_level
            ) |>
            mutate(split_id = current_split_id) |>
            select(split_id, cluster_id, onset = cluster_onset, offset = cluster_offset) |>
            mutate(method = "cluster_tfce") |>
            pivot_longer(cols = onset:offset) |>
            select(split_id, method, cluster_id, onset_offset = name, value) |>
            data.frame()
        
        # binding all results together
        temp_reliability_results <- bind_rows(
            onset_offset_brms, onset_offset_raw_p, onset_offset_pval_bh,
            onset_offset_pval_by, onset_offset_pval_holm, onset_offset_cpt,
            onset_offset_cluster_mass, onset_offset_cluster_tfce
            ) |>
            mutate(sesoi = sesoi)
        
        # cleaning up temporary files
        # unlink(file.path(Sys.getenv("TMPDIR"), "*"), recursive = TRUE)
        
        # cleaning up Stan temporary files only
        Sys.setenv(TMPDIR = "./tmp")
        tmpdir <- Sys.getenv("TMPDIR")
        message_parallel(paste("\nCleaning files in tmpdir:", tmpdir) )
        unlink(list.files(tmpdir, pattern = "\\.(cpp|o|hpp)$", full.names = TRUE), recursive = TRUE)
        
        # returning the results
        temp_reliability_results

	}, error = function(e) {
		
		# returning a row with NA and error message
		data.frame(
		    split_id = split_group_combinations$split_id[i],
		    error = paste("Error:", e$message),
		    stringAsFactors = FALSE
		    )

		 })
        
    }

# saving the results
saveRDS(object = reliability_results, file = glue::glue("./results/reliability_results_array_{array_id}.rds") )

# stopping the cluster
stopCluster(cl)
