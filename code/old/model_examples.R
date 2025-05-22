# averaging across participants
# ppt_df <- raw_df %>%
#     group_by(participant, condition, time) %>%
#     summarise(eeg = mean(eeg) ) %>%
#     ungroup()

# defining a contrast for condition
# contrasts(ppt_df$condition) <- c(-0.5, 0.5)

# fitting the BGAM
# gam <- brm(
#     # thin-plate regression splines with k-1 basis functions
#     eeg ~ condition + s(time, bs = "tp", k = 20, by = condition),
#     data = ppt_df,
#     family = gaussian(),
#     warmup = 2000,
#     iter = 5000,
#     chains = 8,
#     cores = 8,
#     file = "models/gam.rds"
#     )
# 
# meta_gam <- brm(
#     # using by-participant SD of ERPs across trials
#     eeg_mean | se(eeg_sd) ~
#         condition + s(time, bs = "tp", k = 20, by = condition) +
#         (1 + condition | participant),
#     data = summary_df,
#     family = gaussian(),
#     warmup = 2000,
#     iter = 5000,
#     chains = 8,
#     cores = 8,
#     file = "models/meta_gam.rds"
#     )
# 
# full_meta_gamm <- brm(
#     # using by-participant SD of ERPs across trials
#     eeg_mean | se(eeg_sd) ~ 1 + condition +
#         (1 + condition | participant) +
#         s(time, by = condition, bs = "tp", k = 20) +
#         # s(participant, by = condition, bs = "re", xt = list(bs = "tp") ),
#         s(time, participant, by = condition, bs = "fs", xt = list(bs = "tp") ),
#     data = summary_df,
#     family = gaussian(),
#     warmup = 2000,
#     iter = 5000,
#     chains = 8,
#     cores = 8,
#     control = list(adapt_delta = 0.95),
#     file = "models/full_meta_gamm.rds"
#     )

# errors on the HPC cluster
# Error in unserialize(socklist[[n]]) : error reading from connection
# Calls: %dopar% ... recvOneData -> recvOneData.SOCKcluster -> unserialize
# Execution halted
# slurmstepd: error: Detected 4 oom-kill event(s) in StepId=9052260.batch. Some of your processes may have been killed by the cgroup out-of-memory handler.

# Error in unserialize(socklist[[n]]) : error reading from connection
# Calls: %dopar% ... recvOneData -> recvOneData.SOCKcluster -> unserialize
# Execution halted
# slurmstepd: error: Detected 35 oom-kill event(s) in StepId=9060399.batch. Some of your processes may have been killed by the cgroup out-of-memory handler.

#!/bin/sh
#SBATCH -J reliability_simulation
#SBATCH -p skylake
#SBATCH -N 8
#SBATCH -n 8
#SBATCH --mem-per-cpu=12G # requesting 12GB per CPU core
#SBATCH -A b429
#SBATCH -t 6-12
#SBATCH -o /scratch/lnalborczyk/%x_%j.out
#SBATCH -e /scratch/lnalborczyk/%x_%j.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=ladislas.nalborczyk@cnrs.fr

# loading modules
# module purge
# module load userspace/all
# module load R/4.4.3
# module load python3/3.12.0

# fixing messages encoding
# export LC_ALL=en_US.UTF-8
# export LANG=en_US.UTF-8

# moving on scratch (working directory)
# cd /scratch/$SLURM_JOB_USER/
    
# sanity checks
# echo "Working directory: $(pwd)"
# echo "Running on: $SLURM_NODELIST"
# echo "Python version: $(which python)"

# executing the R script
# Rscript split_half_reliability_cluster.R
