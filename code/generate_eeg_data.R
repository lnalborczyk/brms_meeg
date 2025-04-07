# importing R version of Matlab code from Yeung et al. (2004)
source("code/eeg_noise.R")

# importing the ERP template with true onset = 160 ms, F=81, and max at F=126
source("code/erp_template.R")

# importing helper functions
# source("code/functions.R")

# to use with the eeg_noise function
meanpower <- unlist(read.table("code/meanpower.txt") )

# defining the number of simulations to be performed
# nsims <- 100

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
# significance_thresholds <- c(0.05, 0.01, 0.005, 0.001)

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
