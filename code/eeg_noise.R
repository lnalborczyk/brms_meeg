##############################################################################
# Generate noise with the power spectrum of human EEG
# R version + control variance: GAR, University of Glasgow, Oct 2023
# From https://github.com/GRousselet/onsetsim/blob/main/code/eeg_noise.R
# Implemented by: Rafal Bogacz and Nick Yeung, Princeton Univesity, December 2002
# Adapted from Matlab code:
# https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator
############################################################################

##############################################################################
# Input:
#   frames - number of signal frames per each trial
#   epochs - number of simulated trials
#   srate - sampling rate of simulated signal
# Output:
#   signal - simulated EEG signal; vector: 1 by frames*epochs containing concatenated trials
##############################################################################

eeg_noise <- function (frames = 51, srate = 100, outvar = 1, meanpower) {
    
    sumsig <- 50 # number of sinusoids from which each simulated signal is composed of
    out <- vector(mode = "numeric", length = frames)
    freq <- 0
    
    for (i in 1:sumsig) {
        
        freq <- freq + (4 * runif(1) )
        freqamp <- meanpower[min(ceiling(freq), 125)] / meanpower[1]
        phase <- runif(1) * 2 * pi
        out <- out + sin((1:frames)/srate * 2 * pi * freq + phase) * freqamp
        
    }
    
    out <- out - mean(out)
    out <- out / sd(out)
    out <- out * sqrt(outvar)
    
    return(out)
    
}
