library(reticulate)

# use_virtualenv("./reticulate", required = TRUE)
py_require(packages = c("numpy", "mne"), python_version = "3.12.0")
np <- import("numpy")
mne <- import("mne")
X <- np$random$randn(as.integer(10), as.integer(50) )
mne$set_log_level("WARNING") # or "ERROR"

results <- tryCatch({
    mne$stats$spatio_temporal_cluster_1samp_test(
        X = X,
        out_type = "mask",
        threshold = NULL,
        n_permutations = as.integer(2^12),
        n_jobs = as.integer(1),
        verbose = TRUE
        )
    }, error = function(e) {
    message("Error: ", e$message)
    reticulate::py_last_error()
    return(NULL)
    })

print(results)
