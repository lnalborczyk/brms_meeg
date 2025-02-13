# Modelling M/EEG data with Bayesian generalised additive multilevel models

## Abstract

Time-resolved electrophysiological measurements such as those offered by magneto- or electro-encephalography (M/EEG) provide a unique window onto neural activity underlying cognitive process and how they unfold over time. Typically, we are interested in testing whether such measures differ across conditions and/or groups. The conventional approach consists in conducting mass-univariate statistics through followed by some form of multiplicity correction (e.g., FDR, FWER) or cluster-based inference. However, these cluster-based methods have an important downside: they shift the focus of inference from the timepoint to the cluster level, thus preventing any conclusion to be made about the onset and offset of effects (e.g., differences across conditions). Here, we introduce a novel *model-based approch* for analysing one-dimensional M/EEG timeseries such as ERPs or decoding timecourses and their differences across conditions or group. This approach relies on Bayesian generalised additive multilevel models, which output the posterior probabilility of the effect being above 0 (or above chance) at every timestep, while naturally taking into account the temporal dependencies and between-subject variability present in such data.

## Main simulation results

The figure below shows a summary of the simulation results, revealing that the proposed approach (`brms`) has the lowest MAE and variance for both the onset and offset estimates.

![Simulation results](figures/simulation_results_mae_variance.png)
