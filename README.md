# Precise temporal localisation of M/EEG effects with Bayesian generalised additive multilevel models

## Abstract

Time-resolved electrophysiological measurements such as those offered by magneto- or electro-encephalography (M/EEG) provide a unique window onto neural activity underlying cognitive processes. Typically, researchers are interested in testing whether and when such measures differ across conditions and/or groups. The conventional approach consists in conducting mass-univariate statistics through time followed by some form of multiplicity correction (e.g., FDR, FWER) or cluster-based inference. However, these cluster-based methods have an important downside: they shift the focus of inference from the timepoint to the cluster level, thus preventing any conclusion to be made about the onset or offset of effects (e.g., differences across conditions or groups). Here, we introduce a *model-based* approch for analysing M/EEG timeseries such as ERPs or timecourses of decoding performance and their differences across conditions or groups. This approach relies on Bayesian generalised additive multilevel models, which output the posterior probabilility of the effect being above 0 (or above chance) at every timestep, while naturally taking into account the temporal dependencies and between-subject variability present in such data. Using both simulation and actual EEG data, we show that the proposed approach largely outperforms conventional methods in determining both the onset and offset of M/EEG effects (e.g., ERPs difference, decoding performance), producting more precise and more reliable estimates. We provide an R package implementing the approach and illustrate how to integrate it into M/EEG statistical pipelines in MNE-Python.

## Main simulation results

The figure below shows a summary of the simulation results, revealing that the proposed approach (`brms`) has the lowest MAE and variance for both the onset and offset estimates.

![Simulation results](figures/simulation_results_mae_variance.png)
