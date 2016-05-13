Code used in: Fei Deng, Mirko Ledda, Sana Vaziri and Sharon Aviran (2016) Data-directed RNA secondary structure prediction using probabilistic modeling.

## Installation
1. Install R or RStudio for your operating sytem.
2. Unzip the folder `script/data.zip`.
3. In R or RStudio, set the working directory to `Manuscript_code/` before running any scripts.

### Required R libraries:
Use `install.packages(foo)` to install these packages.

* ggplot2
* reshape2
* gridExtra
* grid
* nortest
* PearsonDS
* evd
* ks

## Description
### Folders
* `data/` Contains all data (needs to be unzip from `data.zip`)
* `data/ct_sequence_file/` .ct and .fa files for all reference RNAs
* `data/shape_raw/` .shape files for all reference RNAs
* `data/generated_data/` Data generated in this study - e.g. simulations
* `data/folding_results/` Folding results obtained from RNAstructure or RNAprob (original output files were parsed into .RData files).
* `output/` output folder - where figures/files are written
* `tmp/` Folder for temporary files
* `src/` Source code (contains custom functions)

### Scripts
#### Figure 1
* `plot_all_foldings.r` Create figure 1 and run pairwise $t$-tests
* `kde_based_shape_simulation.R` Simulate SHAPE data based on gaussian kernel density estimation (KDE)

#### Figure 2 / Figure S5
* `reassign_zeros_reactivities.r` Random reassignement of zero and negative SHAPE reactivities data to positive values using the distribution of negative SHAPE reactivities. Also generates Absolute, Ln or Box-Cox transformed SHAPE data. 
* `plot_SHAPE_data.r` Create figure 2. Run normality tests.

#### Figure 3 / Figure S4
* `prepare_noise_and_quantile_cross_validation_from_real.R` Prepare data, based on real data, for the cross-validation analysis as well as the noise simulation.
* `prepare_quantile_cross_validation_from_sim.R` Prepare data, based on simulated data, for the cross-validation analysis as well as the noise simulation.
* `plot_cross_validation.r` Create figure 3 and the quintile range file

#### Figure 4
* `plot_reconstruction.r` Create figure 4
* `real_quantile_reconstruction.R` Generate quantile reconstructed data (for figure 4) based on real SHAPE data
* `sim_quantile_reconstruction.R` Generate quantile reconstructed data (for figure 4) based on simulated SHAPE data

#### Figure 6
* `bayesian_prob_pairing_kde.r` Compute $P(\pi_i|\alpha_i)$ and create figure 6 (+ plots for each RNA).

#### Table 1 / Figure S7
* `mock_probes/simulate_mock_probe_sc1.r` Generate simulation data for mock probe scenario 1.
* `mock_probes/simulate_mock_probe_sc2.r` Generate simulation data for mock probe scenario 2.
* `mock_probes/plot_mock_probe.r` Create figure S7. 
* `mock_probes/mock_probe_stats.r` Stats for Table 1.
* `mock_probes/make_param_files_mock_probe.r` Fit mock probe simulations and creates the parameter files for RNAprob-3



