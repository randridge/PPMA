## R Code for simulation studies in "A Sensitivity Analysis Framework Using the Proxy Pattern-Mixture Model for Generalization of Experimental Results"

Code files, in order, do the following:

1. Determine values of simulation parameters and save as a CSV (01_simulation_parms.R)
2. Generate simulated data (02_simulation_generate_data.R)
3. Perform ML estimation for RCT-PPMM using simulated data (03_simulation_runMethods_PPMM_MLE.R)
4. Perform Bayes estimation for RCT-PPMM using simulated data (04_simulation_runMethods_PPMM_Bayes.R)
5. Summarize simulation results (05_simulation_summarizeResults.R)
6. Create various graphical summaries (06_simulation_plots_trial_om.Rmd, 07_simulation_plots_PPMM_MLE.Rmd, 08_simulation_plots_PPMM_proxystrength.Rmd)
7. Calculate empirical coverage using simulated data (09_simulation_PPMM_coverage.Rmd)

Note that data from the yoga for breast cancer survivors RCT are **not** publicly available and cannot be shared.
