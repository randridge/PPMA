library(tidyverse)
library(rsimsum)

# data generation parameters for each set -----
load("./simdata/simParms.RData")
parms <- parms %>%
  arrange(simset) %>%
  mutate(simset = factor(simset))

# read in truth and trial estimates -----
filelist <- fs::dir_ls("./simresults/truth_trial/", regexp = "\\.csv$")
truth_trial <- filelist %>% 
  map(\(x) read_csv(x, show_col_types=FALSE)) %>%
  list_rbind() %>%
  add_column(estimator = "trial") %>%
  rename(estimate = trial) %>%
  dplyr::select(simset, rep, truth, estimator, estimate) %>%
  mutate(simset = factor(simset))

# Truth -----
truth <- truth_trial %>% 
  dplyr::select(simset, truth) %>%
  group_by(simset) %>%
  summarize(truth = mean(truth))
save(truth, file="./simresults/simSummary_truth.RData")

# Trial -----
sumstats_trial <- truth_trial %>% 
  simsum(estvarname = "estimate", 
         true = "truth", 
         methodvar = "estimator", 
         by = c("simset")) %>%
  summary() %>%
  tidy() %>%
  as_tibble()
save(sumstats_trial, file="./simresults/simSummary_trial.RData")

# Proxy strength (MLE) -----
filelist <- fs::dir_ls("./simresults/ppmm_mle/", regexp = "proxyinfo")
mle_proxystrength <- filelist %>% 
  map(\(x) read_csv(x, show_col_types=FALSE)) %>%
  list_rbind() %>%
  rename(proxy_A1 = rho_XY_s_A1,
         proxy_A0 = rho_XY_s_A0) %>%
  group_by(simset) %>%
  summarize(rho_XY_s_A1 = median(proxy_A1),
            rho_XY_s_A0 = median(proxy_A0),
            lb_rho_XY_s_A1 = quantile(proxy_A1, .025),
            ub_rho_XY_s_A1 = quantile(proxy_A1, .975),
            lb_rho_XY_s_A0 = quantile(proxy_A0, .025),
            ub_rho_XY_s_A0 = quantile(proxy_A0, .975))
save(mle_proxystrength, file="./simresults/simSummary_proxystrength_mle.RData")

# read in PPMM MLE results -----
filelist <- fs::dir_ls("./simresults/ppmm_mle/", regexp = "simResults")
results_mle <- filelist %>% 
  map(\(x) read_csv(x, show_col_types=FALSE)) %>%
  list_rbind() %>%
  arrange(simset) %>%
  mutate(estimator = case_when(phi_A1 == 0   & phi_A0 == 0.  ~ "trt_phi0_pla_phi0",
                               phi_A1 == 0.5 & phi_A0 == 0   ~ "trt_phi0.5_pla_phi0",
                               phi_A1 == 1   & phi_A0 == 0   ~ "trt_phi1_pla_phi0",
                               phi_A1 == 0   & phi_A0 == 0.5 ~ "trt_phi0_pla_phi0.5",
                               phi_A1 == 0   & phi_A0 == 1   ~ "trt_phi0_pla_phi1")) %>%
  rename(estimate = trt_effect_ns) %>%
  select(simset, rep, estimator, estimate) %>%
  mutate(simset = factor(simset))
# merge in truth
results_mle <- truth_trial %>% 
  dplyr::select(simset, rep, truth) %>%
  full_join(results_mle)

# PPMM MLE: calculate simulation summary stats for methods (point ests only) -----
sumstats_mle <- results_mle %>%
  simsum(estvarname = "estimate", 
         true = "truth", 
         methodvar = "estimator", 
         by = c("simset"), 
         ref = "trt_phi0_pla_phi0") %>%
  summary() %>%
  tidy() %>%
  as_tibble()
save(sumstats_mle, file="./simresults/simSummary_PPMM_MLE.RData")

# PPMM MLE: coverage -----
cover_mle <- results_mle %>%
  pivot_wider(values_from = "estimate", names_from = "estimator") %>%
  full_join(dplyr::select(truth_trial, c(simset, rep, truth))) %>%
  rename(bound1 = trt_phi0_pla_phi0) %>%
  pivot_longer(cols=c(trt_phi1_pla_phi0, trt_phi0_pla_phi1),
               names_to = "estimator", values_to = "bound2", names_prefix = "est_") %>%
  select(simset, rep, estimator, bound1, bound2, truth) %>%
  mutate(lb = pmin(bound1, bound2),
         ub = pmax(bound1, bound2)) %>%
  mutate(cover = (lb < truth & ub > truth)) %>%
  group_by(simset, estimator) %>%
  summarize(nsim = n(), 
            cover = mean(cover)) %>%
  ungroup() 
save(cover_mle, file="./simresults/simSummary_PPMM_MLE_cover.RData")

# read in Bayes results -----
filelist <- fs::dir_ls("./simresults/ppmm_bayes/", regexp = "\\.csv$")
results_bayes <- filelist %>%
  map(\(x) read_csv(x, show_col_types=FALSE)) %>%
  list_rbind()

# Bayes: proxy strength -----
bayes_proxystrength <- results_bayes %>%
  filter(phi1 == "ign" & phi0 == "ign") %>%
  select(simset, rep, starts_with("rho")) %>%
  group_by(simset) %>%
  summarize(rho_XY_s_A1 = median(rho_est_A1),
            rho_XY_s_A0 = median(rho_est_A0),
            lb_rho_XY_s_A1 = quantile(rho_est_A1, .025),
            ub_rho_XY_s_A1 = quantile(rho_est_A1, .975),
            lb_rho_XY_s_A0 = quantile(rho_est_A0, .025),
            ub_rho_XY_s_A0 = quantile(rho_est_A0, .975))
save(bayes_proxystrength, file="./simresults/simSummary_proxystrength_bayes.RData")

# Bayes: coverage -----
cover_bayes <- results_bayes %>%
  arrange(simset) %>%
  mutate(simset = as.factor(simset)) %>%
  mutate(estimator = case_when(phi1 == "ign" & phi0 == "ign" ~ "trt:ign/pla:ign",
                               phi1 == "ign" & phi0 == "mod" ~ "trt:ign/pla:mod",
                               phi1 == "ign" & phi0 == "non" ~ "trt:ign/pla:non",
                               phi1 == "mod" & phi0 == "ign" ~ "trt:mod/pla:ign",
                               phi1 == "non" & phi0 == "ign" ~ "trt:non/pla:ign")) %>%
  group_by(simset, estimator) %>%
  summarize(nsim = n(),
            cover = mean(cover)) %>%
  ungroup()

# save -----
save(cover_bayes, file="./simresults/simSummary_PPMM_Bayes_cover.RData")

