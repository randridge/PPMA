# Run methods {PPMM/Bayes}

library(tidyverse)
require(doParallel)
require(foreach)
require(doRNG)

runSet <- function(SIMSET)
{
  # Parallel processing
  n.cores <- parallel::detectCores()
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  # Looping through simulated data sets
  res <- foreach(i=1:length(SIMSET), .combine=rbind, .packages = c("tidyverse")) %dorng% {
    # Functions
    source("./functions/means.R")
    source("./functions/error_check.R")

    # Summary stats for selected (trial, by arm) and non-selected -----
    sumstats_A1 <- SIMSET[[i]]$sumstats_A1
    sumstats_A0 <- SIMSET[[i]]$sumstats_A0
    sumstats_ns <- SIMSET[[i]]$sumstats_ns
    
    # Truth -----
    truth <- SIMSET[[i]]$sumstats$truth

    # PPMM: Estimates of treatment mean -----
    # phi = ignorable (phi fixed at 0)
    draws_A1_phi_ign <- means_bayes(sumstats_A1, sumstats_ns, phi = 0, ndraws = 1000)
    # phi = moderate non-ignorable (phi ~ Uniform(0,0.5))
    draws_A1_phi_mod <- means_bayes(sumstats_A1, sumstats_ns, phi_character = "runif(n=1, min=0, max=.5)", ndraws = 1000)
    # phi = non-ignorable (phi ~ Uniform(0,1))
    draws_A1_phi_non <- means_bayes(sumstats_A1, sumstats_ns, phi_character = "runif(n=1, min=0, max=1)", ndraws = 1000)
    
    # PPMM: Estimates of control mean -----
    # phi = ignorable (phi fixed at 0)
    draws_A0_phi_ign <- means_bayes(sumstats_A0, sumstats_ns, phi = 0, ndraws = 1000)
    # phi = moderate non-ignorable (phi ~ Uniform(0,0.5))
    draws_A0_phi_mod <- means_bayes(sumstats_A0, sumstats_ns, phi_character = "runif(n=1, min=0, max=.5)", ndraws = 1000)
    # phi = non-ignorable (phi ~ Uniform(0,1))
    draws_A0_phi_non <- means_bayes(sumstats_A0, sumstats_ns, phi_character = "runif(n=1, min=0, max=1)", ndraws = 1000)
    
    # Combine estimators ------
    # Estimates of rho in each arm
    rho_per_arm <- rbind(c(quantile(draws_A1_phi_ign$rho_xy_s, c(.025, .5, .975)), quantile(draws_A0_phi_ign$rho_xy_s, c(.025, .5, .975))),
                         c(quantile(draws_A1_phi_ign$rho_xy_s, c(.025, .5, .975)), quantile(draws_A0_phi_mod$rho_xy_s, c(.025, .5, .975))),
                         c(quantile(draws_A1_phi_ign$rho_xy_s, c(.025, .5, .975)), quantile(draws_A0_phi_non$rho_xy_s, c(.025, .5, .975))),
                         c(quantile(draws_A1_phi_mod$rho_xy_s, c(.025, .5, .975)), quantile(draws_A0_phi_ign$rho_xy_s, c(.025, .5, .975))),
                         c(quantile(draws_A1_phi_non$rho_xy_s, c(.025, .5, .975)), quantile(draws_A0_phi_ign$rho_xy_s, c(.025, .5, .975)))) %>%
      as_tibble(.name_repair="unique") %>%
      rename(rho_lb_A1 = `2.5%...1`, rho_est_A1 = `50%...2`, rho_ub_A1 = `97.5%...3`, 
             rho_lb_A0 = `2.5%...4`, rho_est_A0 = `50%...5`, rho_ub_A0 = `97.5%...6`)
    # Estimates in each arm
    results_per_arm <- rbind(c(quantile(draws_A1_phi_ign$muy_ns, c(.025, .5, .975)), quantile(draws_A0_phi_ign$muy_ns, c(.025, .5, .975))),
                             c(quantile(draws_A1_phi_ign$muy_ns, c(.025, .5, .975)), quantile(draws_A0_phi_mod$muy_ns, c(.025, .5, .975))),
                             c(quantile(draws_A1_phi_ign$muy_ns, c(.025, .5, .975)), quantile(draws_A0_phi_non$muy_ns, c(.025, .5, .975))),
                             c(quantile(draws_A1_phi_mod$muy_ns, c(.025, .5, .975)), quantile(draws_A0_phi_ign$muy_ns, c(.025, .5, .975))),
                             c(quantile(draws_A1_phi_non$muy_ns, c(.025, .5, .975)), quantile(draws_A0_phi_ign$muy_ns, c(.025, .5, .975)))) %>%
      as_tibble(.name_repair="unique") %>%
      rename(lb_A1 = `2.5%...1`, est_A1 = `50%...2`, ub_A1 = `97.5%...3`, 
             lb_A0 = `2.5%...4`, est_A0 = `50%...5`, ub_A0 = `97.5%...6`)
    # Estimates of treatment mean in non-selected sample
    results <- rbind(
      # ign/ign
      quantile(draws_A1_phi_ign$muy_ns - draws_A0_phi_ign$muy_ns, c(.025, .5, .975)),
      # ign/mod
      quantile(draws_A1_phi_ign$muy_ns - draws_A0_phi_mod$muy_ns, c(.025, .5, .975)),
      # ign/non
      quantile(draws_A1_phi_ign$muy_ns - draws_A0_phi_non$muy_ns, c(.025, .5, .975)),
      # mod/ign
      quantile(draws_A1_phi_mod$muy_ns - draws_A0_phi_ign$muy_ns, c(.025, .5, .975)),
      # non/ign
      quantile(draws_A1_phi_non$muy_ns - draws_A0_phi_ign$muy_ns, c(.025, .5, .975))
    ) %>% as_tibble(.name_repair="unique") %>%
      add_column(phi1 = c("ign", "ign", "ign", "mod", "non"),
                 phi0 = c("ign", "mod", "non", "ign", "ign")) %>%
      rename(lb = `2.5%`, est = `50%`, ub = `97.5%`) %>%
      mutate(cover = (lb < truth) & (ub > truth)) %>%
      mutate(width = ub - lb) %>%
      add_column(truth = truth, rep = i) %>%
      relocate(c(rep, truth, phi1, phi0))
    # Combine
    results <- bind_cols(results, results_per_arm, rho_per_arm)

    # Output ------
    results
  }
  stopCluster(cl)
  
  # Turn results into tibble
  res <- as_tibble(res)
  return(res)
}

load("./simdata/simParms.Rdata")

setnums <- parms %>% pull(simset)

time.start <- Sys.time()
#for (i in 1:length(setnums))
for (i in 1:length(setnums))
{
  print(i)
  set.seed(531+i)
  load(paste0("./simdata/simData_set_",setnums[i],".RData"))
  temp <- runSet(datalist)
  res <- temp %>%
    add_column(simset = parmset$simset)
  write_csv(res, 
            file=paste0("./simresults/ppmm_bayes/simResults_set_",setnums[i],".csv"))
}
time.stop <- Sys.time()
time.stop - time.start
