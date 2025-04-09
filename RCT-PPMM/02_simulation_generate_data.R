library(tidyverse)
require(doParallel)
require(foreach)
require(doRNG)

# load simulation parms
load("./simdata/simParms.RData")

# Loop through and create/save datasets
for (i in parms$simset)
{
  print(i)
  # Subset to parameters for a specific selection mechanism
  parmset <- parms %>% filter(simset == i)
  # Create data
  set.seed(231517)
#start <- Sys.time()
  # Parallel processing
  n.cores <- parallel::detectCores()
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  # Looping through simulated data sets
  results <- foreach(i=1:1000, .packages = c("tidyverse")) %dorng% {
    
    # parms
    n <- 5000
    SFRAC <- parmset$SFRAC
    B0 <- parmset$B0
    B_Z <- parmset$B_Z
    B_W <- parmset$B_W
    B_U <- parmset$B_U
    C_Z <- parmset$C_Z
    C_W <- parmset$C_W
    C_U <- parmset$C_U
    CORR_U_W <- parmset$CORR_U_W
    S2 <- parmset$S2
    
    # sample size
    n_RCT <- SFRAC * n
    
    # covariates ------
    z <- rnorm(n)
    u <- rnorm(n)
    w <- CORR_U_W*u + rnorm(n, mean = 0, sd = sqrt(1 - CORR_U_W^2))
    
    # selection into trial - fixed sample size -----
    xb <- B0 + B_Z*z + B_W*w + B_U*u
    prob <- exp(xb)/(1+exp(xb))
    pik <- sampling::inclusionprobabilities(prob, n_RCT)
    s <- sampling::UPbrewer(pik)
    # trt assignment - equally split -----
    a <- vector(length=n)
    # first half treatment, second half control
    a[s==1] <- c(rep(1, n_RCT/2), rep(0, n_RCT/2))
    # will assign a=0 for cases with s=0 (not in trial)
    n_RCT_trt <- sum(a)
    n_RCT_pla <- n_RCT - n_RCT_trt
    
    # potential outcomes ------
    y_0 <- 0 + z + w + u + rnorm(n, 0, sqrt(S2))
    y_1 <- 1 + (1+C_Z)*z + (1+C_W)*w + (1+C_U)*u + rnorm(n, 0, sqrt(S2))
    
    # observed outcome ------
    y <- a*y_1 + (1-a)*y_0
    y[s==0] <- NA
    
    # tibble ------
    dat <- tibble(S=s, A=a, Y_0=y_0, Y_1=y_1, Z=z, W=w, U=u) %>%
      mutate(Y = case_when(S == 1 & A == 0 ~ Y_0,
                           S == 1 & A == 1 ~ Y_1,
                           .default = NA))
    dat_A1 <- dat %>%
      filter(S == 1 & A == 1) %>%
      select(Y, Z, W)
    dat_A0 <- dat %>%
      filter(S == 1 & A == 0) %>%
      select(Y, Z, W)
    dat_ns <- dat %>%
      filter(S == 0) %>%
      select(Z, W)

    # summary statistics for selected by arm (in trial), non-selected (not in trial) ------
    sumstats_A1 <- list(mean_YZ = colMeans(dat_A1),
                        var_YZ = var(dat_A1),
                        n_YZ = nrow(dat_A1))
    sumstats_A0 <- list(mean_YZ = colMeans(dat_A0),
                        var_YZ = var(dat_A0),
                        n_YZ = nrow(dat_A0))
    sumstats_ns <- list(mean_Z = colMeans(dat_ns),
                        var_Z = var(dat_ns),
                        n_Z = nrow(dat_ns))
 
    # true effect in NS sample -----
    truth_A1 <- mean(dat$Y_1[dat$S==0])
    truth_A0 <- mean(dat$Y_0[dat$S==0])
    truth <- truth_A1 - truth_A0
    
    # Trial estimator -----
    trial_A1 <- mean(dat$Y_1[dat$S==1 & dat$A==1])
    trial_A0 <- mean(dat$Y_0[dat$S==1 & dat$A==0])
    trial <- trial_A1 - trial_A0
    
    # combine
    sumstats <- tibble(rep=i, truth_A1, truth_A0, truth, trial_A1, trial_A0, trial)
    
    return(list(sumstats_A1 = sumstats_A1,
                sumstats_A0 = sumstats_A0,
                sumstats_ns = sumstats_ns,
                sumstats = sumstats))
  }
  stopCluster(cl)
#stop <- Sys.time()
  
  # Extract results and save
  truth_trial <- purrr::map(results, 4) %>% list_rbind() %>% mutate(simset = i)
  datalist <- results
  
  # Save
  write_csv(truth_trial, paste0("./simresults/truth_trial/simResults_set_",i,".csv"))
  save(parmset, datalist, file=paste("./simdata/simData_set_",i,".RData", sep=""))
}
