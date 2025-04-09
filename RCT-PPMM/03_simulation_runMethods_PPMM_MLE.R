# Run methods {PPMM/MLE}

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
  res <- foreach(i=1:length(SIMSET), .packages = c("tidyverse")) %dorng% {
    # Functions
    source("./functions/ppmm_rct_mle.R")
    
    # Summary stats for selected (trial, by arm) and non-selected -----
    sumstats_A1 <- SIMSET[[i]]$sumstats_A1
    sumstats_A0 <- SIMSET[[i]]$sumstats_A0
    sumstats_ns <- SIMSET[[i]]$sumstats_ns
    
    # PPMM MLEs ------
    set1 <- ppmm_rct_mle(sumstats_A1$mean_YZ, sumstats_A1$var_YZ,
                        sumstats_A0$mean_YZ, sumstats_A0$var_YZ,
                        sumstats_ns$mean_Z,  sumstats_ns$var_Z,
                        PHI_A1 = c(0,0.5,1), PHI_A0 = 0)
    set2 <- ppmm_rct_mle(sumstats_A1$mean_YZ, sumstats_A1$var_YZ,
                         sumstats_A0$mean_YZ, sumstats_A0$var_YZ,
                         sumstats_ns$mean_Z,  sumstats_ns$var_Z,
                         PHI_A1 = 0, PHI_A0 = c(0.5,1))
    ppmm_mles <- bind_rows(set1$results, set2$results) %>%
      add_column(rep = i) %>%
      relocate(rep)

    # output ------
    list(proxyinfo=set1$proxyinfo, ppmm_mles=ppmm_mles)
  }
  stopCluster(cl)
  
  # Turn results into tibbles
  proxyinfo <- purrr::map(res, 1) %>% list_rbind()
  ppmm_mles <- purrr::map(res, 2) %>% list_rbind()
  
  return(list(proxyinfo, ppmm_mles))
}

filelist <- fs::dir_ls("./simdata/", regexp = "simData")
time.start <- Sys.time()
for (i in 1:length(filelist))
{
  print(i)
  load(filelist[i])
  temp <- runSet(datalist)
  proxyinfo <- temp[[1]] %>%
    add_column(simset = parmset$simset) %>% 
    relocate(simset)
  ppmm_mles <- temp[[2]] %>%
    add_column(simset = parmset$simset) %>% 
    relocate(simset)
  write_csv(proxyinfo,
            file=sub("./simdata/simData", "./simresults/ppmm_mle/proxyinfo", sub(".RData", ".csv", filelist[i])))
  write_csv(ppmm_mles,
            file=sub("./simdata/simData", "./simresults/ppmm_mle/simResults", sub(".RData", ".csv", filelist[i])))
}
time.stop <- Sys.time()
time.stop - time.start

