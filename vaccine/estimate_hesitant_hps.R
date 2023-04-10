# Bradley 2021 Nature paper
library(tidyverse)
library(survey)

source("https://github.com/bradytwest/IndicesOfNISB/raw/master/proxyDrawsMUBP_sumstats.R")
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/mubp_functions.R")
source("./functions/mle2stepMUBP_more_more.R")

# Population size
truth <- read.csv(file="https://media.githubusercontent.com/media/vcbradley/ddc-vaccine-US/main/data/CDC/cdc_cleaned_2021-12-03.csv") %>% 
  as_tibble()

####### Probability sample - to be treated as non-prob sample
# HPS microdata
hps <- read_csv("./data/microdata_hps_wave22to29.csv")
hps <- hps %>% mutate(hesitant = vacstatus==3)

####### Population-level estimates from ACS microdata
acs <- read_csv("./data/acs_covmat_hps.csv")
means.pop <- acs %>% filter(`_TYPE_`=="MEAN") %>% dplyr::select(-c(`_TYPE_`,`_NAME_`), -starts_with("income")) %>% as.numeric()
varcov.pop <- acs %>% filter(`_TYPE_`=="COV") %>% dplyr::select(-c(`_TYPE_`,`_NAME_`)) %>% as.matrix()

########################
######### PPMA
########################
set.seed(525)
results <- tibble()
results_parms <- tibble()
results_draws <- tibble()
#a <- Sys.time()
for (week in c(22:29))
{
  print(paste("Week",week))
  # subset survey data to one wave
  dat <- hps %>% filter(WEEK==week) %>% mutate(hesitant=(vacstatus==3))
  
  # probit fit to estimate proxy
  fit <- glm(hesitant ~ 
               male # ref=female
             + educ_somecoll + educ_bach + educ_graddeg   # ref=educ_hsless
             + race_black + race_asian + race_other # ref=race_white
             + hisp # ref=not hispanic
             + age_18_29 + age_30_39 + age_40_49 + age_50_59 + age_60_69 # ref=age_70plus
             , data=dat, binomial(link="probit"))
  # summary(fit)
  
  # Coefficients
  B <- as.vector(coef(fit))
  # Predictors Z
  z0_withInt <- model.matrix(fit)
  z0 <- z0_withInt[,-1] # remove intercept (added in the functions)
  # Create proxy for selected sample
  x0 <- z0_withInt %*% B
  # Proxy sum stats for non-selected sample (note: treated as fixed)
  x1_mean <- c(1,means.pop) %*% B
  x1_var <- t(B[-1]) %*% varcov.pop %*% B[-1]
  
  # mean(x0); x1_mean
  # var(x0); x1_var
  
  # binary outcome
  y0 <- dat$hesitant
  
  # sample sizes
  n0 <- length(x0)  # selected sample
  n <- truth$pop_total[1]
  n1 <- n - n0      # non-selected sample
  
  #########################
  # 2-Step MLE MUBP analysis
  #########################
  more0 <- mle2stepMUBP_more_more(x0, y0, x1_mean, x1_var, n0/(n0+n1), phi=0, verbose=FALSE)
  more0.5 <- mle2stepMUBP_more_more(x0, y0, x1_mean, x1_var, n0/(n0+n1), phi=0.5, verbose=FALSE)
  more1 <- mle2stepMUBP_more_more(x0, y0, x1_mean, x1_var, n0/(n0+n1), phi=1, verbose=FALSE)
  mle <- c(week=week,
           n0=n0,
           sfrac=n0/(n0+n1),
           rho_0=more0$results$rho_0,
           ymean_phi0=more0$results$ymean,
           ymean_phi0.5=more0.5$results$ymean,
           ymean_phi1=more1$results$ymean)
  mle <- as_tibble_row(mle)
  # estimated MLEs of distribution of (X,U|selected)
  mle_parms <- as_tibble_row(c(week=week, n0=n0, sfrac=n0/(n0+n1), unlist(more0$uparms)[1:8]))
  
  #########################
  # Bayesian analysis
  #########################
  draws <- proxyDrawsMSB(y0, z0, n0, means.pop, varcov.pop, n1, phi=0, drawphi=TRUE, scaleX=TRUE, nreps=2000)
  # hist(draws$muY)
  bayes <- c(quantile(draws$muY, c(.5, .025, .975)), mean(draws$muY), sd(draws$muY))
  names(bayes) <- c("ymean_p50","ymean_lb","ymean_ub","ymean_mean","ymean_sd")
  bayes <- as_tibble_row(bayes)
  draws$week <- week
  draws <- as_tibble(draws)
  
  ### Combine
  all <- bind_cols(mle, bayes)
  results <- bind_rows(results, all)
  results_parms <- bind_rows(results_parms, mle_parms)
  results_draws <- bind_rows(results_draws, draws)
}
#b <- Sys.time()
#Time difference of 18.34053 mins

#########################
# Combine with survey estimate and truth
#########################
# HPS estimates (from SAS) of percent hesitant --- using replicate weights
hps_estimates <- read_csv("./results/survey_estimates/hps_vacstatus.csv", na=c("NA","_")) %>% 
  filter(vacstatus==3) %>%
  mutate(survey_est=Percent/100, survey_se=StdErr/100, survey_lb=survey_est-2*survey_se, survey_ub=survey_est+2*survey_se) %>%
  rename(week=WEEK) %>% 
  dplyr::select(week, starts_with("survey"))

# Add to results
results <- list(hps_estimates,results) %>% reduce(left_join) %>% arrange(week)

save(results, file="./results/results_hesitant_hps.RData")
save(results_parms, file="./results/results_hesitant_parms_hps.RData")
save(results_draws, file="./results/results_hesitant_draws_hps.RData")
write_csv(results, file="./results/results_hesitant_hps.csv")
