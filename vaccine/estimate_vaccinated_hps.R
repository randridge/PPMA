# Bradley 2021 Nature paper
library(tidyverse)
library(survey)

source("https://github.com/bradytwest/IndicesOfNISB/raw/master/proxyDrawsMUBP_sumstats.R")
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/mubp_functions.R")
source("./functions/mle2stepMUBP_more_more.R")

# CDC Benchmark (from Nature paper GitHub repo)
truth <- read.csv(file="https://media.githubusercontent.com/media/vcbradley/ddc-vaccine-US/main/data/CDC/cdc_cleaned_2021-12-03.csv") %>% 
  as_tibble() %>% 
  filter(date %in% c("2021-01-18","2021-02-01","2021-02-15","2021-03-01","2021-03-15","2021-03-29","2021-04-26","2021-05-10"))
truth$WEEK <- c(29:22)

####### Probability sample - to be treated as non-prob sample
# HPS microdata
hps <- read_csv("./data/microdata_hps_wave22to29.csv")

####### Population-level estimates from ACS microdata
acs <- read_csv("./data/acs_covmat_hps.csv")
means.pop <- acs %>% filter(`_TYPE_`=="MEAN") %>% dplyr::select(-c(`_TYPE_`,`_NAME_`)) %>% as.numeric()
varcov.pop <- acs %>% filter(`_TYPE_`=="COV") %>% dplyr::select(-c(`_TYPE_`,`_NAME_`)) %>% as.matrix()

########################
######### PPMA
########################
set.seed(525)
results <- tibble()
results_parms <- tibble()
results_draws <- tibble()
for (week in c(22:29))
{
  print(paste("Week",week))
  # subset survey data to one wave
  dat <- hps %>% filter(WEEK==week)

  # probit fit to estimate proxy
  fit <- glm(vaccinated ~ 
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
  y0 <- dat$vaccinated
  
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

#########################
# Combine with survey estimate and truth
#########################
# HPS estimates (from SAS) of percent vaccinated
hps_estimates <- read_csv("./results/survey_estimates/hps_pct_vacc.csv", na=c("NA","_")) %>% 
  filter(vaccinated==1) %>%
  mutate(survey_est=Percent/100, survey_se=StdErr/100, survey_lb=survey_est-2*survey_se, survey_ub=survey_est+2*survey_se) %>%
  rename(week=WEEK) %>% 
  dplyr::select(week, starts_with("survey"))
# CDC benchmark
true_ymean <- truth %>% dplyr::select(WEEK, date, pct_pop_vaccinated) %>% rename(week=WEEK, truth=pct_pop_vaccinated)
# Add to results
results <- list(true_ymean, hps_estimates,results) %>% reduce(left_join) %>% arrange(week)
# Add to parms
results_parms <- left_join(true_ymean, results_parms) %>% arrange(week)

save(results, file="./results/results_hps.RData")
save(results_parms, file="./results/results_parms_hps.RData")
save(results_draws, file="./results/results_draws_hps.RData")
write_csv(results, file="./results/results_hps.csv")

results %>% group_by(week) %>% summarise(week=week,
                                         date=date, 
                                         rho_0=rho_0, 
                                         cover_ppma=(ymean_lb <= truth) & (ymean_ub >= truth),
                                         cover_survey=(survey_lb <= truth) & (survey_ub >= truth))

results %>% dplyr::select(week, truth, survey_est, ymean_p50)

#    week date       rho_0 cover_ppma cover_survey
#   <dbl> <chr>      <dbl> <lgl>      <lgl>       
# 1    22 2021-01-18 0.244 TRUE       TRUE        
# 2    23 2021-02-01 0.267 TRUE       TRUE        
# 3    24 2021-02-15 0.372 TRUE       FALSE       
# 4    25 2021-03-01 0.458 TRUE       FALSE       
# 5    26 2021-03-15 0.517 TRUE       FALSE       
# 6    27 2021-03-29 0.524 TRUE       FALSE       
# 7    28 2021-04-26 0.452 TRUE       FALSE       
# 8    29 2021-05-10 0.439 TRUE       FALSE             
