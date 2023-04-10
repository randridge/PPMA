# Bradley 2021 Nature paper
library(tidyverse)
library(lubridate)
library(survey)

source("https://github.com/bradytwest/IndicesOfNISB/raw/master/proxyDrawsMUBP_sumstats.R")
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/mubp_functions.R")
source("./functions/mle2stepMUBP_more_more.R")

# CDC Benchmark (from Nature paper GitHub repo)
truth <- read_csv(file="https://media.githubusercontent.com/media/vcbradley/ddc-vaccine-US/main/data/CDC/cdc_cleaned_2021-12-03.csv") %>%
  mutate(week=epiweek(date), weekday=wday(date)) %>%
  filter(week>=2 & week<=18 & weekday==7) %>%
  arrange(week)

####### Probability sample - to be treated as non-prob sample
# Facebook microdata
fb <- read_csv("./data/microdata_facebook_week2to18.csv")
a <- fb %>% group_by(week) %>% summarize(n=n(), nmiss=sum(is.na(vaccinated)), pctmiss=nmiss/n*100)
b <- fb %>% dplyr::select(week, male, starts_with("educ"), starts_with("race"), starts_with("age")) %>% drop_na() %>% group_by(week) %>% summarize(ncovars=n())
c <- fb %>% dplyr::select(week, vaccinated, male, starts_with("educ"), starts_with("race"), starts_with("age")) %>% drop_na() %>% group_by(week) %>% summarize(nused=n())
abc <- left_join(a,b) %>% left_join(c) %>% mutate(pctcovars=ncovars/n*100, pctused=nused/n*100)
abc

####### Population-level estimates from ACS microdata
acs <- read_csv("./data/acs_covmat_facebook.csv")
means.pop <- acs %>% filter(`_TYPE_`=="MEAN") %>% dplyr::select(-c(`_TYPE_`,`_NAME_`)) %>% as.numeric()
varcov.pop <- acs %>% filter(`_TYPE_`=="COV") %>% dplyr::select(-c(`_TYPE_`,`_NAME_`)) %>% as.matrix()

########################
######### PPMA
########################
set.seed(1117)
results <- tibble()
results_parms <- tibble()
results_draws <- tibble()
for (w in 2:18)
{
  print(paste("Week",w))
  # subset survey data to one wave
  dat_full <- fb %>% filter(week==w)
  dim(dat_full)
  dat <- dat_full %>% dplyr::select(vaccinated, male, starts_with("educ"), starts_with("race"), starts_with("age")) %>% drop_na()
  
  # probit fit to estimate proxy
  fit <- glm(vaccinated ~ 
               male # ref=female
             + educ_hs + educ_somecoll + educ_bach + educ_graddeg   # ref=educ_lthsgrad
             + raceth_blacknh + raceth_aiannh + raceth_asianaapinh + raceth_othernh + raceth_hisp  #ref=raceth_whitenh
             + age_18_24 + age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_74 # ref=age_75plus
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
  mle <- c(week=w,
           n0=n0,
           sfrac=n0/(n0+n1),
           rho_0=more0$results$rho_0,
           ymean_phi0=more0$results$ymean,
           ymean_phi0.5=more0.5$results$ymean,
           ymean_phi1=more1$results$ymean)
  mle <- as_tibble_row(mle)
  # estimated MLEs of distribution of (X,U|selected)
  mle_parms <- as_tibble_row(c(week=w, n0=n0, sfrac=n0/(n0+n1), unlist(more0$uparms)[1:8]))

  #########################
  # Bayesian analysis
  #########################
  # a <- Sys.time()
  draws <- proxyDrawsMSB(y0, z0, n0, means.pop, varcov.pop, n1, phi=0, drawphi=TRUE, scaleX=TRUE, nreps=2000)
  # b <- Sys.time()
  # b-a
  # hist(draws$muY)
  bayes <- c(quantile(draws$muY, c(.5, .025, .975)), mean(draws$muY), sd(draws$muY))
  names(bayes) <- c("ymean_p50","ymean_lb","ymean_ub","ymean_mean","ymean_sd")
  bayes <- as_tibble_row(bayes)
  draws$week <- w
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
# Survey estimates (from API) of percent vaccinated
fb_estimates <- read_csv("./results/survey_estimates/facebook/weekly_nation_all_indicators_overall.csv") %>%
  mutate(year=floor(period_start/10000)) %>%
  filter(period_val>=2 & period_val<=18 & year==2021) %>%
  mutate(survey_est=val_pct_vaccinated/100, survey_se=se_pct_vaccinated/100, 
         survey_lb=survey_est-2*survey_se, survey_ub=survey_est+2*survey_se,
         week=period_val) %>% 
  dplyr::select(week, starts_with("survey"))
# CDC benchmark
true_ymean <- truth %>% dplyr::select(week, date, pct_pop_vaccinated) %>% rename(truth=pct_pop_vaccinated)
# Add to results
results <- list(true_ymean, fb_estimates, results) %>% reduce(left_join) %>% arrange(week)
# Add to parms
results_parms <- left_join(true_ymean, results_parms) %>% arrange(week)

save(results, file="./results/results_facebook.RData")
save(results_parms, file="./results/results_parms_facebook.RData")
save(results_draws, file="./results/results_draws_facebook.RData")
write_csv(results, file="./results/results_facebook.csv")

results %>% group_by(week) %>% summarise(week=week,
                                         date=date, 
                                         rho_0=rho_0, 
                                         cover_ppma=(ymean_lb <= truth) & (ymean_ub >= truth),
                                         cover_survey=(survey_lb <= truth) & (survey_ub >= truth))

#     week date       rho_0 cover_ppma cover_survey
#    <dbl> <date>     <dbl> <lgl>      <lgl>       
# 1      2 2021-01-16 0.237 TRUE       FALSE       
# 2      3 2021-01-23 0.272 TRUE       FALSE       
# 3      4 2021-01-30 0.331 TRUE       FALSE       
# 4      5 2021-02-06 0.381 TRUE       FALSE       
# 5      6 2021-02-13 0.431 TRUE       FALSE       
# 6      7 2021-02-20 0.472 TRUE       FALSE       
# 7      8 2021-02-27 0.505 TRUE       FALSE       
# 8      9 2021-03-06 0.529 TRUE       FALSE       
# 9     10 2021-03-13 0.532 TRUE       FALSE       
# 10    11 2021-03-20 0.525 TRUE       FALSE       
# 11    12 2021-03-27 0.502 TRUE       FALSE       
# 12    13 2021-04-03 0.479 TRUE       FALSE       
# 13    14 2021-04-10 0.442 TRUE       FALSE       
# 14    15 2021-04-17 0.423 TRUE       FALSE       
# 15    16 2021-04-24 0.401 TRUE       FALSE       
# 16    17 2021-05-01 0.386 TRUE       FALSE       
# 17    18 2021-05-08 0.377 TRUE       FALSE
