# Bradley 2021 Nature paper
library(tidyverse)
library(survey)

source("https://github.com/bradytwest/IndicesOfNISB/raw/master/proxyDrawsMUBP_sumstats.R")
source("https://github.com/bradytwest/IndicesOfNISB/raw/master/mubp_functions.R")

# Population size
truth <- read.csv(file="https://media.githubusercontent.com/media/vcbradley/ddc-vaccine-US/main/data/CDC/cdc_cleaned_2021-12-03.csv") %>% 
  as_tibble()

####### Probability sample - to be treated as non-prob sample
# DFP microdata
dfp <- read_csv("./data/microdata_dfp_wave20to25.csv") %>% rename(WAVE=wave)
dfp <- dfp %>% mutate(hesitant = vacstatus==3)

####### Population-level estimates from ACS microdata
acs <- read_csv("./data/acs_covmat_hps.csv")
means.pop <- acs %>% filter(`_TYPE_`=="MEAN") %>% dplyr::select(-c(`_TYPE_`,`_NAME_`)) %>% as.numeric()
varcov.pop <- acs %>% filter(`_TYPE_`=="COV") %>% dplyr::select(-c(`_TYPE_`,`_NAME_`)) %>% as.matrix()

########################
######### PPMA
########################
set.seed(531217)
results <- tibble()
for (wave in 20:25)
{
  print(paste("Wave",wave))
  # subset survey data to one wave
  dat <- dfp %>% filter(WAVE==wave)
  
  # probit fit to estimate proxy using survey data (unweighted)
  fit <- glm(hesitant ~ 
               male # ref=female
             + educ_somecoll + educ_bach + educ_graddeg   # ref=educ_hsless
             + race_black + race_asian + race_other # ref=race_white
             + hisp # ref=not hispanic
             + age_18_29 + age_30_39 + age_40_49 + age_50_59 + age_60_69 # ref=age_70plus
             + incomemiss_0 + incomemiss_1 + incomemiss_2 + incomemiss_3 + incomemiss_4 + incomemiss_5 + incomemiss_6 + incomemiss_7 # ref=incomemiss_8
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
  y0 <- fit$model$hesitant  # only the cases used in the model fit, i.e., complete cases

  # sample sizes
  n0 <- length(x0)  # selected sample
  n <- truth$pop_total[1]
  n1 <- n - n0      # non-selected sample
  
  #########################
  # 2-Step MLE MUBP analysis
  #########################
  more0 <- mle2stepMUBP_more(x0, y0, x1_mean, x1_var, n0/(n0+n1), phi=0, verbose=FALSE)
  more0.5 <- mle2stepMUBP_more(x0, y0, x1_mean, x1_var, n0/(n0+n1), phi=0.5, verbose=FALSE)
  more1 <- mle2stepMUBP_more(x0, y0, x1_mean, x1_var, n0/(n0+n1), phi=1, verbose=FALSE)
  mle <- c(wave=wave,
           rho_0=more0$rho_0,
           ymean_phi0=more0$ymean,
           ymean_phi0.5=more0.5$ymean,
           ymean_phi1=more1$ymean)
  mle <- as_tibble_row(mle)
  
  #########################
  # Bayesian analysis
  #########################
  draws <- proxyDrawsMSB(y0, z0, n0, means.pop, varcov.pop, n1, phi=0, drawphi=TRUE, scaleX=TRUE, nreps=2000)
  # hist(draws$muY)
  bayes <- c(quantile(draws$muY, c(.5, .025, .975)), mean(draws$muY), sd(draws$muY))
  names(bayes) <- c("ymean_p50","ymean_lb","ymean_ub","ymean_mean","ymean_sd")
  bayes <- as_tibble_row(bayes)
  
  ### Combine
  all <- bind_cols(mle, bayes)
  results <- bind_rows(results, all)
}


#########################
# Combine with survey estimate and truth
#########################
# Axios-Ipsos estimates (from SAS) of percent hesitant
dfp_estimates <- read_csv("./results/dfp_vacstatus.csv", na=c("NA","_")) %>% 
  filter(vacstatus==3) %>%
  mutate(survey_est=Percent/100, survey_se=StdErr/100, survey_lb=survey_est-2*survey_se, survey_ub=survey_est+2*survey_se) %>% 
  dplyr::select(wave, starts_with("survey"))
# Add to results
results <- left_join(results,dfp_estimates)
#      wave rho_0 ymean_phi0 ymean_phi0.5 ymean_phi1 ymean_p50 ymean_lb ymean_ub ymean_mean ymean_sd survey_est survey_se survey_lb survey_ub
#     <dbl> <dbl>      <dbl>        <dbl>      <dbl>     <dbl>    <dbl>    <dbl>      <dbl>    <dbl>      <dbl>     <dbl>     <dbl>     <dbl>
#   1    20 0.509      0.299        0.284     0.258      0.286   0.221     0.329      0.283   0.0279      0.391    0.0170     0.357     0.425
#   2    21 0.441      0.294        0.287     0.0749     0.286   0.147     0.352      0.279   0.0510      0.350    0.0158     0.318     0.382
#   3    22 0.526      0.306        0.306     0.313      0.304   0.266     0.348      0.304   0.0206      0.366    0.0166     0.333     0.399
#   4    23 0.470      0.275        0.246     0.0257     0.251   0.0848    0.311      0.239   0.0541      0.343    0.0193     0.304     0.381
#   5    24 0.472      0.263        0.289     0.349      0.286   0.236     0.387      0.294   0.0389      0.299    0.0169     0.265     0.332
#   6    25 0.379      0.235        0.231     0.230      0.235   0.165     0.309      0.236   0.0355      0.271    0.0140     0.243     0.299
