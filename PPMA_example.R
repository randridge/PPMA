####################################
# Examples of applying the proxy pattern-mixture model of Andridge & Little (2011)
# Author: Rebecca R Andridge
# Last Modified: 11/19/2020
####################################

#####################
# Load Functions    #
#####################
source("./PPMA_functions.R")

#####################
# Load sample data  #
#####################
dat <- read.csv("https://raw.githubusercontent.com/randridge/PPMA/master/exampledata/ofhs.csv")
# Make categorical variables into factors
dat$racecat <- factor(dat$racecat)
dat$educcat <- factor(dat$educcat)
dat$region <- factor(dat$region)
# NOTE: These data are actually from a complex sampling design,
#       thus using MI with weights used for post-MI estimation should be used
#       for population-level estimates.

########################################################################
# Calculate proxy (x) based on linear regression using respondent data #
########################################################################
# Proxy x
# outcome = logincome
fit <- lm(logincome ~ region + numadults + numchildren + age + female + racecat + educcat + medicaid + insured + ownhome + WT_A, data=dat)
x <- predict(fit, newdata=dat)

########################
# MLE                  #
########################
### Point estimates and large sample variance estimates
### Note the proxy x being sent into the function (not the individual predictor variables)
mle(x, dat$logincome, 0)    # lambda = 0
mle(x, dat$logincome, 1)    # lambda = 1
mle(x, dat$logincome, Inf)  # lambda = Inf

########################
# Posterior Draws      #
########################
# Get the model matrix for ALL subjects
z <- model.matrix(~ region + numadults + numchildren + age + female + racecat + educcat + medicaid + insured + ownhome + WT_A, data=dat)
# Remove intercept (added automatically in proxyDraws())
z <- z[,-1]
#######
## Compute draws
# Note that the covariate matrix Z is sent into the function (not the proxy X)
set.seed(217)
lam0   <- proxyDraws(dat$logincome, z, lambda=0,   nreps=2000)   # lambda=0
set.seed(217)
lam1   <- proxyDraws(dat$logincome, z, lambda=1,   nreps=2000)   # lambda=1
set.seed(217)
lamInf <- proxyDraws(dat$logincome, z, lambda=Inf, nreps=2000)   # lambda=Inf
#######
## Percentile-based credible intervals
# lambda = 0
quantile(lam0, c(0.025,0.975))
# lambda = 1
quantile(lam1, c(0.025,0.975))
# lambda = Inf
quantile(lamInf, c(0.025,0.975))

########################
# Multiple Imputation  
########################
## Perform imputation (50 imputations)
# Note that the covariate matrix Z is sent into the function (not the proxy X)
set.seed(727984)
mult0   <- mi(dat$logincome, z, 0, 50)
set.seed(727984)
mult1   <- mi(dat$logincome, z, 1, 50)
set.seed(727984)
multInf <- mi(dat$logincome, z, Inf, 50)
## Rubin's combining rules as implemented in mice package
library(mice)
# Estimated mean and variance of Y in each imputed dataset
# lambda=0
est0 <- colMeans(mult0)                                # vector of estimated means
var0 <- apply(mult0, 2, function(x) var(x)/length(x))  # vector of estimated variances
pool.scalar(est0,var0)                                 # $qbar = MI mean, $t = MI variance
# lambda=1
est1 <- colMeans(mult1)                                # vector of estimated means
var1 <- apply(mult1, 2, function(x) var(x)/length(x))  # vector of estimated variances
pool.scalar(est1,var1)                                 # $qbar = MI mean, $t = MI variance
# lambda=Inf
estInf <- colMeans(multInf)                                # vector of estimated means
varInf <- apply(multInf, 2, function(x) var(x)/length(x))  # vector of estimated variances
pool.scalar(estInf,varInf)                                 # $qbar = MI mean, $t = MI variance

