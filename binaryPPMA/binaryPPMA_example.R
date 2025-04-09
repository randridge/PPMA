####################################
# Examples of applying the binary proxy pattern-mixture model of Andridge & Little (2019)
# Author: Rebecca R Andridge
# Last Modified: 05/22/2019
####################################

rm(list=ls())

#####################
# Load Functions    #
#####################
source("./binaryPPMA_functions.R")

#####################
# Load sample data  #
#####################
# Each example data set has the following variables (for n=400 observations):
# y = binary outcome - missing for 50% of observations
# z1, z2 = fully observed covariates
dat <- read.csv("https://raw.githubusercontent.com/randridge/PPMA/master/exampledata/binaryY_strong.csv")
# or
dat <- read.csv("https://raw.githubusercontent.com/randridge/PPMA/master/exampledata/binaryY_weak.csv")

########################################################################
# Calculate proxy (x) based on probit regression using respondent data #
########################################################################
# Proxy x
glm.fit <- glm(y ~ z1 + z2, data=dat, family=binomial(link=probit))
x <- predict(glm.fit, newdata=dat)

########################
# ML (unmodified)      #
########################
### Point estimates and large sample variance estimates
### Note the proxy x being sent into the function (not the individual Z variables)
mleFull(x, dat$y, 0)    # phi = 0
mleFull(x, dat$y, 0.5)  # phi = 0.5
mleFull(x, dat$y, 1)    # phi = 1
# Estimates of biserial correlation are returned as $rho_0

########################
# MODIFIED ML          #
########################
### Point estimates
### Note the proxy x being sent into the function (not the individual Z variables)
mle2step(x, dat$y, 0)    # phi = 0
mle2step(x, dat$y, 0.5)  # phi = 0.5
mle2step(x, dat$y, 1)    # phi = 1
### Variance estimates (using bootstrap)
require(boot)
# Define mle f'n to go into boot f'n
mle2step_boot <- function(d,i,phi)
{
  res <- mle2step(d$x[i],d$y[i],phi)
  return(res$muY)
}
### Perform bootstrap
set.seed(531); boots_phi0   <- boot(as.data.frame(cbind(x=x, y=dat$y)), mle2step_boot, R=1000, phi=0)
set.seed(531); boots_phi0.5 <- boot(as.data.frame(cbind(x=x, y=dat$y)), mle2step_boot, R=1000, phi=0.5)
set.seed(531); boots_phi1   <- boot(as.data.frame(cbind(x=x, y=dat$y)), mle2step_boot, R=1000, phi=1)
### Bootstrap variance estimates
var(boots_phi0$t)    # phi = 0
var(boots_phi0.5$t)  # phi = 0.5
var(boots_phi1$t)    # phi = 1

########################
# Posterior Draws      #
########################
# Draws from posterior for different values of phi
burnin <- 20 # number of burn-in draws for Gibbs sampler
#######
## Perform draws
# Note that the covariate matrix Z is sent into the function (not the proxy X)
set.seed(217)
phi0   <- proxyDraws(dat$y, dat[,2:3], phi=0,   nreps=2000+burnin)   # phi = 0
set.seed(217)
phi0.5 <- proxyDraws(dat$y, dat[,2:3], phi=0.5, nreps=2000+burnin)   # phi = 0.5
set.seed(217)
phi1   <- proxyDraws(dat$y, dat[,2:3], phi=1,   nreps=2000+burnin)   # phi = 1
#######
## Percentile-based credible intervals
# In the resulting matrices of parameter draws, the notation distinguishing the methods is:
# "_a"=unmodifed posterior draws; "_b"=modification 1 (redraw); "_c"=modification 2 (predprob)
# Note: BURN-IN draws should be discarded for posterior stats
# phi = 0
quantile(phi0$muY_a[-(1:burnin)], c(0.025,0.975))  # Unmodified Bayes method
quantile(phi0$muY_b[-(1:burnin)], c(0.025,0.975))  # Modification 1 (redraw) to Bayes method
quantile(phi0$muY_c[-(1:burnin)], c(0.025,0.975))  # Modification 2 (predprob) to Bayes method
# phi = 0.5
quantile(phi0.5$muY_a[-(1:burnin)], c(0.025,0.975))  # Unmodified Bayes method
quantile(phi0.5$muY_b[-(1:burnin)], c(0.025,0.975))  # Modification 1 (redraw) to Bayes method
quantile(phi0.5$muY_c[-(1:burnin)], c(0.025,0.975))  # Modification 2 (predprob) to Bayes method
# phi = 1
quantile(phi1$muY_a[-(1:burnin)], c(0.025,0.975))  # Unmodified Bayes method
quantile(phi1$muY_b[-(1:burnin)], c(0.025,0.975))  # Modification 1 (redraw) to Bayes method
quantile(phi1$muY_c[-(1:burnin)], c(0.025,0.975))  # Modification 2 (predprob) to Bayes method

########################
# Multiple Imputation
########################
## Perform imputation (50 imputations)
# Note that the covariate matrix Z is sent into the function (not the proxy X)
set.seed(217531)
mult0   <- mi(dat$y, dat[,2:3], phi=0,   drawphi=FALSE, D=50, burnin=20, thin=100)  # phi=0
set.seed(217531)
mult0.5 <- mi(dat$y, dat[,2:3], phi=0.5, drawphi=FALSE, D=50, burnin=20, thin=100)  # phi=0.5
set.seed(217531)
mult1   <- mi(dat$y, dat[,2:3], phi=1,   drawphi=FALSE, D=50, burnin=20, thin=100)  # phi=1
## Rubin's combining rules as implemented in mice package
library(mice)
# Estimated mean and variance of Y in each imputed dataset
# phi=0
est0 <- colMeans(mult0)          # vector of estimated proportions
var0 <- est0*(1-est0)/nrow(dat)  # vector of estimated variances
pool.scalar(est0,var0)           # $qbar = MI mean, $t = MI variance
# phi=0.5
est0.5 <- colMeans(mult0.5)            # vector of estimated proportions
var0.5 <- est0.5*(1-est0.5)/nrow(dat)  # vector of estimated variances
pool.scalar(est0.5,var0.5)             # $qbar = MI mean, $t = MI variance
# phi=1
est1 <- colMeans(mult1)          # vector of estimated proportions
var1 <- est1*(1-est1)/nrow(dat)  # vector of estimated variances
pool.scalar(est1,var1)           # $qbar = MI mean, $t = MI variance
