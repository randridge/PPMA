####################################
# PPMA - BINARY Y
# Using OMAS 2015 Data
# Author: Rebecca R Andridge
# Last Modified: 10/17/2019
####################################

rm(list=ls())

# Libraries
require(mvtnorm)
require(msm)
require(survey)

# Functions
source("./binaryPPMA_functions.R")

# Load data set
dat <- read.csv("./omasdata/omas2015_ppma.csv")
# NOTE: This dataset is a subset of variables from the publicly available 2015 OMAS data.
#       A few data manipulations were done to prep the data for the PPMA analysis,
#       including combining categories of auxiliary variables (e.g., insurance: "ins_cat").

# Make categorical variables into factors
facts <- c("strata","mc_region","Region","gender","age_a_imp","educ_imp","race4_a_imp","ins_cat")
dat[facts] <- lapply(dat[facts], factor)
rm(facts)

# Create survey design object
omasdesign <- svydesign(id=~1, strata=~strata, weights=~WT_A, data=dat)
omasdesign

##############################################
# Choose the outcome here                    #
# (only run one of these two chunks of code) #
##############################################
####################################
# OPTION 1:                        #
# Outcome: Low Income (<138 FPL)   #
# Final Model to Create the Proxy X (from backwards selection)
fit <- glm(fpl_yr_lt_138 ~ strata + WT_A + Region + gender + age_a_imp + educ_imp + race4_a_imp + ins_cat + fairpoor + S11_IMP + S12_IMP + 
                           WT_A*age_a_imp + Region*age_a_imp + gender*age_a_imp + WT_A*educ_imp + gender*educ_imp + age_a_imp*educ_imp + 
                           WT_A*race4_a_imp + gender*race4_a_imp + age_a_imp*race4_a_imp + gender*ins_cat + age_a_imp*ins_cat + 
                           educ_imp*ins_cat + race4_a_imp*ins_cat + fairpoor*Region + fairpoor*age_a_imp + fairpoor*ins_cat + WT_A*S11_IMP + 
                           S11_IMP*educ_imp + S11_IMP*ins_cat + WT_A*S12_IMP + S12_IMP*Region + S12_IMP*gender + S12_IMP*age_a_imp + 
                           S12_IMP*educ_imp + S12_IMP*ins_cat + S11_IMP*S12_IMP,
                           data=dat, family=binomial(link="probit"))

# Pull out design matrix and outcome
z <- model.matrix( ~ strata + WT_A + Region + gender + age_a_imp + educ_imp + race4_a_imp + ins_cat + fairpoor + S11_IMP + S12_IMP + 
                     WT_A*age_a_imp + Region*age_a_imp + gender*age_a_imp + WT_A*educ_imp + gender*educ_imp + age_a_imp*educ_imp + 
                     WT_A*race4_a_imp + gender*race4_a_imp + age_a_imp*race4_a_imp + gender*ins_cat + age_a_imp*ins_cat + 
                     educ_imp*ins_cat + race4_a_imp*ins_cat + fairpoor*Region + fairpoor*age_a_imp + fairpoor*ins_cat + WT_A*S11_IMP + 
                     S11_IMP*educ_imp + S11_IMP*ins_cat + WT_A*S12_IMP + S12_IMP*Region + S12_IMP*gender + S12_IMP*age_a_imp + 
                     S12_IMP*educ_imp + S12_IMP*ins_cat + S11_IMP*S12_IMP,
                     data=dat)
# Remove intercept; mi() function needs the column for intercept removed from the dataset
z <- z[,-1]
# Pull off binary outcome variable (with missingness)
y <- dat$fpl_yr_lt_138
# Indicator for missing Y
m <- is.na(y)
####################################
# OPTION 2:                        #
# Outcome: High Income (>300 FPL)  #
# Final Model
fit <- glm(fpl_yr_gt_300 ~ strata + WT_A + Region + gender + age_a_imp + educ_imp + race4_a_imp + ins_cat + fairpoor + S11_IMP + S12_IMP + 
             WT_A*gender + WT_A*age_a_imp + Region*age_a_imp + gender*age_a_imp + WT_A*educ_imp + gender*educ_imp + 
             age_a_imp*educ_imp + gender*race4_a_imp + age_a_imp*race4_a_imp + gender*ins_cat + age_a_imp*ins_cat + 
             educ_imp*ins_cat + WT_A*S11_IMP + S11_IMP*gender + S11_IMP*age_a_imp + S11_IMP*educ_imp + S11_IMP*ins_cat + 
             S12_IMP*gender + S12_IMP*age_a_imp + S12_IMP*race4_a_imp + S12_IMP*ins_cat + fairpoor*S12_IMP + S11_IMP*S12_IMP,
           data=dat, family=binomial(link="probit"))

# Pull out design matrix and outcome
z <- model.matrix( ~ strata + WT_A + Region + gender + age_a_imp + educ_imp + race4_a_imp + ins_cat + fairpoor + S11_IMP + S12_IMP + 
                     WT_A*gender + WT_A*age_a_imp + Region*age_a_imp + gender*age_a_imp + WT_A*educ_imp + gender*educ_imp + 
                     age_a_imp*educ_imp + gender*race4_a_imp + age_a_imp*race4_a_imp + gender*ins_cat + age_a_imp*ins_cat + 
                     educ_imp*ins_cat + WT_A*S11_IMP + S11_IMP*gender + S11_IMP*age_a_imp + S11_IMP*educ_imp + S11_IMP*ins_cat + 
                     S12_IMP*gender + S12_IMP*age_a_imp + S12_IMP*race4_a_imp + S12_IMP*ins_cat + fairpoor*S12_IMP + S11_IMP*S12_IMP, 
                   data=dat)
z <- z[,-1]  # Remove intercept
y <- dat$fpl_yr_gt_300
m <- is.na(y)

##################################################
# Proceed with all code below for chosen outcome #
##################################################

############################
# Calculate proxy X*, RHO, D
############################
# These values (% missing, rho_0, and d*) are reported in the paper in the text and also on Figure 2 (top of plots)
# Create proxy (and proxy transformed to probability scale) for MLEs
x <- predict.glm(fit, type="link", newdata=dat)
# Percent missing
100 * mean(m)
# rho_0 = Proxy strength
#         Note that it doesn't matter what you choose for phi (sensitivity parameter),
#         since the estimate of the correlation (rho_0) is not dependent on phi
mleFull(x,y,0)$rho_0   # ML
mle2step(x,y,0)$rho_0  # Modified ML ("two-step" method)
# d* = Standardized difference between mean of proxy for whole sample compard to respondents
(mean(x)-mean(x[m==0]))/sqrt(var(x[m==0]))

########################
# Complete Case Estimates
########################
ccest <- svymean(~y, omasdesign, na.rm=T)
ccWt <- data.frame(cbind(mean=as.numeric(coef(ccest)),
                         se=as.numeric(sqrt(vcov(ccest)))))
N_obs <- sum(1-m)
H <- length(unique(dat$strata))
ccWt$lb <- ccWt$mean - qt(0.025, N_obs-H, lower.tail=F)*ccWt$se
ccWt$ub <- ccWt$mean + qt(0.025, N_obs-H, lower.tail=F)*ccWt$se
rm(N_obs,H,ccest)
# Print complete case estimates (plotted in the paper in Figure 2)
ccWt

########################
# Multiple Imputation
########################
# Perform multiple imputation with the binary PPM model
# Result are matrices of imputed data (rows=subjects, columns=multiple imputations)
# Note that these may take a while to run
system.time({set.seed(5254772); mult0   <- mi(y,z,phi=0,   drawphi=FALSE,D=20,burnin=20,thin=100)})
system.time({set.seed(5254772); mult1   <- mi(y,z,phi=0.5, drawphi=FALSE,D=20,burnin=20,thin=100)})
system.time({set.seed(5254772); multInf <- mi(y,z,phi=1,   drawphi=FALSE,D=20,burnin=20,thin=100)})

# Function to calculate stats from MI data using sample design
miStatsWts <- function(MIDATA, DESIGN)
{
  D <- ncol(MIDATA)
  thetaHats <- apply(MIDATA, 2, function(x) coef(svymean(~x, DESIGN)))
  Ws <- apply(MIDATA, 2, function(x) vcov(svymean(~x, DESIGN)))
  thetaHat <- mean(thetaHats)
  W <- mean(Ws)
  B <- (1/(D-1))*sum((thetaHats-thetaHat)^2)
  T <- W + ((D+1)/D)*B
  df <- (D-1)*(1 + (D/(D+1))*(W/B))^2
  FMI <- ((D+1)/D)*B/T
  lb <- thetaHat - qt(0.025, df, lower.tail=F)*sqrt(T)
  ub <- thetaHat + qt(0.025, df, lower.tail=F)*sqrt(T)
  stats <- as.data.frame(cbind(mean=thetaHat, var=T, lb, ub, df, FMI))
  return(stats)
}

temp <- as.data.frame(rbind(miStatsWts(mult0,omasdesign),
                            miStatsWts(mult1,omasdesign),
                            miStatsWts(multInf,omasdesign)))
temp$se <- sqrt(temp$var)
multWt <- subset(temp, select=c(mean,se,lb,ub,df,FMI))
# Print MI PPMA estimates (plotted in the paper in Figure 2)
multWt

