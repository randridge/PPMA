####################################################################################
# R functions implementing binary proxy pattern-mixture analysis
#
# Author: 
# Rebecca Andridge (andridge.1@osu.edu)
#
# Version: May 22 2019
####################################################################################

#####################################################
# ML (unmodified) estimate of mean of Y
# for a given phi in [0,1]
#
# Authors: Rebecca Andridge
# Last Modified: 05/15/19
#
# Inputs:
#    x = fully observed proxy vector
#    y = partially observed binary outcome vector
#    phi = sensitivity parameter in [0,1]
#
# Returns:
#    rho_0 = biserial correlation estimate
#    muY_0 = estimated mean of Y for respondents
#    muY_1 = estimated mean of Y for nonrespondents
#    muY = estimated mean of Y (overall)
#    muYvar = estimated variance of the mean of Y (overall)
#####################################################
require(msm)
mleFull <- function(x, y, phi)
{
  # Indicator for missingness
  m <- ifelse(is.na(y), 1, 0)
  # Number of respondents, nonrespondents
  n <- length(y)
  r <- sum(1-m)
  
  # Separate out respondent and non-respondent vectors
  x_0 <- x[m==0]
  y_0 <- y[m==0]
  x_1 <- x[m==1]
  
  # Identified parameters
  # E[X|M=0]
  muX_0 <- mean(x_0)
  # Var[X|M=0]
  sigmaXX_0 <- sum((x_0 - muX_0)^2)/r
  # E[X|M=1]
  muX_1 <- mean(x_1)
  # Var[X|M=1]
  sigmaXX_1 <- sum((x_1 - muX_1)^2)/(n-r)
  # Pr(M=1)
  pi <- r/n
  
  # Polyserial correlation and threshold
  ## FULL MAXIMUM LIKELIHOOD
  # Uses transform of (w,p) --> (a,b)
  # Maximize likelihood wrt (a,b)
  f <- function(pars)
  {
    a <- pars[1]
    b <- pars[2]
    logPhi <- pnorm(a + b*x_0, log.p=TRUE)
    log1minusPhi <- pnorm(a + b*x_0, log.p=TRUE, lower.tail=FALSE)
    -sum(y_0*logPhi + (1-y_0)*log1minusPhi)
  }
  result <- optim(c(0, 1), f, method="BFGS", hessian=TRUE)
  a <- result$par[1]
  b <- result$par[2]
  vcov.ab <- solve(result$hessian)
  ## other variances
  v.muX_0 <- sigmaXX_0/r
  v.sigmaXX_0 <- 2*sigmaXX_0^2/r
  v.muX_1 <- sigmaXX_1/(n-r)
  v.sigmaXX_1 <- 2*sigmaXX_1^2/(n-r)
  vcov.X <- diag(c(v.muX_0, v.sigmaXX_0, v.muX_1, v.sigmaXX_1))
  # Full variance-coariance matrix for (muX_0,sigmaXX_0,muX_1,sigmaXX_1,a,b)
  upcorner <- matrix(data=0, nrow=4, ncol=2)
  vcov <- rbind(cbind(vcov.X, upcorner),
                cbind(t(upcorner), vcov.ab))
  # Full mean vector
  meanests <- c(muX_0, sigmaXX_0, muX_1, sigmaXX_1, a, b)
  
  # Transform back to (muU_0,rho_o)
  rho_0 <- sqrt(b^2*sigmaXX_0/(1+b^2*sigmaXX_0))
  muU_0 <- sqrt(1-rho_0^2)*(a + b*muX_0)
  
  ## MLEs for distribution of U
  g <- (phi+(1-phi)*rho_0)/(phi*rho_0+(1-phi))
  sigmaXU_0 <- rho_0*sqrt(sigmaXX_0)
  muU_1 <- muU_0 + g*(muX_1 - muX_0)/sqrt(sigmaXX_0)
  sigmaUU_1 <- 1 + (1/sigmaXX_0)*g^2*(sigmaXX_1-sigmaXX_0)
  sigmaXU_1 <- sigmaXU_0 + sqrt(1/sigmaXX_0)*g*(sigmaXX_1-sigmaXX_0)
  # If sigmaUU_1 < 0 replace with boundary value .Machine$double.eps
  # This will cause denominator of muY_1 estimate to be very small,
  # hence muY_1=pnorm(very large magnitude)=1 or -1 depending on sign of muU_1
  sigmaUU_1 <- ifelse(sigmaUU_1<0, .Machine$double.eps, sigmaUU_1)
  
  
  ## MLEs for E[Y|M=0], E[Y|M=1], E[Y]
  ### E[Y|M=0]
  muY_0 <- pnorm(muU_0)
  ### E[Y|M=1]
  muY_1 <- pnorm(muU_1/sqrt(sigmaUU_1))
  ### MU_Y
  muY <- pi*muY_0 + (1-pi)*muY_1
  
  ## Variance of MLEs
  # This code requires the sensitivity parameter to be parameterized as lambda 
  # (as in Andridge & Little 2011) instead of as phi (as in Andridge & Little 2019).
  lambda <- phi/(1-phi)
  
  # Notation for 1st deltamethod function
  #   x1=muX_0; x2=sigmaXX_0; x3=muX_1; x4=sigmaXX_1; x5=a; x6=b
  # Notation for remaining deltamethod functions
  #   x1=w=-muU_0; x2=muU_1; x3=sigmaUU_1  
  
  # Transform from (muX_0,sigmaXX_0,muX_1,sigmaXX_1,a,b) --> (w,muU_1,sigmaUU_1)  where w = -muU_0
  wMuSig <- c(-muU_0, muU_1, sigmaUU_1)
  if (lambda==Inf){
    formMuU <- sprintf("~(x5+x6*x1)/sqrt(1+x6^2*x2) + sqrt((1+x6^2*x2)/(x6^2*x2)) * (x3-x1)/sqrt(x2)")
    formSigU <- sprintf("~1+(1+x6^2*x2)/(x6^2*x2)*(x4/x2-1)")
  } else {
    formMuU <- sprintf("~(x5+x6*x1)/sqrt(1+x6^2*x2) + (%f+sqrt((x6^2*x2)/(1+x6^2*x2)))/(%f*sqrt((x6^2*x2)/(1+x6^2*x2)) + 1) * (x3-x1)/sqrt(x2)",lambda,lambda)
    formSigU <- sprintf("~1+((%f+sqrt((x6^2*x2)/(1+x6^2*x2)))/(%f*sqrt((x6^2*x2)/(1+x6^2*x2)) + 1))^2*(x4/x2-1)",lambda,lambda)
  }
  v.wMuSig <- deltamethod(list(~-(x5+x6*x1)/sqrt(1+x6^2*x2),
                               as.formula(formMuU),
                               as.formula(formSigU)
  ), meanests, vcov, ses=F)
  muY_0var <- deltamethod(~pnorm(-x1), wMuSig, v.wMuSig, ses=F)
  muY_1var <- deltamethod(~pnorm(x2/sqrt(x3)), wMuSig, v.wMuSig, ses=F)
  form <- sprintf("~%f*pnorm(-x1) + (1-%f)*pnorm(x2/sqrt(x3))", pi,pi)
  muYvar <- deltamethod(as.formula(form), wMuSig, v.wMuSig, ses=F)
  
  return(list(rho_0=rho_0,
              muY_0=muY_0,
              muY_1=muY_1,
              muY=muY,
              muYvar=as.numeric(muYvar)))
}

#####################################################
# Modified ML (two-step) estimate of mean of Y
# for a given phi in [0,1]
#
# Authors: Rebecca Andridge
# Last Modified: 05/15/19
#
# Inputs:
#    x = fully observed proxy vector
#    y = partially observed binary outcome vector
#    phi = sensitivity parameter in [0,1]
#
# Returns:
#    rho_0 = biserial correlation estimate
#    muY_0 = estimated mean of Y for respondents
#    muY_1 = estimated mean of Y for nonrespondents
#    muY = estimated mean of Y (overall)
#####################################################
mle2step <- function(x, y, phi)
{
  # Indicator for missingness
  m <- ifelse(is.na(y), 1, 0)
  # Number of respondents, nonrespondents
  n <- length(y)
  r <- sum(1-m)
  
  # Separate out respondent and non-respondent vectors
  x_0 <- x[m==0]
  y_0 <- y[m==0]
  x_1 <- x[m==1]
  
  # Identified parameters
  # E[X|M=0]
  muX_0 <- mean(x_0)
  # Var[X|M=0]
  sigmaXX_0 <- sum((x_0 - muX_0)^2)/r
  # E[X|M=1]
  muX_1 <- mean(x_1)
  # Var[X|M=1]
  sigmaXX_1 <- sum((x_1 - muX_1)^2)/(n-r)
  # Pr(M=1)
  pi <- r/n
  
  # Polyserial correlation and threshold
  ## TWO-STEP METHOD
  # Cutpoint fixed at inverse probit of respondent sample mean
  muU_0 <- qnorm(mean(y_0))

  # Maximize likelihood wrt p, holding w constant
  # Likelihood containing (p)
  f <- function(pars)
  {
    p <- pars[1]
    a <- -(-muU_0 + muX_0*p/sqrt(sigmaXX_0))/sqrt(1-p^2)
    b <- (p/sqrt(sigmaXX_0))/sqrt(1-p^2)
    logPhi <- pnorm(a + b*x_0, log.p=TRUE)
    log1minusPhi <- pnorm(a + b*x_0, log.p=TRUE, lower.tail=FALSE)
    -sum(y_0*logPhi + (1-y_0)*log1minusPhi)
  }
  result <- optimize(f, interval=c(-0.99, 0.99))
  rho_0 <- result$minimum
  
  ## MLEs for distribution of U - for selected phi value
  g <- (phi+(1-phi)*rho_0)/(phi*rho_0+(1-phi))
  sigmaXU_0 <- rho_0*sqrt(sigmaXX_0)
  muU_1 <- muU_0 + g*(muX_1 - muX_0)/sqrt(sigmaXX_0)
  sigmaUU_1 <- 1 + (1/sigmaXX_0)*g^2*(sigmaXX_1-sigmaXX_0)
  sigmaXU_1 <- sigmaXU_0 + sqrt(1/sigmaXX_0)*g*(sigmaXX_1-sigmaXX_0)
  # If sigmaUU_1 < 0 replace with boundary value .Machine$double.eps
  # This will cause denominator of muY_1 estimate to be very small,
  # hence muY_1=pnorm(very large magnitude)=1 or -1 depending on sign of muU_1
  sigmaUU_1 <- ifelse(sigmaUU_1<0, .Machine$double.eps, sigmaUU_1)
  
  ## MLEs for E[Y|M=0], E[Y|M=1], E[Y]
  ### E[Y|M=0]
  muY_0 <- pnorm(muU_0)
  ### E[Y|M=1]
  muY_1 <- pnorm(muU_1/sqrt(sigmaUU_1))
  ### MU_Y
  muY <- pi*muY_0 + (1-pi)*muY_1

  return(list(rho_0=rho_0,
              muY_0=muY_0,
              muY_1=muY_1,
              muY=muY))
}

#####################################################
# Bayesian posterior draws (unmodified and both modifications)
# for a given phi in [0,1]
#
# Authors: Rebecca Andridge
# Last Modified: 05/15/19
#
# Inputs:
#    y = partially observed binary outcome vector
#    z = matrix of fully observed covariates **WITH NO INTERCEPT TERM**
#    phi = sensitivity parameter in [0,1] (defaults to 0.5)
#    drawphi = T/F draw phi from Uniform(0,1)
#              Defaults to FALSE
#              If TRUE then draws overwrite the specified phi value
#    nreps = number of draws (defaults to 1000)
#
# Returns:
#    dataframe w/matrix of drawn parameter values
#####################################################
require(mvtnorm)
proxyDraws <- function(y,z,phi=0.5,drawphi=FALSE,nreps=1000)
{
  ### Some functions to do truncated normals in R
  rnorm.lt <- function(n, lv=rep(0,n), mv=rep(0.5,n), sv=rep(1,n))
  {
    lstd <- (lv-mv)/sv 		
    sv*qnorm(runif(n,pnorm(lstd),rep(1,n)))+mv
  } 
  rnorm.rt <- function(n, rv=rep(0,n), mv=rep(0.5,n), sv=rep(1,n))
  {
    rstd <- (rv-mv)/sv 		
    sv*qnorm(runif(n,rep(0,n),pnorm(rstd)))+mv
  }
  # Indicator for missingness
  m <- ifelse(is.na(y), 1, 0)
  
  # Count total n and number of respondents
  n <- length(y)
  r <- n - sum(m)
  # Count where Y=1
  r1 <- sum(y, na.rm=T)
  
  # Make sure Z is a matrix
  if (!is.matrix(z)) z <- as.matrix(z)
  
  # Probit regression of Y|Z to find starting point for Gibbs sampler
  fit <- glm(y ~ z, family=binomial(link="probit"))
  betaHat <- as.vector(fit$coef)
  Z <- cbind(rep(1,nrow(z)), z)           # attaches column of 1s (R+NR)
  Zobs <- model.matrix(fit)               # attaches column of 1s (R only)
  ZobsTZobsInv <- solve(t(Zobs)%*%Zobs)   # so only have to calculate once
  
  # Starting point from probit fit
  B <- as.matrix(betaHat)
  
  # Starting proxy value
  X <- Z %*% B    # both R + NR
  
  # Initialize matrix to hold results
  draws <- matrix(nrow=nreps, ncol=19)
  
  # Start looping
  for (j in 1:nreps)
  {
    #    if ((j %% 100)==0) print(paste("Iteration",j))
    # Intialize latent U vector
    u <- rep(NA, n)
    
    ## (0) Draw phi from Beta prior with alpha=1, beta=1 if requested, and overwrite input value
    if (drawphi)
    {
      phi <- runif(1)
    }
    
    ## (1) Draw latent U|Y,B  (aka U|Y,X)
    u[y==1 & m==0] <- rnorm.lt(r1,   mv=X[y==1 & m==0]) # Draws for Y=1 --> truncated at left by 0
    u[y==0 & m==0] <- rnorm.rt(r-r1, mv=X[y==0 & m==0]) # Draws for Y=0 --> truncated at right by 0
    
    ## (2) Draw B|Z,U
    B <- t(rmvnorm(1, ZobsTZobsInv %*% t(Zobs) %*% u[m==0], ZobsTZobsInv))
    
    ## (3) Create proxy given current B
    X <- Z %*% B    # both R + NR
    
    ## (4) Scale proxy X to have same variance as latent U among respondents
    # Draw the population variances of X, Y* from posterior
    varXdraw <- sum((X[m==0] - mean(X[m==0]))^2) / rchisq(1, r-1)
    varUdraw <- sum((u[m==0] - mean(u[m==0]))^2) / rchisq(1, r-1)
    # Use draws to scale the proxy
    x <- X * sqrt(varUdraw/varXdraw)

    ## (5) Draw from PPM dependent on value of phi, using (X,U)
    if (phi==0){
      y1 <- x
      y2 <- u
      # Calculate means, sums of squares
      # Respondent data (M=0)
      y1Bar_0 <- mean(y1[m==0])
      y2Bar_0 <- mean(y2[m==0])
      s11_0 <- sum((y1[m==0] - y1Bar_0)^2)/r
      s22_0 <- sum((y2[m==0] - y2Bar_0)^2)/r
      s12_0 <- sum((y1[m==0] - y1Bar_0) * (y2[m==0] - y2Bar_0))/r
      b21.1_0 <- s12_0 / s11_0
      s22.1_0 <- s22_0 - (s12_0)^2 / s11_0
      # Nonrespondent data
      y1Bar_1 <- mean(y1[m==1])
      s11_1 <- sum((y1[m==1] - y1Bar_1)^2)/(n-r)
      # Draws
      # (1) PI
      PI <- rbeta(1, r + 0.5, n-r + 0.5)
      # (2) SIGMA11_0
      SIGMA11_0 <- r * s11_0 / rchisq(1, r-1)
      # (3) MU1_0 | SIGMA11_0
      MU1_0 <- rnorm(1, y1Bar_0, sqrt(SIGMA11_0/r))
      # (4) SIGMA22.1_0
      SIGMA22.1_0 <- r * s22.1_0 / rchisq(1, r-2)
      # (5) BETA21.1_0 | SIGMA22.1_0
      BETA21.1_0 <- rnorm(1, b21.1_0, sqrt(SIGMA22.1_0/(r*s11_0)))
      # (6) BETA20.1_0 | BETA21.1_0, SIGMA22.1_0
      BETA20.1_0 <- rnorm(1, y2Bar_0 - BETA21.1_0 * y1Bar_0, sqrt(SIGMA22.1_0/r))
      # (7) SIGMA11_1
      SIGMA11_1 <- (n-r) * s11_1 / rchisq(1, n-r-1)
      # (8) MU1_1 | SIGMA11_1
      MU1_1 <- rnorm(1, y1Bar_1, sqrt(SIGMA11_1/(n-r)))
      # Transform draws to get other parameters
      # (a) MU2_0
      MU2_0 <- BETA20.1_0 + BETA21.1_0 * MU1_0
      # (b) MU2_1
      MU2_1 <- BETA20.1_0 + BETA21.1_0 * MU1_1
      # (c) SIGMA12_0
      SIGMA12_0 <- BETA21.1_0 * SIGMA11_0
      # (d) SIGMA22_0
      SIGMA22_0 <- SIGMA22.1_0 + BETA21.1_0^2 * SIGMA11_0
      # (e) SIGMA22_1
      SIGMA22_1 <- SIGMA22.1_0 + BETA21.1_0^2 * SIGMA11_1
      # (f) SIGMA12_1
      SIGMA12_1 <- SIGMA11_1 * BETA21.1_0
      # All Draws
      drawsPPM <- list(pi=PI, mu1_0=MU1_0, mu2_0=MU2_0, mu1_1=MU1_1, mu2_1=MU2_1,
                       sigma11_0=SIGMA11_0, sigma12_0=SIGMA12_0, sigma22_0=SIGMA22_0,
                       sigma11_1=SIGMA11_1, sigma12_1=SIGMA12_1, sigma22_1=SIGMA22_1)      
    } else {
      if (phi==1){
        y1 <- x
        y2 <- u
      } else {
        y1 <- x
        y2 <- (1-phi)*x + phi*u
      }
      # Calculate means, sums of squares
      # Respondent data (M=0)
      y1Bar_0 <- mean(y1[m==0])
      y2Bar_0 <- mean(y2[m==0])
      s11_0 <- sum((y1[m==0] - y1Bar_0)^2)/r
      s22_0 <- sum((y2[m==0] - y2Bar_0)^2)/r
      s12_0 <- sum((y1[m==0] - y1Bar_0) * (y2[m==0] - y2Bar_0))/r
      b12.2_0 <- s12_0 / s22_0
      s11.2_0 <- s11_0 - (s12_0)^2 / s22_0
      # Nonrespondent data
      y1Bar_1 <- mean(y1[m==1])
      s11_1 <- sum((y1[m==1] - y1Bar_1)^2)/(n-r)
      # Draws
      # (1) PI
      PI <- rbeta(1, r + 0.5, n-r + 0.5)
      # (2) SIGMA22_0
      SIGMA22_0 <- r * s22_0 / rchisq(1, r-1)
      # (3) MU2_0 | SIGMA22_0
      MU2_0 <- rnorm(1, y2Bar_0, sqrt(SIGMA22_0/r))
      # (4, 5) SIGMA11.2_0, SIGMA11_1 with constraint
      goodDraw <- FALSE
      ct <- 1
      while (!goodDraw)
      {
        # Repeat these 2 draws until SIGMA11_1 > SIGMA11.2_0
        # (4) SIGMA11.2_0
        SIGMA11.2_0 <- r * s11.2_0 / rchisq(1, r-2)
        # (5) SIGMA11_1
        SIGMA11_1 <- (n-r) * s11_1 / rchisq(1, n-r-1)
        # Check to see if draws meet the condition
        goodDraw <- (SIGMA11_1 >= SIGMA11.2_0)
        if (ct > 20){
          goodDraw <- TRUE
          SIGMA11.2_0 <- SIGMA11_1
        }
        ct <- ct + 1
      }
      # (6) BETA12.2_0 | SIGMA11.2_0
      BETA12.2_0 <- rnorm(1, b12.2_0, sqrt(SIGMA11.2_0/(r*s22_0)))
      # (7) BETA10.2_0 | BETA12.2_0, SIGMA11.2_0
      BETA10.2_0 <- rnorm(1, y1Bar_0 - BETA12.2_0*y2Bar_0, sqrt(SIGMA11.2_0/r))
      # (8) MU1_1 | SIGMA11_1
      MU1_1 <- rnorm(1, y1Bar_1, sqrt(SIGMA11_1/(n-r)))
      # Transform draws to get other parameters
      # (a) MU2_1
      MU2_1 <- (MU1_1 - BETA10.2_0) / BETA12.2_0
      # (b) MU1_0
      MU1_0 <- BETA10.2_0 + BETA12.2_0 * MU2_0
      # (c) SIGMA12_0
      SIGMA12_0 <- BETA12.2_0 * SIGMA22_0
      # (d) SIGMA11_0
      SIGMA11_0 <- SIGMA11.2_0 + BETA12.2_0^2 * SIGMA22_0
      # (e) SIGMA22_1
      SIGMA22_1 <- (SIGMA11_1 - SIGMA11.2_0) / BETA12.2_0^2
      # (f) SIGMA12_1
      SIGMA12_1 <- SIGMA22_1 * BETA12.2_0
      # All Draws
      drawsPPM <- list(pi=PI, mu1_0=MU1_0, mu2_0=MU2_0, mu1_1=MU1_1, mu2_1=MU2_1,
                       sigma11_0=SIGMA11_0, sigma12_0=SIGMA12_0, sigma22_0=SIGMA22_0,
                       sigma11_1=SIGMA11_1, sigma12_1=SIGMA12_1, sigma22_1=SIGMA22_1)     
      
    }
    if (phi != 0 & phi != 1)
    {
      # Transform draws of [X,W] to get draws from [X,U]
      # W = (1-phi)*X + phi*U --> U = (W - (1-phi)*X)/phi
      # Start with draws of [X,W] and then overwrite parms relating to U
      drawsXW <- drawsPPM
      drawsPPM$mu2_0 <- (drawsXW$mu2_0 - (1-phi)*drawsXW$mu1_0)/phi
      drawsPPM$mu2_1 <- (drawsXW$mu2_1 - (1-phi)*drawsXW$mu1_1)/phi
      drawsPPM$sigma22_0 <- (drawsXW$sigma22_0 + (1-phi)^2*drawsXW$sigma11_0 - 2*(1-phi)*drawsXW$sigma12_0)/phi^2
      drawsPPM$sigma22_1 <- (drawsXW$sigma22_1 + (1-phi)^2*drawsXW$sigma11_1 - 2*(1-phi)*drawsXW$sigma12_1)/phi^2
      drawsPPM$sigma12_0 <- (drawsXW$sigma12_0 - (1-phi)*drawsXW$sigma11_0)/phi
      drawsPPM$sigma12_1 <- (drawsXW$sigma12_1 - (1-phi)*drawsXW$sigma11_1)/phi
    }
    
    ## (6) Draw U for nonrespondents given PPM parameters and current value of X (for modifications to PD)
    ## Calculate conditional mean and variance for [U|X,M=1]
    cmean <- drawsPPM$mu2_1 + (drawsPPM$sigma12_1/drawsPPM$sigma11_1)*(x[m==1] - drawsPPM$mu1_1)
    cvar <- drawsPPM$sigma22_1 - drawsPPM$sigma12_1^2/drawsPPM$sigma11_1
    ## Draw U for nonrespondents
    u[m==1] <- rnorm(n-r, mean=cmean, sd=sqrt(cvar))
    
    ## (7) Draw U for respondents given PPM parameters and current value of X (for modifications to PD)
    ## Calculate conditional mean and variance for [U|X,M=0]
    cmean <- drawsPPM$mu2_0 + (drawsPPM$sigma12_0/drawsPPM$sigma11_0)*(x[m==0] - drawsPPM$mu1_0)
    cvar <- drawsPPM$sigma22_0 - drawsPPM$sigma12_0^2/drawsPPM$sigma11_0
    ## Draw U for respondents
    u[m==0] <- rnorm(r, mean=cmean, sd=sqrt(cvar))
    
    # Three methods of obtaining draws of E[Y]:
    # (a) Unmodified method: transform the parameters
    # (b) Modification 1 (redraw):
    #            E[Y|M=1] = average of I(U>0) using draws of U for nonrespondents
    #            E[Y|M=0] = average of I(U>0) using draws of U for respondents ("redrawing" U, unconditional on Y)
    # (c) Modification 2 (predprob)
    #            E[Y|M=1] = average of I(U>0) using draws of U for nonrespondents
    #            E[Y|M=0] = average of predicted probabilities using draws of unscaled X for respondents
    #### Estimate of mean of Y|M=0
    # (a) Unmodified Bayes:
    drawsPPM$muY_0a <- pnorm(drawsPPM$mu2_0/sqrt(drawsPPM$sigma22_0))
    # (b) Mod 1 (redraw)
    drawsPPM$muY_0b <- mean(u[m==0]>0)
    # (c) Mod 2 (predprob)    
    drawsPPM$muY_0c <- mean(sample(pnorm(X[m==0]), replace=TRUE))  # note: uses UNSCALED proxy
    
    #### Estimate of mean of Y|M=1 (only 2 estimates, since (b) and (c) use same estimate)
    # (a) Unmodified Bayes:
    drawsPPM$muY_1a <- pnorm(drawsPPM$mu2_1/sqrt(drawsPPM$sigma22_1))
    # (b) Mod 1 (redraw) / (c) Mod 2 (predprob)    
    drawsPPM$muY_1bc <- mean(u[m==1]>0)
    
    #### Estimate of mean of Y
    # (a) Unmodified Bayes:
    drawsPPM$muY_a <- drawsPPM$pi*drawsPPM$muY_0a + (1-drawsPPM$pi)*drawsPPM$muY_1a
    # (b) Mod 1 (redraw)
    drawsPPM$muY_b <- drawsPPM$pi*drawsPPM$muY_0b + (1-drawsPPM$pi)*drawsPPM$muY_1bc
    # (c) Mod 2 (predprob)
    drawsPPM$muY_c <- drawsPPM$pi*drawsPPM$muY_0c + (1-drawsPPM$pi)*drawsPPM$muY_1bc
    
    # Save draws
    draws[j,] <- unlist(drawsPPM)
    
  }
  # End looping
  draws <- as.data.frame(draws)
  names(draws) <- names(unlist(drawsPPM))
  
  return(draws)
}


#####################################################
# Multiple Imputation for a given phi in [0,1]
#
# Authors: Rebecca Andridge
# Last Modified: 05/15/19
#
# Inputs:
#    y = partially observed binary outcome vector
#    z = matrix of fully observed covariates **WITH NO INTERCEPT TERM**
#    phi = sensitivity parameter in [0,1] (defaults to 0.5)
#    drawphi = T/F draw phi from Uniform(0,1)
#              Defaults to FALSE
#              If TRUE then draws overwrite the specified phi value
#    D = number of multiply imputed data sets (defaults to 10)
#    burnin = # of burn-in replicates (defaults to 500)
#    thin = # of replicates between imputations (defaults to 100)
#
# Returns:
#    matrix of multiply imputed Y vectors
#####################################################
mi <- function(y,z,phi=0.5,drawphi=FALSE,D=10,burnin=500,thin=100)
{
  ### Some functions to do truncated normals in R
  rnorm.lt <- function(n, lv=rep(0,n), mv=rep(0.5,n), sv=rep(1,n))
  {
    lstd <- (lv-mv)/sv 		
    sv*qnorm(runif(n,pnorm(lstd),rep(1,n)))+mv
  } 
  rnorm.rt <- function(n, rv=rep(0,n), mv=rep(0.5,n), sv=rep(1,n))
  {
    rstd <- (rv-mv)/sv 		
    sv*qnorm(runif(n,rep(0,n),pnorm(rstd)))+mv
  }
  # Indicator for missingness
  m <- ifelse(is.na(y), 1, 0)
  
  # Count total n and number of respondents
  n <- length(y)
  r <- n - sum(m)
  # Count where Y=1
  r1 <- sum(y, na.rm=T)
  
  # Make sure Z is a matrix
  if (!is.matrix(z)) z <- as.matrix(z)
  
  # Probit regression of Y|Z to find starting point for Gibbs sampler
  fit <- glm(y ~ z, family=binomial(link="probit"))
  betaHat <- as.vector(fit$coef)
  Z <- cbind(rep(1,nrow(z)), z)           # attaches column of 1s (R+NR)
  Zobs <- model.matrix(fit)               # attaches column of 1s (R only)
  ZobsTZobsInv <- solve(t(Zobs)%*%Zobs)   # so only have to calculate once
  
  # Starting point from probit fit
  B <- as.matrix(betaHat)
  
  # Starting proxy value
  X <- Z %*% B    # both R + NR
  
  # Initialize vector(matrix) to hold multiply-imputed Y
  imp <- matrix(nrow=n, ncol=D)
  imp[m==0,] <- y[m==0]  # Fill in observed values of Y
  
  # Total number of replicates needed is burnin+1+(D-1)*thin
  # Draw U, B on all replicates
  # Only draw PPM parameters on the replicates to perform MI
  # Keep track of which iterate you're on
  d=0
  # Burn-in iterations for drawing (U,B)
  for (i in 1:(burnin+1+(D-1)*thin))
  {
    #print(i)
    # Intialize latent U vector
    u <- rep(NA, n)
    
    ## (0) Draw phi from Beta prior with alpha=1, beta=1 if requested, and overwrite input value
    if (drawphi)
    {
      phi <- runif(1)
    }
    
    ## (1) Draw latent U|Y,B  (aka U|Y,X)
    u[y==1 & m==0] <- rnorm.lt(r1,   mv=X[y==1 & m==0]) # Draws for Y=1 --> truncated at left by 0
    u[y==0 & m==0] <- rnorm.rt(r-r1, mv=X[y==0 & m==0]) # Draws for Y=0 --> truncated at right by 0
    
    ## (2) Draw B|Z,U
    B <- t(rmvnorm(1, ZobsTZobsInv %*% t(Zobs) %*% u[m==0], ZobsTZobsInv))
    
    ## (3) Create proxy given current B
    X <- Z %*% B    # both R + NR
    
    ##### IF REPLICATE IS TO BE USED FOR IMPUTATION ####
    # Proceed with drawing PPM parameters conditional on drawn U,B
    if ((i>burnin) & (((i - 1 - burnin) %% thin)==0))
    {
      d = d + 1
      #print(c("INNER LOOP :",d))
      ############ Step 1 ############
      
      ## (4) Scale proxy X to have same variance as latent U among respondents
      # Draw the population variances of X, Y* from posterior
      varXdraw <- sum((X[m==0] - mean(X[m==0]))^2) / rchisq(1, r-1)
      varUdraw <- sum((u[m==0] - mean(u[m==0]))^2) / rchisq(1, r-1)
      # Use draws to scale the proxy
      x <- X * sqrt(varUdraw/varXdraw)

      ## (5) Draw from PPM dependent on value of phi, using (X,U)
      if (phi==0){
        y1 <- x
        y2 <- u
        # Calculate means, sums of squares
        # Respondent data (M=0)
        y1Bar_0 <- mean(y1[m==0])
        y2Bar_0 <- mean(y2[m==0])
        s11_0 <- sum((y1[m==0] - y1Bar_0)^2)/r
        s22_0 <- sum((y2[m==0] - y2Bar_0)^2)/r
        s12_0 <- sum((y1[m==0] - y1Bar_0) * (y2[m==0] - y2Bar_0))/r
        b21.1_0 <- s12_0 / s11_0
        s22.1_0 <- s22_0 - (s12_0)^2 / s11_0
        # Nonrespondent data
        y1Bar_1 <- mean(y1[m==1])
        s11_1 <- sum((y1[m==1] - y1Bar_1)^2)/(n-r)
        # Draws
        # (1) PI
        PI <- rbeta(1, r + 0.5, n-r + 0.5)
        # (2) SIGMA11_0
        SIGMA11_0 <- r * s11_0 / rchisq(1, r-1)
        # (3) MU1_0 | SIGMA11_0
        MU1_0 <- rnorm(1, y1Bar_0, sqrt(SIGMA11_0/r))
        # (4) SIGMA22.1_0
        SIGMA22.1_0 <- r * s22.1_0 / rchisq(1, r-2)
        # (5) BETA21.1_0 | SIGMA22.1_0
        BETA21.1_0 <- rnorm(1, b21.1_0, sqrt(SIGMA22.1_0/(r*s11_0)))
        # (6) BETA20.1_0 | BETA21.1_0, SIGMA22.1_0
        BETA20.1_0 <- rnorm(1, y2Bar_0 - BETA21.1_0 * y1Bar_0, sqrt(SIGMA22.1_0/r))
        # (7) SIGMA11_1
        SIGMA11_1 <- (n-r) * s11_1 / rchisq(1, n-r-1)
        # (8) MU1_1 | SIGMA11_1
        MU1_1 <- rnorm(1, y1Bar_1, sqrt(SIGMA11_1/(n-r)))
        # Transform draws to get other parameters
        # (a) MU2_0
        MU2_0 <- BETA20.1_0 + BETA21.1_0 * MU1_0
        # (b) MU2_1
        MU2_1 <- BETA20.1_0 + BETA21.1_0 * MU1_1
        # (c) SIGMA12_0
        SIGMA12_0 <- BETA21.1_0 * SIGMA11_0
        # (d) SIGMA22_0
        SIGMA22_0 <- SIGMA22.1_0 + BETA21.1_0^2 * SIGMA11_0
        # (e) SIGMA22_1
        SIGMA22_1 <- SIGMA22.1_0 + BETA21.1_0^2 * SIGMA11_1
        # (f) SIGMA12_1
        SIGMA12_1 <- SIGMA11_1 * BETA21.1_0
        # All Draws
        drawsPPM <- list(pi=PI, mu1_0=MU1_0, mu2_0=MU2_0, mu1_1=MU1_1, mu2_1=MU2_1,
                         sigma11_0=SIGMA11_0, sigma12_0=SIGMA12_0, sigma22_0=SIGMA22_0,
                         sigma11_1=SIGMA11_1, sigma12_1=SIGMA12_1, sigma22_1=SIGMA22_1)      
      } else {
        if (phi==1){
          y1 <- x
          y2 <- u
        } else {
          y1 <- x
          y2 <- (1-phi)*x + phi*u
        }
        # Calculate means, sums of squares
        # Respondent data (M=0)
        y1Bar_0 <- mean(y1[m==0])
        y2Bar_0 <- mean(y2[m==0])
        s11_0 <- sum((y1[m==0] - y1Bar_0)^2)/r
        s22_0 <- sum((y2[m==0] - y2Bar_0)^2)/r
        s12_0 <- sum((y1[m==0] - y1Bar_0) * (y2[m==0] - y2Bar_0))/r
        b12.2_0 <- s12_0 / s22_0
        s11.2_0 <- s11_0 - (s12_0)^2 / s22_0
        # Nonrespondent data
        y1Bar_1 <- mean(y1[m==1])
        s11_1 <- sum((y1[m==1] - y1Bar_1)^2)/(n-r)
        # Draws
        # (1) PI
        PI <- rbeta(1, r + 0.5, n-r + 0.5)
        # (2) SIGMA22_0
        SIGMA22_0 <- r * s22_0 / rchisq(1, r-1)
        # (3) MU2_0 | SIGMA22_0
        MU2_0 <- rnorm(1, y2Bar_0, sqrt(SIGMA22_0/r))
        # (4, 5) SIGMA11.2_0, SIGMA11_1 with constraint
        goodDraw <- FALSE
        while (!goodDraw)
        {
          # Repeat these 2 draws until SIGMA11_1 > SIGMA11.2_0
          # (4) SIGMA11.2_0
          SIGMA11.2_0 <- r * s11.2_0 / rchisq(1, r-2)
          # (5) SIGMA11_1
          SIGMA11_1 <- (n-r) * s11_1 / rchisq(1, n-r-1)
          # Check to see if draws meet the condition
          goodDraw <- (SIGMA11_1 > SIGMA11.2_0)
        }
        # (6) BETA12.2_0 | SIGMA11.2_0
        BETA12.2_0 <- rnorm(1, b12.2_0, sqrt(SIGMA11.2_0/(r*s22_0)))
        # (7) BETA10.2_0 | BETA12.2_0, SIGMA11.2_0
        BETA10.2_0 <- rnorm(1, y1Bar_0 - BETA12.2_0*y2Bar_0, sqrt(SIGMA11.2_0/r))
        # (8) MU1_1 | SIGMA11_1
        MU1_1 <- rnorm(1, y1Bar_1, sqrt(SIGMA11_1/(n-r)))
        # Transform draws to get other parameters
        # (a) MU2_1
        MU2_1 <- (MU1_1 - BETA10.2_0) / BETA12.2_0
        # (b) MU1_0
        MU1_0 <- BETA10.2_0 + BETA12.2_0 * MU2_0
        # (c) SIGMA12_0
        SIGMA12_0 <- BETA12.2_0 * SIGMA22_0
        # (d) SIGMA11_0
        SIGMA11_0 <- SIGMA11.2_0 + BETA12.2_0^2 * SIGMA22_0
        # (e) SIGMA22_1
        SIGMA22_1 <- (SIGMA11_1 - SIGMA11.2_0) / BETA12.2_0^2
        # (f) SIGMA12_1
        SIGMA12_1 <- SIGMA22_1 * BETA12.2_0
        # All Draws
        drawsPPM <- list(pi=PI, mu1_0=MU1_0, mu2_0=MU2_0, mu1_1=MU1_1, mu2_1=MU2_1,
                         sigma11_0=SIGMA11_0, sigma12_0=SIGMA12_0, sigma22_0=SIGMA22_0,
                         sigma11_1=SIGMA11_1, sigma12_1=SIGMA12_1, sigma22_1=SIGMA22_1)     
        
      }
      if (phi != 0 & phi != 1)
      {
        # Transform draws of [X,W] to get draws from [X,U]
        # W = (1-phi)*X + phi*U --> U = (W - (1-phi)*X)/phi
        # Start with draws of [X,W] and then overwrite parms relating to U
        drawsXW <- drawsPPM
        drawsPPM$mu2_0 <- (drawsXW$mu2_0 - (1-phi)*drawsXW$mu1_0)/phi
        drawsPPM$mu2_1 <- (drawsXW$mu2_1 - (1-phi)*drawsXW$mu1_1)/phi
        drawsPPM$sigma22_0 <- (drawsXW$sigma22_0 + (1-phi)^2*drawsXW$sigma11_0 - 2*(1-phi)*drawsXW$sigma12_0)/phi^2
        drawsPPM$sigma22_1 <- (drawsXW$sigma22_1 + (1-phi)^2*drawsXW$sigma11_1 - 2*(1-phi)*drawsXW$sigma12_1)/phi^2
        drawsPPM$sigma12_0 <- (drawsXW$sigma12_0 - (1-phi)*drawsXW$sigma11_0)/phi
        drawsPPM$sigma12_1 <- (drawsXW$sigma12_1 - (1-phi)*drawsXW$sigma11_1)/phi
      }
      
      ############ Step 2 ############
      
      ## (6) Draw U for nonrespondents given PPM parameters and current value of X
      ### Calculate conditional mean and variance for [U|X,M=1]
      cmean <- drawsPPM$mu2_1 + (drawsPPM$sigma12_1/drawsPPM$sigma11_1)*(x[m==1] - drawsPPM$mu1_1)
      cvar <- drawsPPM$sigma22_1 - drawsPPM$sigma12_1^2/drawsPPM$sigma11_1
      ### Draw U for nonrespondents
      u[m==1] <- rnorm(n-r, mean=cmean, sd=sqrt(cvar))
      ### imputed Y=1 if U>0 ##
      imp[m==1, d] <- (u[m==1]>0)
    }
    
  }
  # End looping
  
  return(imp)
  
}
