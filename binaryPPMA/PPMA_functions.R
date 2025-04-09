####################################################################
# Functions for Proxy Pattern-Mixture Analysis (Continuous Outcomes)
# Author: Rebecca Andridge (andridge.1@osu.edu)
# Last Modified: 2019.11.19
####################################################################

#####################################################
# MLE for mean of Y for given lambda value
#
# Inputs:
#    x = fully observed proxy vector
#    y = partially observed outcome vector
#    lambda = sensitivity parameter in [0,infinity)
#
# Returns:
#    muY = estimated mean of Y (overall)
#    muYvar = estimated variance of the mean of Y (overall)
#####################################################
mle <- function(x, y, lambda)
{
  # Indicator for missingness
  m <- ifelse(is.na(y), 1, 0)
  
  # Number of respondents, nonrespondents
  n <- length(x)
  r <- sum(1-m)
  
  # Calculate means, sums of squares
  # Respondent data
  xBar_0 <- mean(x[m==0])
  yBar_0 <- mean(y[m==0])
  sxx_0 <- sum((x[m==0] - xBar_0)^2)/r
  syy_0 <- sum((y[m==0] - yBar_0)^2)/r
  sxy_0 <- sum((x[m==0] - xBar_0) * (y[m==0] - yBar_0))/r
  byx.x_0 <- sxy_0 / sxx_0
  by0.x_0 <- yBar_0 - byx.x_0*xBar_0
  syy.x_0 <- syy_0 - (sxy_0)^2 / sxx_0
  bxy.y_0 <- sxy_0 / syy_0
  sxx.y_0 <- sxx_0 - (sxy_0)^2 / syy_0
  rhoHat_0 <- sxy_0/sqrt(sxx_0*syy_0)
  # Nonrespondent data
  xBar_1 <- mean(x[m==1])
  sxx_1 <- sum((x[m==1] - xBar_1)^2)/(n-r)
  
  # Value of g_lambda
  if (lambda==Inf) {
    g.lambda <- 1/bxy.y_0
  } else {
    g.lambda <- sqrt(syy_0/sxx_0)*(lambda + rhoHat_0)/(lambda*rhoHat_0 + 1)
    # Note that this reduces to byx.x_0=sxy_0/sxx_0 when lambda=0
  }
  
  # Variance of g_lambda
  if (lambda==Inf) {
    num <- (sxx_0*syy_0 - sxy_0^2)*syy_0^2
    denom <- r*sxy_0^4
    g.lambda.var <- num/denom
  } else {
    a <- sxx_0^2*syy_0^2*(1-lambda^2+lambda^4)
    b <- 2*sxx_0*syy_0*sxy_0*lambda*(3*lambda*sxy_0 + sqrt(sxx_0*syy_0)*(1+lambda^2))
    c <- lambda*sxy_0^3*(lambda*sxy_0 + 2*sqrt(sxx_0*syy_0)*(1+lambda^2))
    num <- (sxx_0*syy_0 - sxy_0^2)*(a + b + c)
    denom <- r*sxx_0^2*(sqrt(sxx_0*syy_0) + lambda*sxy_0)^4
    g.lambda.var <- num/denom
  }
  
  # MLE of E(x)
  muX <- mean(x)
  
  # MLE of Var(x)
  sigmaXX <- (r/n)*sxx_0 + ((n-r)/n)*sxx_1 + (r/n)*((n-r)/n)*(xBar_0 - xBar_1)^2
  
  # MLE of E(y)
  muY <- yBar_0 + g.lambda*(muX - xBar_0)
  
  # MLE of Var(y)
  sigmaYY <- syy_0 + (g.lambda^2)*(sigmaXX - sxx_0)
  
  # MLE of Cov(x,y)
  sigmaXY <- sxy_0 + g.lambda*(sigmaXX - sxx_0)
  
  # Variance of MLE of E(y)
  one <- sigmaYY/n
  two <- g.lambda.var*(muX - xBar_0)^2
  three <- ((n-r)/(r*n)) *(syy_0 - 2*g.lambda*sxy_0 + g.lambda^2*sxx_0)
  muYvar <- one + two + three
  
  return(list(muY=muY,
              muYvar=muYvar))
}

#####################################################
# Bayesian posterior draws
# for a given lambda in [0,Infinity)
#
# Inputs:
#    y = partially observed outcome vector
#    z = matrix of fully observed covariates **WITH NO INTERCEPT TERM**
#    lambda = sensitivity parameter in [0,infinity)
#    nreps = number of draws
#
# Returns:
#    muY = vector of drawn means of outcome Y
#####################################################
require(mvtnorm)
proxyDraws <- function(y,z,lambda,nreps)
{
  # Indicator for missingness
  m <- ifelse(is.na(y), 1, 0)
  
  # Count total n and number of respondents
  n <- length(y)
  r <- n - sum(m)
  
  # Make sure X is a matrix
  if (!is.matrix(z)) z <- as.matrix(z)
  
  # Regression of Y|X
  fit <- lm(y ~ z)
  betaHat <- as.vector(fit$coef)
  Vb <- summary(fit)$cov.unscaled
  s2 <- summary(fit)$sigma^2
  dfResid <- fit$df.residual
  
  # Initialize vector(matrix) to hold draws
  draws <- vector()
  
  # Loop through draws
  for (rep in 1:nreps)
  {
    ## (1) Draw SIGMA^2 from inverse-ChiSq
    SIGMA2 <- dfResid * s2 / rchisq(1, dfResid)
      
    ## (2) Draw BETA | SIGMA, Y
    BETA <- rmvnorm(1, betaHat, Vb*SIGMA2)
      
    ## (3) Calculate proxy x from these drawn parameters
    x <- cbind(rep(1,n),z) %*% t(BETA)

    ## (4) Scale proxy x to have same variance as y respondents
    # Draw the population variances of Y, Y* from posterior
    varYdraw <- sum((y[m==0] - mean(y[m==0]))^2) / rchisq(1, r-1)
    varXdraw <- sum((x[m==0] - mean(x[m==0]))^2) / rchisq(1, r-1)
    # Use draws to scale the proxy x
    x <- x * sqrt(varYdraw/varXdraw)

    ## (5) Draw from PPM dependent on value of lambda
    if (lambda==0){
      y1 <- x
      y2 <- y
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
      if (lambda==Inf){
        y1 <- x
        y2 <- y
      } else {
        y1 <- x
        y2 <- x + lambda*y   # Create W = X + LAMBDAxY
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
    
    # Transform parameter draws to get draws of E(Y)
    if (lambda==0 || lambda==Inf)
    {
      muY <- drawsPPM$pi*drawsPPM$mu2_0 + (1-drawsPPM$pi)*drawsPPM$mu2_1
    } else {
      # Transform draws of [X,W] to get draws from [X,Y]
      # W = X + LAMBDA*Y --> Y = (W-X)/LAMBDA
      muX <- drawsPPM$pi*drawsPPM$mu1_0 + (1-drawsPPM$pi)*drawsPPM$mu1_1
      muW <- drawsPPM$pi*drawsPPM$mu2_0 + (1-drawsPPM$pi)*drawsPPM$mu2_1
      muY <- (muW - muX)/lambda
    }
    
    # Save draws
    draws[rep] <- muY
  }
  return(draws)
}

#####################################################
# Multiple imputation
# for a given lambda in [0,Infinity)
#
# Inputs:
#    y = partially observed outcome vector
#    z = matrix of fully observed covariates **WITH NO INTERCEPT TERM**
#    lambda = sensitivity parameter in [0,infinity)
#    D = number of imputations
#
# Returns:
#    matrix of multiply imputed Y vectors
#####################################################
require(mvtnorm)
mi <- function(y,z,lambda,D)
{
  # Indicator for missingness
  m <- ifelse(is.na(y), 1, 0)
  
  # Count total n and number of respondents
  n <- length(y)
  r <- n - sum(m)
  
  # Make sure X is a matrix
  if (!is.matrix(z)) z <- as.matrix(z)
  
  # Regression of Y|X
  fit <- lm(y ~ z)
  betaHat <- as.vector(fit$coef)
  Vb <- summary(fit)$cov.unscaled
  s2 <- summary(fit)$sigma^2
  dfResid <- fit$df.residual
  
  # Initialize vector(matrix) to hold multiply-imputed Y
  imp <- vector()
  
  # Loop through number of MI data sets
  for (rep in 1:D)
  {
    ############ Step 1 ############
    # Draw parameter values from posterior dependent on value of lambda
    
    ## (1) Draw SIGMA^2 from inverse-ChiSq
    SIGMA2 <- dfResid * s2 / rchisq(1, dfResid)
    
    ## (2) Draw BETA | SIGMA, Y
    BETA <- rmvnorm(1, betaHat, Vb*SIGMA2)
    
    ## (3) Calculate proxy x from these drawn parameters
    x <- cbind(rep(1,n),z) %*% t(BETA)
    
    ## (4) Scale proxy x to have same variance as y respondents
    # Draw the population variances of Y, Y* from posterior
    varYdraw <- sum((y[m==0] - mean(y[m==0]))^2) / rchisq(1, r-1)
    varXdraw <- sum((x[m==0] - mean(x[m==0]))^2) / rchisq(1, r-1)
    # Use draws to scale the proxy x
    x <- x * sqrt(varYdraw/varXdraw)
    
    ## (5) Draw from PPM dependent on value of lambda
    if (lambda==0){
      y1 <- x
      y2 <- y
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
      if (lambda==Inf){
        y1 <- x
        y2 <- y
      } else {
        y1 <- x
        y2 <- x + lambda*y   # Create W = X + LAMBDAxY
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

    ############ Step 2 ############
    # Get parameters of conditional dist'n of Y|X,M=1
    if (lambda==0 || lambda==Inf)
    {
      cmean <- drawsPPM$mu2_1 + (drawsPPM$sigma12_1/drawsPPM$sigma11_1)*(x[m==1] - drawsPPM$mu1_1)
      cvar <- drawsPPM$sigma22_1 - (drawsPPM$sigma12_1^2)/drawsPPM$sigma11_1
    } else {
      # Transform draws of [X,W] to get draws from [X,Y]
      # W = X + LAMBDA*Y --> Y = (W-X)/LAMBDA
      mu2_1 <- (drawsPPM$mu2_1 - drawsPPM$mu1_1) / lambda
      sigma12_1 <- (1/lambda) * (drawsPPM$sigma12_1 - drawsPPM$sigma11_1)
      sigma22_1 <- (1/lambda^2) * (drawsPPM$sigma22_1 + drawsPPM$sigma11_1 - 2*drawsPPM$sigma12_1)
      cmean <- mu2_1 + (sigma12_1/drawsPPM$sigma11_1)*(x[m==1] - drawsPPM$mu1_1)
      cvar <- sigma22_1 - (sigma12_1^2)/drawsPPM$sigma11_1
    }
    # Impute Y from conditional dist'n of Y|X
    yimp <- rnorm(n-r, mean=cmean, sd=sqrt(cvar))
    imp <- cbind(imp, y)
    imp[m==1, rep] <- yimp
  }
  colnames(imp) <- NULL
  return(imp)
}
