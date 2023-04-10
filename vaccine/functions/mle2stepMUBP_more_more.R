### VERSION OF FUNCTION THAT RETURNS MORE VALUES (in addition to the MUBP)
#  phi
#  rho_0
#  mubp (MUBP for chosen phi value)
#  ymean (estimated overall mean of Y)
#  ymean_0 (selected sample mean of Y)
#  ymean_1 (estimated nonselected sample mean of Y)
#  xmean_0 (selected sample mean of X)
#  xmean_1 (nonselected sample mean of X)
#### NOTE: Notation throughout code uses 0=sampled(not missing), 1=not sampled(missing)
mle2stepMUBP_more_more <- function(x_0, y_0, xmean_1, xvar_1, sfrac, phi, verbose = F)
{
  
  # Proxy vector, mean, and variance for selected sample
  xmean_0 <- mean(x_0)
  xvar_0 <- sum((x_0-xmean_0)^2)/length(x_0)
  
  # Polyserial correlation and threshold
  # TWO-STEP METHOD
  # Cutpoint fixed
  w <- qnorm(1-mean(y_0))  
  # Maximize likelihood wrt p, holding w constant
  # Likelihood containing (p)
  f <- function(pars)
  {
    p <- pars[1]
    a <- -(w + xmean_0*p/sqrt(xvar_0))/sqrt(1-p^2)
    b <- (p/sqrt(xvar_0))/sqrt(1-p^2)
    logPhi <- pnorm(a + b*x_0, log.p = T)
    log1minusPhi <- pnorm(a + b*x_0, log.p = T, lower.tail = F)
    -sum(y_0*logPhi + (1-y_0)*log1minusPhi)
  }
  result <- optimize(f, interval=c(-0.99, 0.99))
  rho_0 <- result$minimum
  if(verbose) {
    cat("Two-Step Biserial Correlation: ",rho_0,"\n") # two-step biserial correlation
  }
  
  # MLEs for distribution of U
  if (phi==1) {
    g <- 1/rho_0
  } else {
    g <- (phi+(1-phi)*rho_0)/(phi*rho_0+(1-phi))
  }
  umean_0 <- -w
  uvar_0 <- 1
  xucov_0 <- rho_0*sqrt(xvar_0)
  umean_1 <- umean_0 + g*(xmean_1 - xmean_0)/sqrt(xvar_0)
  uvar_1 <- 1 + g^2*(xvar_1-xvar_0)/xvar_0
  # If uvar_1 < 0 replace with boundary value .Machine$double.eps
  # This will cause pnorm(umean_1/sqrt(uvar_1)) = +1 or -1 depending on sign of umean_1
  uvar_1 <- ifelse(uvar_1<0, .Machine$double.eps, uvar_1)
  xucov_1 <- xucov_0 + sqrt(xvar_0)*g*(xvar_1-xvar_0)/xvar_0
  rho_1 <- xucov_1/sqrt(xvar_1*uvar_1)
  uparms <- list(xmean_0=xmean_0, xvar_0=xvar_0, umean_0=umean_0, uvar_0=uvar_0, xucov_0=xucov_0, rho_0=rho_0, 
                 xmean_1=xmean_1, xvar_1=xvar_1, umean_1=umean_1, uvar_1=uvar_1, xucov_1=xucov_1, rho_1=rho_1)
  # MLEs for regression of U|X
  # intercepts
  beta_0 <- umean_0 - xucov_0/xvar_0*xmean_0
  beta_1 <- umean_1 - xucov_1/xvar_1*xmean_1
  # slopes
  alpha_0 <- xucov_0/xvar_0 
  alpha_1 <- xucov_1/xvar_1
  # cond'l variances
  condvar_0 <- (1-rho_0^2)*uvar_0
  condvar_1 <- (1-rho_1^2)*uvar_1
  
  # MLEs for distribution of Y and the MSB
  ## E[Y|M=0]
  ymean_0 <- mean(y_0)  # same as pnorm(umean_0) b/c 2-step estimator
  ## E[Y|M=1]
  ymean_1 <- pnorm(umean_1/sqrt(uvar_1))
  ## E[Y]
  ymean <- sfrac*ymean_0 + (1-sfrac)*ymean_1
  ## MUBP(phi)
  if(verbose) {
    cat("MUBP(",phi,"):",ymean_0-ymean,"\n");
  }
  
  # Parameters of corresponding selection model
  # Intercept
  gamma_int <- log((1-sfrac)/sfrac) + xmean_0^2/(2*xvar_0) - xmean_1^2/(2*xvar_1) + 1/2*log(xvar_0/xvar_1) + beta_0^2/(2*condvar_0) - beta_1^2/(2*condvar_1) + 1/2*log(condvar_0/condvar_1)
  # X effect
  gamma_x <- xmean_1/xvar_1 - xmean_0/xvar_0 + beta_0*alpha_0/condvar_0 - beta_1*alpha_1/condvar_1
  # X^2 effect
  gamma_x2 <- 1/(2*xvar_0) - 1/(2*xvar_1) + alpha_0^2/(2*condvar_0) - alpha_1^2/(2*condvar_1)
  # U effect
  gamma_u <- beta_1/condvar_1 - beta_0/condvar_0
  # U^2 effect
  gamma_u2 <- 1/(2*condvar_0) - 1/(2*condvar_1)
  # XY effect
  gamma_xu <- alpha_1/condvar_1 - alpha_0/condvar_0
  # list
  sm <- list(int=gamma_int, x=gamma_x, x2=gamma_x2, u=gamma_u, u2=gamma_u2, xu=gamma_xu)
  
  return(list(results=list(phi=phi, rho_0=rho_0, mubp=ymean_0-ymean, 
                           ymean=ymean, ymean_0=ymean_0, ymean_1=ymean_1, 
                           xmean_0=xmean_0, xmean_1=xmean_1),
              uparms=uparms,
              sm=sm))
}