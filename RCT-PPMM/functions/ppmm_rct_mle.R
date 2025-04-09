# Inputs:
# mean_YZ_s_A1 = vector of means of (Y,Z|S=1,A=1) (selected sample, treatment)
# var_YZ_s_A1 = variance matrix of (Y,Z|S=1,A=1) (selected sample, treatment)
# mean_YZ_s_A0 = vector of means of (Y,Z|S=1,A=0) (selected sample, control)
# var_YZ_s_A0 = variance matrix of (Y,Z|S=1,A=0) (selected sample, control)
# mean_Z_ns = vector of means of (Z|S=0) (non-selected sample)
# var_Z_ns = variance matrix of (Z|S=0) (non-selected sample)
# PHI_A1 = sensitivity parameter phi for treatment group (A=1) - can be a vector
# PHI_A0 = sensitivity parameter phi for treatment group (A=0) - can be a vector

ppmm_rct_mle <- function(mean_YZ_s_A1, var_YZ_s_A1,
                         mean_YZ_s_A0, var_YZ_s_A0,
                         mean_Z_ns, var_Z_ns,
                         PHI_A1, PHI_A0)
{
  # Separate selected sum stats for ease of use -----
  ## Treatment -----
  mean_Y_s_A1 <- mean_YZ_s_A1[1]
  mean_Z_s_A1 <- mean_YZ_s_A1[-1]
  var_Y_s_A1 <- var_YZ_s_A1[1,1]
  var_Z_s_A1 <- var_YZ_s_A1[-1,-1]
  covar_YZ_s_A1 <- var_YZ_s_A1[1,-1]
  ## Control -----
  mean_Y_s_A0 <- mean_YZ_s_A0[1]
  mean_Z_s_A0 <- mean_YZ_s_A0[-1]
  var_Y_s_A0 <- var_YZ_s_A0[1,1]
  var_Z_s_A0 <- var_YZ_s_A0[-1,-1]
  covar_YZ_s_A0 <- var_YZ_s_A0[1,-1]
  
  # Calculate proxy X means and variances -----
  ## Selected sample, treatment -----
  beta_YZ.Z_s_A1 <- drop(covar_YZ_s_A1 %*% solve(var_Z_s_A1)) # slopes
  beta_YZ.Z_s_A1 <- c(mean_Y_s_A1 - beta_YZ.Z_s_A1 %*% mean_Z_s_A1, beta_YZ.Z_s_A1) # add intercept
  mean_X_s_A1 <- drop(beta_YZ.Z_s_A1 %*% c(1, mean_Z_s_A1)) # proxy mean
  var_X_s_A1 <- drop(tcrossprod(beta_YZ.Z_s_A1[-1], var_Z_s_A1) %*% beta_YZ.Z_s_A1[-1]) # proxy variance
  ## Selected sample, control -----
  beta_YZ.Z_s_A0 <- drop(covar_YZ_s_A0 %*% solve(var_Z_s_A0)) # slopes
  beta_YZ.Z_s_A0 <- c(mean_Y_s_A0 - beta_YZ.Z_s_A0 %*% mean_Z_s_A0, beta_YZ.Z_s_A0) # add intercept
  mean_X_s_A0 <- drop(beta_YZ.Z_s_A0 %*% c(1, mean_Z_s_A0)) # proxy mean
  var_X_s_A0 <- drop(tcrossprod(beta_YZ.Z_s_A0[-1], var_Z_s_A0) %*% beta_YZ.Z_s_A0[-1]) # proxy variance
  ## Non-selected sample, under treatment -----
  mean_X_ns_A1 <- drop(beta_YZ.Z_s_A1 %*% c(1, mean_Z_ns)) # proxy mean
  var_X_ns_A1 <- drop(tcrossprod(beta_YZ.Z_s_A1[-1],var_Z_ns) %*% beta_YZ.Z_s_A1[-1]) # proxy variance
  ## Non-selected sample, under treatment -----
  mean_X_ns_A0 <- drop(beta_YZ.Z_s_A0 %*% c(1, mean_Z_ns)) # proxy mean
  var_X_ns_A0 <- drop(tcrossprod(beta_YZ.Z_s_A0[-1],var_Z_ns) %*% beta_YZ.Z_s_A0[-1]) # proxy variance
  
  # Correlation between proxy X and outcome Y -----
  cor_XY_s_A1 <- drop(beta_YZ.Z_s_A1[-1] %*% covar_YZ_s_A1 / sqrt(var_X_s_A1 * var_Y_s_A1))
  cor_XY_s_A0 <- drop(beta_YZ.Z_s_A0[-1] %*% covar_YZ_s_A0 / sqrt(var_X_s_A0 * var_Y_s_A0))
  
  # Mean of Y for non-selected for given phi values -----
  mean_Y_ns_A1 <- mean_Y_s_A1 + 
    sqrt(var_Y_s_A1/var_X_s_A1) * ((PHI_A1 + (1-PHI_A1)*cor_XY_s_A1)/(PHI_A1*cor_XY_s_A1 + 1 - PHI_A1)) * (mean_X_ns_A1 - mean_X_s_A1)
  mean_Y_ns_A0 <- mean_Y_s_A0 + 
    sqrt(var_Y_s_A0/var_X_s_A0) * ((PHI_A0 + (1-PHI_A0)*cor_XY_s_A0)/(PHI_A0*cor_XY_s_A0 + 1 - PHI_A0)) * (mean_X_ns_A0 - mean_X_s_A0)
  phis <- expand_grid(phi_A1 = PHI_A1, phi_A0 = PHI_A0)
  means <- expand_grid(mean_Y_ns_A1, mean_Y_ns_A0)
  parms <- tibble(rho_XY_s_A1 = cor_XY_s_A1, rho_XY_s_A0 = cor_XY_s_A0, mean_X_s_A1, mean_X_ns_A1, mean_X_s_A0, mean_X_ns_A0)
  results <- bind_cols(phis, means) %>%
    mutate(trt_effect_ns = mean_Y_ns_A1 - mean_Y_ns_A0)
  return(list(results=results, proxyinfo = parms))
}