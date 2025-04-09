######### Bayesian ############
require(MCMCpack)
require(mnormt)

means_bayes <- function(data_YZ, #  selected/non-missing sample
                            data_Z, # non-selected/missing sample
                            phi=0,# if phi_character=NULL, fix phi to a specific value in [0,1]    
                            phi_character=NULL, # draw phi from specific distribution (i.e."runif(1)"), the value of phi input to the function is ignored; 
                            # if phi_character=NULL, fix at value phi
                            ndraws=1500) # number of draws
{
  sumry_list = ErrorCheck(data_YZ, data_Z, prop_check = FALSE)
  data_YZ = sumry_list$data_YZ
  sumry_YZ = sumry_list$sumry_YZ
  sumry_Z = sumry_list$sumry_Z
  
  # Selected sample size
  n_s <- sumry_YZ$n_YZ
  # Non-Selected sample size
  n_ns <- sumry_Z$n_Z

  # Make sure data_YZ$Z (Z from selected sample), sumry_Z$var_Z (covariates of Z for non-selected sample) are matrices (in case they are scalars)
#  if (!is.matrix(data_YZ$Z)) data_YZ$Z <- as.matrix(data_YZ$Z)
  if (!is.matrix(sumry_Z$var_Z)) sumry_Z$var_Z <- as.matrix(sumry_Z$var_Z)  
  
  
  var_YZ_s = sumry_YZ$var_YZ
  mean_YZ_s = sumry_YZ$mean_YZ
  var_Z_s = var_YZ_s[-1,-1]
  mean_Z_s = mean_YZ_s[-1]
  zparams = length(mean_YZ_s)-1
  
  mean_Z_ns = sumry_Z$mean_Z
  var_Z_ns = sumry_Z$var_Z

  
  D.s <- solve(var_Z_s*(n_s-1))
  C.s <- drop(crossprod(mean_Z_s, D.s)) 
  F.s <- drop(C.s %*% mean_Z_s)
  
  mult.mat.s <- matrix(nrow=1+zparams, ncol=1+zparams)
  mult.mat.s[1,1] <- F.s + 1/n_s
  mult.mat.s[1,-1] <- mult.mat.s[-1,1] <- -C.s
  mult.mat.s[-1,-1] <- D.s
  ## Or from Brady
  # mult.mat.s <- solve(rbind(c(n_s,n_s*mean_YZ_s[-1]),cbind(n_s*mean_YZ_s[-1],(var_YZ_s[-1,-1] * (n_s-1) + n_s * tcrossprod(mean_YZ_s[-1])))))
  
  # MLEs of regression coefficients for Y|Z,S=1
  beta_YZ.Z_s <- drop(var_YZ_s[1,-1]%*%solve(var_YZ_s[-1,-1]))               # slopes
  beta_YZ.Z_s <- c(mean_YZ_s[1] - beta_YZ.Z_s%*%mean_YZ_s[-1], beta_YZ.Z_s)  # intercept
  # Residual variance
  var_Y.Z_s <- ((n_s-1)/(n_s-(zparams+1)))*drop(var_YZ_s[1,1]-var_YZ_s[1,-1]%*%solve(var_YZ_s[-1,-1])%*%var_YZ_s[-1,1])
  
  draws <- matrix(nrow = ndraws,ncol = 10)
  
  DRAWS_var_Y.Z_s <- (n_s-(zparams+1)) * var_Y.Z_s / rchisq(ndraws, n_s-(zparams+1))
  DRAWS_beta_YZ.Z_s <- t(apply(matrix(DRAWS_var_Y.Z_s),1, function(s) rmnorm(mean=beta_YZ.Z_s, varcov=s*mult.mat.s)))
  
  
  #Compute predicted values of X for all cases in subpopulation
  DRAWS_var_X_s <- rowSums(tcrossprod(DRAWS_beta_YZ.Z_s[,-1],var_YZ_s[-1,-1])*DRAWS_beta_YZ.Z_s[,-1])
  DRAWS_var_X_ns <- double(ndraws)
  DRAWS_var_XY_s <- matrix(0,nrow = ndraws, ncol = 3)
  DRAWS_S_X_ns <- double(ndraws)
  DRAWS_mean_X_s <- double(ndraws)
  DRAWS_mean_XY_s <- matrix(0,nrow = ndraws, ncol = 2)
  DRAWS_mean_X_ns <- double(ndraws)
  
  DRAWS_SMUB0_muy_ns <- double(ndraws)
  DRAWS_SMUB0_sigmayy_ns <- double(ndraws)
  DRAWS_SMUB0_cov_XY_ns <- double(ndraws)
  DRAWS_SMUB0_muy <- double(ndraws)
  
  for(d in 1:ndraws)
  {
    pd <- 0
    while(pd == 0)
    {
      #var_XY_invertible <- 0 #check for select cov XY
      #while(var_XY_invertible==0)
      #{
        var_XY_s_d <- matrix(0,2,2)
        var_XY_s_d[1,1] <- DRAWS_var_X_s[d]
        var_XY_s_d[1,2] <- var_XY_s_d[2,1] <- DRAWS_beta_YZ.Z_s[d,-1]%*%var_YZ_s[1,-1]
        var_XY_s_d[2,2] <- var_YZ_s[1,1]
        
        # if(var_XY_s_d[1,1] * var_XY_s_d[2,2] - var_XY_s_d[1,2]^2 > 0){
        #   var_XY_invertible <-1 # not 0, pd check pass! GOOD! 
        # }else{
        #   DRAWS_beta_YZ.Z_s[d,] <- rmnorm(mean=beta_YZ.Z_s, varcov=DRAWS_var_Y.Z_s[d]*mult.mat.s) #redraw beta
        #   DRAWS_var_X_s[d] <- sum(tcrossprod(DRAWS_beta_YZ.Z_s[d,-1],var_YZ_s[-1,-1])*DRAWS_beta_YZ.Z_s[d,-1]) #recalculate varx in cov(x,y) matrix
        # }
      #}
        
      # Calculate S_xxns^(0)(d)
      S_X_ns <- sum(tcrossprod(DRAWS_beta_YZ.Z_s[d,-1],var_Z_ns)*DRAWS_beta_YZ.Z_s[d,-1])
      # Draw cov(x,y) matrix for selected
      DRAW_var_XY_s<- riwish(n_s-1,var_XY_s_d)*(n_s-1)
      # Draw var_X_ns
      DRAW_var_X_ns <- (n_ns -1) * S_X_ns/rchisq(1,n_ns-1)
      # Calculate mean_X_s = beta_y0.z^(d)+beta_yz.z^(d)*mean_Z_s
      mean_X_s <- DRAWS_beta_YZ.Z_s[d,]%*%c(1,mean_YZ_s[-1])
      # Draw mean_XY_s from joint N(,)
      DRAW_mean_XY_s <- rmnorm(mean = c(mean_X_s, mean_YZ_s[1]), varcov = DRAW_var_XY_s/n_s)
      # Draw mean_X_ns from N(mean_X_ns^(d),/)
      DRAW_mean_X_ns <- rnorm(1,DRAWS_beta_YZ.Z_s[d,]%*%c(1,mean_Z_ns), sqrt(DRAW_var_X_ns/n_ns))
      
      # Draw phi from prior specified by phi_character
      if (!is.null(phi_character))
      {
        #Draw phi using provided string
        phi = eval(parse(text=phi_character));
        phi = pmax(.Machine$double.eps, phi);
      }
      
      # Draw
      DRAW_rho_XY_s <- DRAW_var_XY_s[1,2]/sqrt(DRAW_var_XY_s[1,1]*DRAW_var_XY_s[2,2])
      # g(phi) multiplier
      g_d <- ((phi + (1-phi)*DRAW_rho_XY_s) / ((1-phi) + phi*DRAW_rho_XY_s))
      # ratio of variance components for Y|Z, X|Z
      vratio_d <- DRAW_var_XY_s[2,2]/DRAW_var_XY_s[1,1]
      
      # Regression coefficients for non-selected cases, Y|Z,S=0
      DRAW_mean_Y_ns <- DRAW_mean_XY_s[2] + g_d * sqrt(vratio_d) *(DRAW_mean_X_ns - DRAW_mean_XY_s[1])
      DRAW_var_Y_ns <- DRAW_var_XY_s[2,2] + g_d^2 * vratio_d * (DRAW_var_X_ns - DRAW_var_XY_s[1,1])
      DRAW_cov_XY_ns <- DRAW_var_XY_s[1,2] +  g_d * sqrt(vratio_d) * (DRAW_var_X_ns - DRAW_var_XY_s[1,1])
      
      if(DRAW_var_X_ns * DRAW_var_Y_ns - DRAW_cov_XY_ns^2 > 0){ #check pd for non-select cov XY
        pd <- 1 #pd check pass! GOOD!!
      }else{
        DRAWS_beta_YZ.Z_s[d,] <- rmnorm(mean=beta_YZ.Z_s, varcov=DRAWS_var_Y.Z_s[d]*mult.mat.s)
        DRAWS_var_X_s[d] <- sum(tcrossprod(DRAWS_beta_YZ.Z_s[d,-1],var_YZ_s[-1,-1])*DRAWS_beta_YZ.Z_s[d,-1])
      }
      
    }
    
    
    # save into the DRAWS matrix
    DRAWS_S_X_ns[d] <- S_X_ns
    DRAWS_var_XY_s[d,] <- DRAW_var_XY_s[c(1,2,4)]
    DRAWS_var_X_ns[d] <- DRAW_var_X_ns
    DRAWS_mean_X_s[d] <- mean_X_s
    DRAWS_mean_XY_s[d,] <- DRAW_mean_XY_s
    DRAWS_mean_X_ns[d] <- DRAW_mean_X_ns
    
    DRAW_Y_all <- (n_s * DRAW_mean_XY_s[2] + n_ns * DRAW_mean_Y_ns) / (n_ns + n_s)
    DRAW_SMUB <- (mean_YZ_s[1] - DRAW_Y_all)/sqrt(DRAW_var_XY_s[2,2])
    
    # Calculation of SMUB(0)
    DRAW_SMUB0_scale <- DRAW_rho_XY_s * sqrt(vratio_d)
    DRAW_SMUB0_mean_Y_ns <- DRAW_mean_XY_s[2] + DRAW_SMUB0_scale * (DRAW_mean_X_ns - DRAW_mean_XY_s[1])
    DRAW_SMUB0_var_Y_ns <- DRAW_var_XY_s[2,2] + DRAW_SMUB0_scale^2 * (DRAW_var_X_ns - DRAW_var_XY_s[1,1])
    DRAW_SMUB0_cov_XY_ns <- DRAW_var_XY_s[1,2] +  DRAW_SMUB0_scale * (DRAW_var_X_ns - DRAW_var_XY_s[1,1])
    DRAW_SMUB0_Y_all <- (n_s * DRAW_mean_XY_s[2] + n_ns * DRAW_SMUB0_mean_Y_ns) / (n_ns + n_s)
    DRAW_SMUB0 <- (mean_YZ_s[1] - DRAW_SMUB0_Y_all)/sqrt(DRAW_var_XY_s[2,2])

    DRAWS_SMUB0_muy_ns[d] <- DRAW_SMUB0_mean_Y_ns
    DRAWS_SMUB0_sigmayy_ns[d] <- DRAW_SMUB0_var_Y_ns
    DRAWS_SMUB0_cov_XY_ns[d] <- DRAW_SMUB0_cov_XY_ns
    DRAWS_SMUB0_muy[d] <- DRAW_SMUB0_Y_all
    
    drawsPPM <- list(phi = phi,
                     muy_s = DRAW_mean_XY_s[2],
                     muy_ns = DRAW_mean_Y_ns,
                     sigmayy_s = DRAW_var_XY_s[2,2],
                     sigmayy_ns = DRAW_var_Y_ns,
                     rho_xy_s = DRAW_rho_XY_s,
                     muy = DRAW_Y_all,
                     smub = DRAW_SMUB,
                     smub0 = DRAW_SMUB0,
                     smab = DRAW_SMUB - DRAW_SMUB0)
    draws[d,] <- unlist(drawsPPM)
    
  }
  draws <- as.data.frame(draws)
  names(draws) <- names(unlist(drawsPPM))
  
  return(draws)
}

########### MLE ############
means_mle <- function(data_YZ, # selected/non-missing sample
                      data_Z, # non-selected/missing sample
                      intervals_at = c(0, 0.5, 1), # specific phi value
                      sig_digits_return = 5){
  require(mnormt);require(MCMCpack);require(ggplot2);require(nlme);require(magrittr);
  
  sumry_list = ErrorCheck(data_YZ, data_Z)
  sumry_YZ = sumry_list$sumry_YZ
  sumry_Z = sumry_list$sumry_Z

  #Sufficient statistics: non-selected/missing population
  mean_Z_no_selected = sumry_Z$mean_Z;
  var_Z_no_selected = sumry_Z$var_Z;
  n0 = sumry_Z$n_Z;
  
  #Sufficient statistics: selected/non-missing population
  mean_YZ_selected = sumry_YZ$mean_YZ;
  var_YZ_selected = sumry_YZ$var_YZ;
  n1 = sumry_YZ$n_YZ;
  #First step calculates the slopes
  beta_YZ.Z_selected = drop(var_YZ_selected[1,-1]%*%solve(var_YZ_selected[-1,-1]));
  #Second step calculates the intercept
  beta_YZ.Z_selected = c(mean_YZ_selected[1] - beta_YZ.Z_selected%*%mean_YZ_selected[-1],beta_YZ.Z_selected);
  npreds_as_dummy = length(mean_YZ_selected);
  var_Y.Z_selected = (n1 / (n1-npreds_as_dummy+1)) * drop(var_YZ_selected[1,1] - var_YZ_selected[1,-1]%*%solve(var_YZ_selected[-1,-1])%*%var_YZ_selected[-1,1])
  ZtZinv_selected = solve(rbind(c(n1,n1*mean_YZ_selected[-1]),cbind(n1*mean_YZ_selected[-1],(var_YZ_selected[-1,-1] * (n1-1) + n1 * tcrossprod(mean_YZ_selected[-1])))))
  
  ##NonBayes-calculation: point estimates of SMUB
  mean_X_selected = drop(beta_YZ.Z_selected%*%c(1,mean_YZ_selected[-1]));
  mean_X_no_selected = drop(beta_YZ.Z_selected%*%c(1,mean_Z_no_selected));
  var_X_selected = drop(tcrossprod(beta_YZ.Z_selected[-1],var_YZ_selected[-1,-1])%*%beta_YZ.Z_selected[-1]);
  mean_X_pop = (mean_X_selected*n1 + mean_X_no_selected*n0)/(n0+n1);
  cor_XY_selected = drop(beta_YZ.Z_selected[-1]%*%var_YZ_selected[1,-1])/sqrt(var_X_selected * var_YZ_selected[1,1]);
  smub_point_est = (intervals_at + (1 - intervals_at) * cor_XY_selected) / (intervals_at * cor_XY_selected + (1 - intervals_at)) * (mean_X_selected - mean_X_pop) / sqrt(var_X_selected);
  names(smub_point_est) = intervals_at;
  
  ##NonBayes-calculation: point estimates of SMAB
  smab_point_est = (intervals_at * (1 - cor_XY_selected^2)) / (intervals_at * cor_XY_selected + (1 - intervals_at)) * (mean_X_selected - mean_X_pop) / sqrt(var_X_selected);
  names(smab_point_est) = intervals_at;
  
  ## point estimates of MUB
  mub_point_est = smub_point_est*sqrt(var_YZ_selected[1,1])
  names(mub_point_est) = intervals_at;
  
  ## population Y
  mean_Y_pop = mean_YZ_selected[1] - mub_point_est
  names(mean_Y_pop) = intervals_at;
  
  ## final result
  sumry_table = data.frame(cbind(phi = intervals_at,
                                 mub_point_est,
                                 smub_point_est,
                                 smab_point_est,
                                 mean_Y_pop))
  
  ## return final result
  return(list(sumry_table = round(sumry_table,sig_digits_return),
              proxy_strength = cor_XY_selected))
}