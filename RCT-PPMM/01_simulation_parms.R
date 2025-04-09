# Numerically solve for B0 for various combos of n_RCT/n (sampling fraction) and B_Z, B_W, B_U
library(tidyverse)

# function to find B0 [INTERCEPT in SELECTION MODEL]
findB0 <- function(SFRAC, CORR_U_W, B_Z, B_W, B_U)
{
  set.seed(7711)
  
  # covariates
  N <- 10000000
  z <- rnorm(N)
  u <- rnorm(N)
  w <- CORR_U_W*u + rnorm(N, mean = 0, sd = sqrt(1 - CORR_U_W^2))
  
  # Find value of intercept for specified selection fraction
  if (SFRAC==0.5) {
    chosen.b0 <- 0
  } else {
    f <- function(b0) abs(SFRAC - mean(1 - 1/(1+exp(b0 + B_Z*z + B_W*w + B_U*u))))
    chosen.b0 <- optimize(f, c(-8, 2))$minimum
  }
  return(tibble(SFRAC=SFRAC, CORR_U_W=CORR_U_W, B_Z=B_Z, B_W=B_W, B_U=B_U, B0=chosen.b0))
}

# parameters of SELECTION model -----
sel <- vector()
i <- 0
betas <- tibble(B_Z=c(0,1,1,1,1), 
                B_W=c(0,0,1,0,1), 
                B_U=c(1,1,1,-1,-1))
time.start <- Sys.time()
for (CORR_U_W in c(0, 0.3))
{
  i <- i + 1
  print(paste("i=",i,sep=""))
  for (j in 1:nrow(betas))
  {
    print(j)
    parms <- findB0(SFRAC=0.04, CORR_U_W, betas$B_Z[j], betas$B_W[j], betas$B_U[j])
    sel <- rbind(sel, parms)
  }
}
time.end <- Sys.time()
time.end - time.start
sel$B0 <- round(sel$B0,2)
table(sel$B0, sel$CORR_U_W)
rm(i, j, time.start, time.end, parms, betas)
sel <- sel %>%
  mutate(selection = case_when(B_Z == 0 & B_W == 0 & B_U == 1 ~ "SEL{U}",
                               B_Z == 1 & B_W == 0 & B_U == 1 ~ "SEL{Z,U}",
                               B_Z == 1 & B_W == 1 & B_U == 1 ~ "SEL{Z,W,U}",
                               B_Z == 1 & B_W == 0 & B_U == -1 ~ "SEL{Z,-U}",
                               B_Z == 1 & B_W == 1 & B_U == -1 ~ "SEL{Z,W,-U}"),
         sel_neg_u = (B_U < 0))

# parameters of EFFECT MODIFICATION model -----
# y_0 <- 0 + z + w + u + rnorm(n, 0, sqrt(S2))
# y_1 <- 1 + (1+C_Z)*z + (1+C_W)*w + (1+C_U)*u + rnorm(n, 0, sqrt(S2))
# y_1 - y_0 = 1 + C_Z*z + C_W*w + C_U*u + e_1-e_0
em <- tibble(C_Z=c(0,1,1), C_W=c(0,0,1), C_U=1) %>%
  mutate(effectmod = case_when(C_Z == 0 & C_W == 0 & C_U == 1 ~ "EM{U}",
                               C_Z == 1 & C_W == 0 & C_U == 1 ~ "EM{Z,U}",
                               C_Z == 1 & C_W == 1 & C_U == 1 ~ "EM{Z,W,U}"))

# combine, adding residual error -----
parms <- expand_grid(sel, em)
parms <- expand_grid(S2=c(1, 4, 13), parms)

# Set number
parms <- parms %>%
  mutate(simset = 1:n())

# Save
save(parms, file="./simdata/simParms.RData")
