############### stop function ################
stopfn<-function(message) {stop(message, call.=FALSE)}

############### Turn Microdata to Summary Statistics ###############
MicroToSummary <- function(Z, # microdata of Z
                           Y = NULL){ 
  if (is.null(Y)){
    # Y is NULL
    dt <- Z
  }else{
    # combine Z and Y if Y has value
    dt <- cbind(Y, Z)
  }
  # Mean for (Y,Z)
  dt_mean <- colMeans(dt)
  # Covariance matrix for (Y,Z)
  dt_var <- var(dt)
  # Sample size
  n <- nrow(dt)
  # Summary statistics for non-selected/missing
  if(is.null(Y)){
    sumry <- list(mean_Z = dt_mean,
                  var_Z = dt_var,
                  n_Z = n)
  }else{# Summary statistics for selected/non-missing
    sumry <- list(mean_YZ = dt_mean,
                  var_YZ = dt_var,
                  n_YZ = n)
  }
  return(sumry)
}

################## Error Check #############################
ErrorCheck <-function(data_YZ, # selected/non-missing sample
                      data_Z,  # non-selected/missing sample
                      prop_check = FALSE, # Y default not binary
                      mi_check = FALSE){ # default not multiple imputation
  
  ## test if list of Y or Z missing in data_YZ
  if(all(c("Y","Z") %in% names(data_YZ))){# microdata
    # check if data_YZ has missing (need Y no missing & Z no missing)
    if (any(is.na(data_YZ$Y)) | any(is.na(data_YZ$Z))){# if Y missing or Z missing
      stopfn("Requires complete Microdata data_YZ (list of Y and Z)")
    }
    if(!is.matrix(data_YZ$Z)){# check if Z is a matrix, if not, turn Z into matrix
      data_YZ$Z <- as.matrix(data_YZ$Z)
      colnames(data_YZ$Z) = c("Z")
    }
    if(is.matrix(data_YZ$Y)){ # check if Y is a matrix
      if(ncol(data_YZ$Y) != 1){ # check if the column of the matrix is not 1
        stopfn("Microdata data_YZ need to have only one column of Y (both vector and matrix work)")
      }else{
        data_YZ$Y = c(data_YZ$Y)
      }
    }
    # check if Y is not a number -> change to 0/1 indicator
    if(!is.numeric(data_YZ$Y)){
      if(length(unique(data_YZ$Y))!=2){
        stopfn("Requires binary Y")
      }else{
        # change factor Y -> binary Y
        Binary_Y = fastDummies::dummy_cols(data.frame(Y = data_YZ$Y),
                                           "Y",
                                           remove_first_dummy = T,
                                           remove_selected_columns = T)
        # change Y to vector Y
        data_YZ$Y = Binary_Y[,1]
        cat("Input Y in data_YZ is string -> make Y to be indicator of",
            str_sub(names(Binary_Y),3),"\n");
      }
    }
    # if prop -> binary
    if(prop_check){
      if(!identical(c(0,1), as.double(sort(unique(data_YZ$Y))))){
        stopfn("Requires binary Y")
      }
    }
    
    #############################
    # Make sure data_YZ$Z (Z from selected sample) is matrix (in case it is a scalar)
    # if (!is.matrix(data_YZ$Z)){
    #   data_YZ$Z <- as.matrix(data_YZ$Z)
    #   colnames(data_YZ$Z) = c("Z")
    #}
    #############################
    
    # check if Y and Z have same dimension
    if (length(data_YZ$Y) != nrow(data_YZ$Z) ){
      stopfn("Microdata data_YZ need to have same observation numbers for Y and Z")
    }
    # turn Microdata data_YZ to Summary Statistics sumry_YZ
    sumry_YZ = MicroToSummary(data_YZ$Z,data_YZ$Y)
  }else if(all(c("mean_YZ", "var_YZ", "n_YZ") %in% names(data_YZ))){ # summary
    if(mi_check){ # check for multiple imputation
      stopfn("Requires Microdata data_YZ (list of Y and Z)")
    }
    if(prop_check){ # check for proportion
      stopfn("Requires Microdata data_YZ (list of Y and Z)")
    }
    # check if data_YZ has missing
    if (any(is.na(data_YZ$mean_YZ)) | any(is.na(data_YZ$var_YZ)) | is.na(data_YZ$n_YZ)){# if mean/var/n missing
      stopfn("Requires complete Summary Statistics data_YZ (list of mean_YZ, var_YZ and n_YZ)")
    }
    
    # check if length(mean) = ncol(covariance matrix) = nrows(covariance matrix)
    if ((length(data_YZ$mean_YZ) != nrow(data_YZ$var_YZ)) | (length(data_YZ$mean_YZ) != ncol(data_YZ$var_YZ))){
      stopfn("Summary Statistics data_YZ need to have same number of covariates for Y and Z")
    }
    # return Summary Statistics
    sumry_YZ = data_YZ 
  }else{
    stopfn("Requires data_YZ to be Microdata (list of Y and Z) or Summary Statistics(list of mean_YZ, var_YZ and n_YZ)")
  }
  
  ## test if list of Z missing in data_Z
  if("Z" %in% names(data_Z)){# microdata data_Z
    # check if data_Z has missing (need Z no missing)
    if (any(is.na(data_Z$Z))){# if Z missing
      stopfn("Requires complete Microdata data_Z (list of Z)")
    }
    if(!is.matrix(data_Z$Z)){# check if Z is a matrix, if not, turn Z into matrix
      data_Z$Z <- as.matrix(data_Z$Z)
      colnames(data_Z$Z) = c("Z")
    }
    # turn Microdata data_Z to Summary Statistics sumry_Z
    sumry_Z = MicroToSummary(data_Z$Z)
  }else if(all(c("mean_Z", "var_Z", "n_Z") %in% names(data_Z))){ # summary
    if(mi_check){ # check for multiple imputation
      stopfn("Requires Microdata data_Z (list of Z)")
    }
    # check if data_Z has missing
    if (any(is.na(data_Z$mean_Z)) | any(is.na(data_Z$var_Z)) | is.na(data_Z$n_Z)){# if mean/var/n missing
      stopfn("Requires complete Summary Statistics data_Z (list of mean_Z, var_Z and n_Z)")
    }
    if(!is.matrix(data_Z$var_Z)){# check if var_Z is a matrix, if not, turn var_Z into matrix
      data_Z$var_Z <- as.matrix(data_Z$var_Z)
      colnames(data_Z$var_Z) = c("Z")
      rownames(data_Z$var_Z) = c("Z")
    }
    
    # check if length(mean) = ncol(covariance matrix) = nrows(covariance matrix)
    if ((length(data_Z$mean_Z) != nrow(data_Z$var_Z)) | (length(data_Z$mean_Z) != ncol(data_Z$var_Z))){
      stopfn("Summary Statistics data_Z need to have same number of covariates for Z")
    }
    # return Summary Statistics
    sumry_Z = data_Z 
  }else{
    stopfn("Requires data_Z to be Microdata (list of Z) or Summary Statistics(list of mean_Z, var_Z and n_Z)")
  }
  
  ## check if Z match in data_YZ and data_Z
  # check if covariates Z matched in data_YZ and data_Z
  if(!identical((colnames(sumry_YZ$var_YZ)[-1]), colnames(sumry_Z$var_Z))){
    stopfn("data_YZ and data_Z need to have exactly same covariates for Z")
  }
  # return data_YZ, data_Z, sumry_YZ, sumry_Z
  return(list(data_YZ = data_YZ,
              data_Z = data_Z,
              sumry_YZ = sumry_YZ,
              sumry_Z = sumry_Z))
} 
