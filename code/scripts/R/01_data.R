# 01_data.R
# Description: Functions to simulate datasets and generate missing values
# Author: Pierre CATOIRE
# Date: 2025-12-02
# Purpose:
#   - Generate train/test sets with continuous variables,
#   - induce missing values,
#   - approximate the oracle and pragmatic targets

# --- Preamble ---
library(tidyverse)
library(parallel)

# --- Functions ---

#' Set the intercept of phi
#'
#' This function sets the intercept of phi such as the expectancy of M1
#' is approximately missCoef. This allows increasing the proportion of
#' missingness without affecting the other parameters of the generation model.
#'
#' @param missCoef missingness proportion target
#' @param phi Numeric vector of logistic regression coefficients for M1.
#' @param theta List of parameters for X1, X2, and Y.
#'
#' @return the value of intercept for the logistic regression of phi such that
#' E[M1=1] is approximately missCoef.
TuneInterceptM1 = function(missCoef, theta, phi,
                           ns = 1e6, tol = 1e-4, maxIter = 200, lim = c(-12,4)){
  get_pM1_given_intercept = function(intercept, model, theta, phi, ns = 1e6){
    phi[1] = intercept
    X1 = rnorm(ns,theta[["X1"]][["mu"]],theta[["X1"]][["sigma"]])
    X2 = rnorm(ns,theta[["X1"]][["mu"]],theta[["X1"]][["sigma"]])
    Y = cbind(1,X1,X2) %*% theta[["Y"]][["beta"]] + rnorm(ns,0,theta[["Y"]][["sigma"]])
    pM1 = 1/(1+exp(-(cbind(1,X1,X2,Y) %*% phi)))
    return(mean(pM1))
  }
  min = lim[1]
  max = lim[2]
  for(iter in 1:maxIter){
    candidate = mean(c(min,max))
    delta = get_pM1_given_intercept(candidate, model, theta, phi, ns = 1e6)-missCoef
    if(abs(delta) < tol){
      return(candidate)
    }else if(delta>=0){
      max = candidate
    }else if(delta < 0){
      min = candidate
    }
  }
  warning("Tuning phi intercept failed. Returning the last iteration.")
  return(candidate)
}

#' Compute E[Y|X2 = x2, M1 = 1]
#'
#' This function computes the expected value of Y given X2 and M1 = 1
#' using a Monte Carlo approximation.
#'
#' @param X2_vec Numeric vector of X2 values.
#' @param theta List of parameters for X1, X2, and Y.
#' @param phi Numeric vector of logistic regression coefficients for M1.
#' @param n_mc Number of Monte Carlo samples (default 1e5).
#'
#' @return Numeric vector of expected Y values, same length as X2_vec.
ComputeEYGivenX2M1_1 = function(X2_vec, theta, phi, n_mc = 1e5) {
  n_X2 = length(X2_vec)
  
  # 1. Sample X1 (same for all X2)
  X1_mc = rnorm(n_mc, mean = theta[["X1"]][["mu"]], sd = theta[["X1"]][["sigma"]])
  
  # 2. Expand X1 and X2 to matrices for vectorized computation
  X1_mat = matrix(X1_mc, nrow = n_mc, ncol = n_X2, byrow = FALSE)  # n_mc × n_X2
  X2_mat = matrix(X2_vec, nrow = n_mc, ncol = n_X2, byrow = TRUE)   # n_mc × n_X2
  
  # 3. Compute mean of Y for each (X1, X2)
  mu_Y_mat = theta[["Y"]][["beta"]][1] + 
    theta[["Y"]][["beta"]][2] * X1_mat + 
    theta[["Y"]][["beta"]][3] * X2_mat
  
  # 4. Sample Y
  Y_mat = mu_Y_mat + matrix(rnorm(n_mc * n_X2, mean = 0,
                                  sd = theta[["Y"]][["sigma"]]), 
                            nrow = n_mc, ncol = n_X2)
  
  # 5. Compute probability M1=1
  eta_mat = phi[1] + phi[2]*X1_mat + phi[3]*X2_mat + phi[4]*Y_mat
  p_M1_mat = 1 / (1 + exp(-eta_mat))
  
  # 6. Weighted expectation column-wise
  EY_vec = colSums(Y_mat * p_M1_mat) / colSums(p_M1_mat)
  
  return(EY_vec)
}

#' Compute E[Y|X1 = x1, X2 = x2, M1 = 1]
#'
#' This function computes the expected value of Y given X1, X2 and M1
#' using a Monte Carlo approximation.
#'
#' @param X1_vec Numeric vector of X1 values.
#' @param X2_vec Numeric vector of X2 values.
#' @param M1_vec Numeric vector of M1 values.
#' @param theta List of parameters for X1, X2, and Y.
#' @param phi Numeric vector of logistic regression coefficients for M1.
#' @param n_mc Number of Monte Carlo samples (default 1e5).
#'
#' @return Numeric vector of expected Y values, same length as X1_vec.
ComputeEYGivenX1X2M1 = function(X1_vec, X2_vec, M1_vec, theta, phi,
                                n_mc = 1e5) {
  stopifnot(length(X1_vec) == length(X2_vec), length(X1_vec) == length(M1_vec))
  n_points = length(X1_vec)
  
  # 1. Sample Y for each (X1, X2) point
  X1_mat = matrix(rep(X1_vec, each = n_mc), nrow = n_mc, ncol = n_points)
  X2_mat = matrix(rep(X2_vec, each = n_mc), nrow = n_mc, ncol = n_points)
  
  # 2. Compute mean of Y
  mu_Y_mat = theta[["Y"]][["beta"]][1] +
    theta[["Y"]][["beta"]][2] * X1_mat +
    theta[["Y"]][["beta"]][3] * X2_mat
  
  # 3. Sample Y
  Y_mat = mu_Y_mat + matrix(rnorm(n_mc * n_points, mean = 0,
                                  sd = theta[["Y"]][["sigma"]]),
                            nrow = n_mc, ncol = n_points)
  
  # 4. Compute probability M1=1
  eta_mat = phi[1] + phi[2]*X1_mat + phi[3]*X2_mat + phi[4]*Y_mat
  p_M1_1_mat = 1 / (1 + exp(-eta_mat))
  
  # 5. Set weights according to observed M1 for each individual
  #    If M1 = 1, weight = P(M1=1); if M1=0, weight = P(M1=0)
  weights_mat = matrix(NA, nrow = n_mc, ncol = n_points)
  weights_mat[, M1_vec == 1] = p_M1_1_mat[, M1_vec == 1]
  weights_mat[, M1_vec == 0] = 1 - p_M1_1_mat[, M1_vec == 0]
  
  # 6. Compute weighted expectation column-wise
  EY_vec = colSums(Y_mat * weights_mat) / colSums(weights_mat)
  
  return(EY_vec)
}

#' Generate a sample of observations
#'
#' This function generates a 
#'
#' @param X1_vec Numeric vector of X1 values.
#' @param X2_vec Numeric vector of X2 values.
#' @param M1_vec Numeric vector of M1 values.
#' @param theta List of parameters for X1, X2, and Y.
#' @param phi Numeric vector of logistic regression coefficients for M1.
#' @param n_mc Number of Monte Carlo samples (default 1e5).
#'
#' @return Numeric vector of expected Y values, same length as X1_vec.
#' @export
SimulateDataContinuous = function(n, theta, phi, missCoef = .5, n_mc = 1e4, type = "train"){
  if(!(type %in% c("train","test"))) stop("type must be either \"train\" or \"test\"")
  
  # --- set phi[1] to match the expected proportion of missingness missCoef ---
  # NB: if missCoef == 0, phi does not matter
  if(missCoef > 0){phi[1] = TuneInterceptM1(missCoef, theta, phi)}
  
  # --- Generate continuous predictors ---
  X1 = rnorm(n, mean = theta[["X1"]][["mu"]], sd = theta[["X1"]][["sigma"]])
  X2 = rnorm(n, mean = theta[["X2"]][["mu"]], sd = theta[["X2"]][["sigma"]])
  
  # --- Generate outcome Y ---
  X_mat = cbind(1,X1,X2)
  # Note: ORACLE_MU is the expected value of Y given X1,x2
  ORACLE_MU = X_mat %*% theta[["Y"]][["beta"]]
  Y = ORACLE_MU +rnorm(n,0,theta[["Y"]][["sigma"]])
  
  # --- Generate missingness indicator M1 ---
  XY_mat = cbind(1,X1,X2,Y)
  probM1 = (1+exp(-XY_mat %*% phi))^-1
  if(missCoef > 0){
    M1 = rbinom(length(probM1), size = 1, prob = probM1)
  } else{
    M1 = rep(0,n)
  }
  
  # --- Generate observed X1 with missing values ---
  X1OBS = ifelse(M1 == 0, X1, NA)
  
  if(type == "train"){
    return(data.frame(X1,X2,Y,M1,X1OBS))
  }else if(type == "test"){
    # --- Generate the expected values of Y ---
    # Note: these quantities are used for estimate missingness-conditioned (MC)
    # and Missingness-Unconditioned (MU) probabilities
    
    # E[Y|X2]
    Y_GIVEN_X2 = theta[["Y"]][["beta"]][1] + 
      theta[["Y"]][["beta"]][2] * theta[["X1"]][["mu"]] + 
      theta[["Y"]][["beta"]][3] * X2
    
    # E[Y|X2,M1=1]
    Y_GIVEN_X2M1 = ComputeEYGivenX2M1_1_parallel(X2,theta,phi,n_mc)
    
    # --- Generate the pragmatic and oracle MU and MC probabilities ---
    
    # Oracle MU: OMU = E[Y|X1,X2]
    # Oracle MC: OMC = E[Y|X1,X2,M]
    ORACLE_MC = ComputeEYGivenX1X2M1_parallel(X1,X2,M1,theta,phi,n_mc)
    
    # Pragmatic MU: E[Y|X1,X2] if M1 = 0, E[Y|X2] if M1 = 1
    PRAGMATIC_MU = ifelse(M1 == 0, ORACLE_MU, Y_GIVEN_X2)
    
    # PRAGMATIC_MC: E[Y|X1,X2,M1=0] if M1 = 0, E[Y|X2,M1=1] if M1 = 1
    PRAGMATIC_MC = ifelse(M1 == 0, ORACLE_MC, Y_GIVEN_X2M1)
    
    return(data.frame(X1,X2,Y,M1,X1OBS,
                      ORACLE_MU,ORACLE_MC,PRAGMATIC_MU,PRAGMATIC_MC))
  }
}

#' Merge the simulated datasets
#'
#' @param datasetList a list of simulated datasets
#' @param missCoefs the vectors of coefficients for missingness proportion
#'
#' @return a data.frame of merged simulated datasets
MergeDatasets = function(datasetList,
                         missCoefs){
  result = list()
  i = 1
  for(model in names(datasetList)){
    for(i in 1:length(datasetList[[model]][["TRAIN"]])){
      dTrain = datasetList[[model]][["TRAIN"]][[i]]
      dTrain[["SET"]] = "TRAIN"
      dTest = datasetList[[model]][["TEST"]][[i]]
      dTest[["SET"]] = "TEST"
      dat = dplyr::bind_rows(dTrain,dTest)
      dat[["MODEL"]] = model
      dat[["MISSCOEF"]] = missCoef
      result[[i]] = dat
    }
  }
  return(do.call(rbind, result))
}

#' Compute E[Y | X1, X2, M1]
#'
#' This function computes the expected value of Y given X1, X2 and M1, which
#' is used for computation of the MC probability reference.
#'
#' @param X1_vec the vector of realised values of X1
#' @param X2_vec the vector of realised values of X2
#' @param M1_vec the vector of realised values of M1
#' @param theta the vector of parameters of Pr(Y, X1, X2)
#' @param phi the vector of parameters of Pr(Y | X1, X2, M1)
#' @param n_mc the number of samples of Monte Carlo simulation
#' @param n_cores the number of processors (for parallelisation)
#' 
#' @return the vector of E[Y | X1, X2, M1]
ComputeEYGivenX1X2M1_parallel = function(X1_vec, X2_vec, M1_vec, theta, phi, n_mc = 1e5, n_cores = detectCores() - 1) {
  stopifnot(length(X1_vec) == length(X2_vec), length(X1_vec) == length(M1_vec))
  n_points = length(X1_vec)
  
  cl = makeCluster(n_cores)
  clusterExport(cl, varlist = c("X1_vec", "X2_vec", "M1_vec", "theta", "phi", "n_mc", "n_points"), envir = environment())
  
  # Split n_mc into chunks
  chunks = split(1:n_mc, cut(1:n_mc, n_cores, labels = FALSE))
  
  results = parLapply(cl, chunks, function(idx) {
    X1_mat = matrix(rep(X1_vec, each = length(idx)), nrow = length(idx), ncol = n_points)
    X2_mat = matrix(rep(X2_vec, each = length(idx)), nrow = length(idx), ncol = n_points)
    
    mu_Y_mat = theta[["Y"]][["beta"]][1] + 
      theta[["Y"]][["beta"]][2] * X1_mat + 
      theta[["Y"]][["beta"]][3] * X2_mat
    Y_mat = mu_Y_mat + matrix(rnorm(length(idx) * n_points, 0, theta[["Y"]][["sigma"]]),
                              nrow = length(idx), ncol = n_points)
    
    eta_mat = phi[1] + phi[2]*X1_mat + phi[3]*X2_mat + phi[4]*Y_mat
    p_M1_1_mat = 1 / (1 + exp(-eta_mat))
    
    # Compute weights
    weights_mat = matrix(NA, nrow = length(idx), ncol = n_points)
    weights_mat[, M1_vec == 1] = p_M1_1_mat[, M1_vec == 1]
    weights_mat[, M1_vec == 0] = 1 - p_M1_1_mat[, M1_vec == 0]
    
    list(Y_mat = Y_mat, weights_mat = weights_mat)
  })
  
  stopCluster(cl)
  
  # Combine results
  Y_all = do.call(rbind, lapply(results, `[[`, "Y_mat"))
  W_all = do.call(rbind, lapply(results, `[[`, "weights_mat"))
  
  EY_vec = colSums(Y_all * W_all) / colSums(W_all)
  return(EY_vec)
}

#' Compute E[Y | X2, M1 = 1]
#'
#' This function computes the expected value of Y given X1, X2 and M1 = 1, which
#' is used for computation of the MC probability reference.
#'
#' @param X2_vec the vector of realised values of X2
#' @param theta the vector of parameters of Pr(Y, X1, X2)
#' @param phi the vector of parameters of Pr(Y | X1, X2, M1)
#' @param n_mc the number of samples of Monte Carlo simulation
#' @param n_cores the number of processors (for parallelisation)
#' 
#' @return the vector of E[Y | X2, M1 = 1]
ComputeEYGivenX2M1_1_parallel = function(X2_vec, theta, phi, n_mc = 1e5, n_cores = detectCores() - 1) {
  n_X2 = length(X2_vec)
  cl = makeCluster(n_cores)
  
  # Split n_mc into roughly equal chunks
  chunks = split(1:n_mc, cut(1:n_mc, n_cores, labels = FALSE))
  
  # Export variables to cluster
  clusterExport(cl, varlist = c("X2_vec", "theta", "phi"), envir = environment())
  
  results = parLapply(cl, chunks, function(idx) {
    X1_mc = rnorm(length(idx), mean = theta[["X1"]][["mu"]], sd = theta[["X1"]][["sigma"]])
    X1_mat = matrix(X1_mc, nrow = length(idx), ncol = n_X2, byrow = FALSE)
    X2_mat = matrix(X2_vec, nrow = length(idx), ncol = n_X2, byrow = TRUE)
    
    mu_Y_mat = theta[["Y"]][["beta"]][1] + 
      theta[["Y"]][["beta"]][2]*X1_mat + 
      theta[["Y"]][["beta"]][3]*X2_mat
    Y_mat = mu_Y_mat + matrix(rnorm(length(idx)*n_X2, 0, theta[["Y"]][["sigma"]]),
                              nrow = length(idx), ncol = n_X2)
    eta_mat = phi[1] + phi[2]*X1_mat + phi[3]*X2_mat + phi[4]*Y_mat
    p_M1_mat = 1 / (1 + exp(-eta_mat))
    
    return(list(Y_mat = Y_mat, p_M1_mat = p_M1_mat))
  })
  
  stopCluster(cl)
  
  # Combine results
  Y_all = do.call(rbind, lapply(results, `[[`, "Y_mat"))
  p_all = do.call(rbind, lapply(results, `[[`, "p_M1_mat"))
  
  EY_vec = colSums(Y_all * p_all) / colSums(p_all)
  return(EY_vec)
}