# 02_methods.R
# Description: Functions to train the models and perform the predictions
# Author: Pierre CATOIRE
# Date: 2025-12-02
# Purpose:
# - train each model on a training set
# - perform the prediction on a testing set

# Preamble
library(tidyverse)
library(mice)

# Functions
FitPatternSubmodels = function(dTrain){
  modelPS = list()
  if(any(dTrain[["M1"]] == 0)){
    modelPS[["M0"]] = lm(formula("Y ~ X1OBS + X2"),
                         data = dTrain[dTrain[["M1"]] == 0,
                                       c("Y","X1OBS","X2")])
  }
  if(any(dTrain[["M1"]] == 1)){
    modelPS[["M1"]] = lm(formula("Y ~ X2"),
                         data = dTrain[dTrain[["M1"]] == 1,
                                       c("Y","X2")])
  }
  return(modelPS)
}

PredictPatternSubmodels = function(modelPS,dTest){
  result = rep(NA,nrow(dTest))
  if(any(dTest[["M1"]] == 0)){
    if(is.null(modelPS[["M0"]])){
      warning("Model does not contain a submodel for M1 == 0, and some rows in dTest have M1 == 0.
            \nPrediction is NA for these cases.")
    }
    else{
      result[dTest[["M1"]] == 0] = predict(modelPS[["M0"]],
                                           newdata = dTest[dTest[["M1"]] == 0,])
    }
  }
  if(any(dTest[["M1"]] == 1)){
    if(is.null(modelPS[["M1"]])){
      warning("Model does not contain a submodel for M1 == 1, and some rows in dTest have M1 == 1.
            \nPrediction is NA for these cases.")
    }
    else{
      result[dTest[["M1"]] == 1] = predict(modelPS[["M1"]],
                                           newdata = dTest[dTest[["M1"]] == 1,])
    }
  }
  return(result)
}

fitMarg = function(dTrain,
                   max_iter = 1000,
                   tol = 1e-6,
                   verbose = FALSE) {
  # dTrain: numeric matrix or dTrain.frame (n x p) with NAs for missing entries
  data = as.matrix(dTrain[,c("X1OBS","X2","Y")])
  n = nrow(data)
  p = ncol(data)
  
  # init
  mu = colMeans(data, na.rm = TRUE)           # p
  Sigma = stats::cov(data, use = "pairwise.complete.obs")  # p x p
  # ensure positive definite (tiny ridge)
  Sigma = Sigma + diag(1e-8, p)
  
  iter = 0
  repeat {
    iter = iter + 1
    S1 = numeric(p)              # sum of expected Z_i
    S2 = matrix(0, p, p)         # sum of expected Z_i Z_i^T
    
    # E-step: for each observation, compute E[Z_i | observed] and E[Z_i Z_i^T | observed]
    for (i in 1:n) {
      xi = data[i, ]
      obs_idx = which(!is.na(xi))
      miss_idx = which(is.na(xi))
      
      if (length(miss_idx) == 0) {
        # fully observed
        zi = xi
        S1 = S1 + zi
        S2 = S2 + (zi %*% t(zi))
      } else {
        # partition parameters
        mu_o = mu[obs_idx]
        mu_m = mu[miss_idx]
        
        Sigma_oo = Sigma[obs_idx, obs_idx, drop = FALSE]
        Sigma_mo = Sigma[miss_idx, obs_idx, drop = FALSE]
        Sigma_mm = Sigma[miss_idx, miss_idx, drop = FALSE]
        
        x_o = xi[obs_idx]
        
        # stabilize inversion if needed
        eps = 1e-8
        Sigma_oo_inv = solve(Sigma_oo + diag(eps, nrow(Sigma_oo)))
        
        # conditional mean of missing given observed (E-step)
        mu_cond = as.numeric(mu_m + Sigma_mo %*% (Sigma_oo_inv %*% (x_o - mu_o)))
        
        # conditional covariance of missing (Var[missing | observed])
        Sigma_cond = Sigma_mm - Sigma_mo %*% (Sigma_oo_inv %*% t(Sigma_mo))
        
        # build expected full vector E[Z_i]
        zi_exp = numeric(p)
        zi_exp[obs_idx] = x_o
        zi_exp[miss_idx] = mu_cond
        
        # build expected outer product E[Z Z^T]
        EZZ = matrix(0, p, p)
        
        # obs-obs block: x_o x_o^T
        EZZ[obs_idx, obs_idx] = x_o %*% t(x_o)
        
        # miss-miss block: Sigma_cond + mu_cond mu_cond^T
        EZZ[miss_idx, miss_idx] = Sigma_cond + mu_cond %*% t(mu_cond)
        
        # miss-obs block: mu_cond x_o^T
        EZZ[miss_idx, obs_idx] = mu_cond %*% t(x_o)
        EZZ[obs_idx, miss_idx] = t(EZZ[miss_idx, obs_idx])
        
        S1 = S1 + zi_exp
        S2 = S2 + EZZ
      }
    } # end loop over i
    
    # M-step: update mu and Sigma
    mu_new = S1 / n
    Sigma_new = S2 / n - (mu_new %*% t(mu_new))
    
    # ensure symmetry and numeric stability
    Sigma_new = (Sigma_new + t(Sigma_new)) / 2
    Sigma_new = Sigma_new + diag(1e-10, p)
    
    # convergence check
    delta = max(abs(mu_new - mu), abs(Sigma_new - Sigma))
    mu = mu_new
    Sigma = Sigma_new
    
    if (verbose && iter %% 50 == 0) {
      cat("iter =", iter, "delta =", format(delta, digits = 6), "\n")
    }
    
    if (delta < tol || iter >= max_iter) break
  }
  
  list(mu = mu, Sigma = Sigma, iter = iter, converged = (delta < tol))
}

cond_mean_vec = function(fit, observed_vals_mat, observed_idx, target_idx) {
  
  mu = fit$mu
  Sigma = fit$Sigma
  
  mu_o = mu[observed_idx]
  mu_t = mu[target_idx]
  
  Sigma_to = Sigma[target_idx, observed_idx, drop = FALSE]   # 1 x k
  Sigma_oo = Sigma[observed_idx, observed_idx, drop = FALSE] # k x k
  
  # stabilised inverse
  Sigma_oo_inv = solve(Sigma_oo + diag(1e-8, nrow(Sigma_oo)))
  
  # This is the key fix:
  coef = Sigma_to %*% Sigma_oo_inv    # 1 x k
  
  # diff is (n × k)
  diff = sweep(observed_vals_mat, 2, mu_o, `-`)
  
  # vectorised conditional means
  out = as.vector(mu_t + diff %*% t(coef))
  return(out)
}

# convenience wrappers for your variables: assume columns are c("X1OBS","X2","Y")
predict_y_given_x1_x2 = function(fit, data) {
  x1 = data[["X1OBS"]]
  x2 = data[["X2"]]
  if (length(x1) != length(x2)){
    stop("x1 and x2 must have same length")
  }
  observed_vals_mat = cbind(x1, x2)   # n × 2
  cond_mean_vec(
    fit = fit,
    observed_vals_mat = observed_vals_mat,
    observed_idx = c(1, 2),
    target_idx = 3
  )
}


predict_y_given_x2 = function(fit, x2) {
  observed_vals_mat = cbind(x2)   # n × 1
  cond_mean_vec(
    fit = fit,
    observed_vals_mat = observed_vals_mat,
    observed_idx = 2,
    target_idx = 3
  )
}

PredictMarg = function(modelMarg,dTest){
  result = rep(NA,nrow(dTest))
  result[dTest[["M1"]] == 1] = predict_y_given_x2(modelMarg,
                                                   dTest[dTest[["M1"]] == 1,
                                                         c("X2")])
  result[dTest[["M1"]] == 0] = predict_y_given_x1_x2(modelMarg,
                                                   dTest[dTest[["M1"]] == 0,
                                                         c("X1OBS","X2")])
  return(result)
}

FitMargMI = function(dTrain){
  message("To code ...")
}

PredictMargMI = function(modelMargMI,dTest){
  message("To code ...")
}

FitUnconditionalImputation = function(dTrain, mode = "mean", constant = 0){
  if(!(mode %in% c("mean","median","constant"))) stop("Mode must be \"mean\", \"median\" or \"constant\"")
  impX1 = case_when(
    mode == "mean" ~ mean(dTrain[["X1OBS"]], na.rm = T),
    mode == "median" ~ quantile(dTrain[["X1OBS"]], .5, na.rm = T),
    mode == "constant" ~ constant)
  dTrain[dTrain[["M1"]] == 1,"X1OBS"] = impX1
  
  list(predModel = lm(formula("Y ~ X1OBS + X2"), data = dTrain),
       impModel = impX1)
}

PredictUnconditionalImputation = function(modelUI, dTest){
  dTest[dTest[["M1"]] == 1,"X1OBS"] = modelUI[["impModel"]]
  predict(modelUI[["predModel"]],newdata = dTest)
}

FitSingleConditionalImputation = function(dTrain){
  impModel = lm(formula("X1OBS ~ X2"), data = dTrain)
  dTrain[dTrain[["M1"]] == 1,"X1OBS"] = predict(impModel,
                                                newdata = dTrain[dTrain[["M1"]] == 1,])
  list(predModel = lm(formula("Y ~ X1OBS + X2"), data = dTrain),
       impModel = impModel)
}

PredictSingleConditionalImputation = function(modelSCI,dTest){
  dTest[dTest[["M1"]] == 1,"X1OBS"] = predict(modelSCI[["impModel"]],
                                                newdata = dTest[dTest[["M1"]] == 1,])
  predict(modelSCI[["predModel"]],newdata = dTest[,c("X1OBS","X2")])
}

FitSingleConditionalImputationMI = function(dTrain){
  impModel = lm(formula("X1OBS ~ X2"), data = dTrain[dTrain[["M1"]] == 0,])
  dTrain[dTrain[["M1"]] == 1,"X1OBS"] = predict(impModel,
                                                newdata = dTrain[dTrain[["M1"]] == 1,])
  list(predModel = lm(formula("Y ~ X1OBS*X2*M1"), data = dTrain),
       impModel = impModel)
}

PredictSingleConditionalImputationMI = function(modelSCIMI, dTest){
  dTest[dTest[["M1"]] == 1,"X1OBS"] = predict(modelSCIMI[["impModel"]],
                                              newdata = dTest[dTest[["M1"]] == 1,])
  predict(modelSCIMI[["predModel"]],newdata = dTest[,c("X1OBS","X2","M1")])
}

FitMultipleImputation = function(dTrain, m = 5, method = "norm"){
  imp = mice(dTrain[,c("X1OBS","X2","Y")],
             m = m, method = method, printFlag = F)
  result = list(m = m,
                impModel = list(),
                predModel = list())
  for(impSet in 1:m){
    result[["impModel"]][[impSet]] = lm(formula("X1OBS ~ X2"),
                                        data = complete(imp,impSet))
    result[["predModel"]][[impSet]] = lm(formula("Y ~ X1OBS*X2"),
                                        data = complete(imp,impSet))
  }
  return(result)
}

PredictMultipleImputation = function(modelMI, dTest){
  result = list()
  for(m in 1:modelMI[["m"]]){
    dTest[dTest[["M1"]] == 1,"X1OBS"] = predict(modelMI[["impModel"]][[m]],
                                                newdata = dTest[dTest[["M1"]] == 1,])
    result[[m]] = predict(modelMI[["predModel"]][[m]],
                          newdata = dTest[,c("X1OBS","X2","M1")])
  }
  return(Reduce(`+`, result) / length(result))
}

FitMultipleImputationMI = function(dTrain, m = 5, method = "norm"){
  pm = make.predictorMatrix(dTrain[,c("X1OBS","X2","Y","M1")])
  pm[,"M1"] = 0
  pm["M1",] = 0
  imp = mice(dTrain[,c("X1OBS","X2","Y","M1")], m = m,
             method = method, predictorMatrix = pm, printFlag = F)
  result = list(m = m,
                impModel = list(),
                predModel = list())
  for(impSet in 1:m){
    result[["impModel"]][[impSet]] = lm(formula("X1OBS ~ X2"),
                                        data = complete(imp,impSet)[dTest[["M1"]] == 0,])
    result[["predModel"]][[impSet]] = lm(formula("Y ~ X1OBS*X2*M1"),
                                         data = complete(imp,impSet))
  }
  return(result)
}

PredictMultipleImputationMI = function(modelMIMI, dTest){
  result = list()
  for(m in 1:modelMIMI[["m"]]){
    dTest[dTest[["M1"]] == 1,"X1OBS"] = predict(modelMIMI[["impModel"]][[m]],
                                                newdata = dTest[dTest[["M1"]] == 1,])
    result[[m]] = predict(modelMIMI[["predModel"]][[m]],
                          newdata = dTest[,c("X1OBS","X2","M1")])
  }
  return(Reduce(`+`, result) / length(result))
}