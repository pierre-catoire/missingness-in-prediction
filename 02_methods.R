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

FitCompleteCasesSubmodels = function(dTrain){
  modelCCS = list()
  if(any(dTrain[["M1"]] == 0)){
    modelCCS[["M0"]] = lm(formula("Y ~ X1OBS + X2"),
                          data = dTrain[dTrain[["M1"]] == 0,])
  }
  if(any(dTrain[["M1"]] == 1)){
    modelCCS[["M1"]] = lm(formula("Y ~ X2"),
                          data = dTrain)
  }
  return(modelCCS)
}

PredictCompleteCasesSubmodels = function(modelCCS,dTest){
  result = rep(NA,nrow(dTest))
  if(any(dTest[["M1"]] == 0)){
    if(is.null(modelCCS[["M0"]])){
      warning("Model does not contain a submodel for M1 == 0, and some rows in dTest have M1 == 0.
            \nPrediction is NA for these cases.")
    }
    else{
      result[dTest[["M1"]] == 0] = predict(modelCCS[["M0"]],
                                           newdata = dTest[dTest[["M1"]] == 0,])
    }
  }
  if(any(dTest[["M1"]] == 1)){
    if(is.null(modelCCS[["M1"]])){
      warning("Model does not contain a submodel for M1 == 1, and some rows in dTest have M1 == 1.
            \nPrediction is NA for these cases.")
    }
    else{
      result[dTest[["M1"]] == 1] = predict(modelCCS[["M1"]],
                                           newdata = dTest[dTest[["M1"]] == 1,])
    }
  }
  return(result)
}

FitMarginalisation = function(dTrain, tol = 1e-8, max_iter = 100) {
  
  X1OBS = dTrain[["X1OBS"]]
  X2 = dTrain[["X2"]]
  Y = dTrain[["Y"]]
  M1 = dTrain[["M1"]]
  # Initialize parameters
  thetaEM = list(
    X1 = list(mu = mean(X1OBS, na.rm = TRUE),
              sigma = sd(X1OBS, na.rm = TRUE)),
    Y  = list(beta = rep(0, 4), sigma = 1)
  )
  
  # Initial regression using complete cases
  obs_idx = which(!is.na(X1OBS))
  modY = lm(Y[obs_idx] ~ X1OBS[obs_idx]*X2[obs_idx])
  thetaEM$Y$beta = coef(modY)
  thetaEM$Y$sigma = sigma(modY)
  
  iter = 1
  converged = FALSE
  
  while(iter <= max_iter && !converged){
    
    thetaEM_old = thetaEM
    
    # --- E-step: expected X1 for missing values ---
    miss_idx = M1 == 1
    if(length(miss_idx) > 0){
      # Posterior of X1 | Y, X2 for missing
      beta0 = thetaEM$Y$beta[1]
      beta1 = thetaEM$Y$beta[2]
      beta2 = thetaEM$Y$beta[3]
      sigma_y = thetaEM$Y$sigma
      mu_X1 = thetaEM$X1$mu
      sigma_X1 = thetaEM$X1$sigma
      
      # Standard normal conditioning formula:
      # X1 | Y, X2 ~ N(mean_post, var_post)
      var_post = 1 / (1/sigma_X1^2 + beta1^2 / sigma_y^2)
      mean_post = var_post * (mu_X1 / sigma_X1^2 + beta1 * (Y[miss_idx] - beta0 - beta2*X2[miss_idx]) / sigma_y^2)
      
      X1OBS[miss_idx] = mean_post
    }
    
    # --- M-step: update parameters using completed data ---
    modY = lm(Y ~ X1OBS*X2)
    thetaEM$Y$beta = coef(modY)
    thetaEM$Y$sigma = sigma(modY)
    
    thetaEM$X1$mu = mean(X1OBS)
    thetaEM$X1$sigma = sd(X1OBS)
    
    # --- Convergence check ---
    max_change = max(
      abs(thetaEM$X1$mu - thetaEM_old$X1$mu) / abs(thetaEM_old$X1$mu),
      abs(thetaEM$X1$sigma - thetaEM_old$X1$sigma) / abs(thetaEM_old$X1$sigma),
      max(abs(thetaEM$Y$beta - thetaEM_old$Y$beta) / abs(thetaEM_old$Y$beta)),
      abs(thetaEM$Y$sigma - thetaEM_old$Y$sigma) / abs(thetaEM_old$Y$sigma)
    )
    
    if(max_change < tol) converged = TRUE
    iter = iter + 1
  }
  
  if(!converged) warning("EM did not converge")
  return(list("Y" = modY,
              "X1" = mean(X1OBS),
              "iter" = iter,
              "converged" = converged))
}

PredictMarginalisation = function(modelMarg,dTest){
  result = rep(NA,nrow(dTest))
  # --- E(Y|X1,X2) ---
  result[dTest[["M1"]] == 0] = predict(modelMarg[["Y"]],dTest[dTest[["M1"]]== 0,])
  
  # --- E(Y|X2) = B0 + B1mu1 + B2X2 + B3X2mu1 ---
  if(any(dTest[["M1"]] == 1)){
    result[dTest[["M1"]] == 1] = cbind(1,
                                       modelMarg[["X1"]],
                                       dTest[dTest[["M1"]] == 1,"X2"],
                                       modelMarg[["X1"]]*dTest[dTest[["M1"]] == 1,"X2"]) %*% coef(modelMarg[["Y"]])
  }
  return(result)
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
  impModelEstimation = lm(X1OBS ~ X2*Y, data = dTrain[dTrain[["M1"]] == 0,])
  dTrain[["X1IMP"]] = dTrain[["X1OBS"]]
  dTrain[dTrain[["M1"]] == 1,"X1IMP"] = predict(impModelEstimation,
                                                newdata = dTrain[dTrain[["M1"]] == 1,])
  
  ## Fit the imputation function
  impModelDeployment = lm(X1IMP ~ X2, data = dTrain)
  
  ## Fit the prediction model
  predModel = lm(Y ~ X1IMP*X2, data = dTrain)
  return(list("impModel" = impModelDeployment,
              "predModel" = predModel))
}



PredictSingleConditionalImputation = function(modelSCI, dTest){
  dTest[["X1IMP"]] = dTest[["X1OBS"]]
  dTest[dTest[["M1"]] == 1,"X1IMP"] = predict(modelSCI[["impModel"]],
                                              newdata = dTest[dTest[["M1"]] == 1,])
  return(predict(modelSCI[["predModel"]], newdata = dTest))
}


FitSingleConditionalImputationMI = function(dTrain){
  impModelEstimationMI = lm(X1OBS ~ X2*Y, data = dTrain[dTrain[["M1"]] == 0,])
  dTrain[["X1IMP"]] = dTrain[["X1OBS"]]
  dTrain[dTrain[["M1"]] == 1,"X1IMP"] = predict(impModelEstimationMI,
                                                newdata = dTrain[dTrain[["M1"]] == 1,])
  
  ## Fit the imputation function
  impModelDeploymentMI = lm(X1IMP ~ X2*M1, data = dTrain,)
  
  ## Fit the prediction model
  predModelMI = lm(Y ~ X1IMP*X2*M1, data = dTrain)
  
  return(list(impModel = impModelDeploymentMI,
              predModel = predModelMI))
}



PredictSingleConditionalImputationMI = function(modelSCIMI, dTest){
  dTest[["X1IMP"]] = dTest[["X1OBS"]]
  dTest[dTest[["M1"]] == 1,"X1IMP"] = predict(modelSCIMI[["impModel"]],
                                              newdata = dTest[dTest[["M1"]] == 1,])
  # Predict Y
  YPREDSCIMI = predict(modelSCIMI[["predModel"]], newdata = dTest)
  return(YPREDSCIMI)
}



FitMultipleImputation = function(dTrain, m = 5, method = "norm"){
  # --- impute missing X1 from X2,Y in training dataset ---
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
                          newdata = dTest[,c("X1OBS","X2")])
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
    if(any(dTrain[["M1"]] == 1)){
      result[["impModel"]][[impSet]] = lm(formula("X1OBS ~ X2"),
                                          data = complete(imp,impSet)[dTrain[["M1"]] == 1,])
    }else{
      result[["impModel"]][[impSet]] = NULL
    }
    result[["predModel"]][[impSet]] = lm(formula("Y ~ X1OBS*X2*M1"),
                                         data = complete(imp,impSet))
  }
  return(result)
}

PredictMultipleImputationMI = function(modelMIMI, dTest){
  result = list()
  for(m in 1:modelMIMI[["m"]]){
    if(any(dTest[["M1"]] == 1)){
      if(is.null(modelMIMI[["impModel"]][[m]])){
        warning("Missing values are present at testing but not at estimation.
              \n No prediction for individuals with missing values is possible.
              \n NA will be returned")
      }
      dTest[dTest[["M1"]] == 1,"X1OBS"] = predict(modelMIMI[["impModel"]][[m]],
                                                  newdata = dTest[dTest[["M1"]] == 1,])
    }
    result[[m]] = predict(modelMIMI[["predModel"]][[m]],
                          newdata = dTest[,c("X1OBS","X2","M1")])
  }
  return(Reduce(`+`, result) / length(result))
}