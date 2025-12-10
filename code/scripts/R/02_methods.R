# 02_methods.R
# Description: Functions to train the models and perform the predictions
# Author: Pierre CATOIRE
# Date: 2025-12-02
# Purpose:
# - train each model on a training set
# - perform the prediction on a testing set

# Preamble
load_package = function(pkg){
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

load_package("tidyverse")
load_package("mice")

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

FitMarginalisation = function(dTrain){
  require(norm)
  dMat = as.matrix(dTrain[,c("X1OBS","X2","Y")])
  s = prelim.norm(dMat)
  theta = em.norm(s,showits=F)
  params = getparam.norm(s,theta)
  names(params[["mu"]]) = c("X1OBS","X2","Y")
  colnames(params[["sigma"]]) = c("X1OBS","X2","Y")
  rownames(params[["sigma"]]) = c("X1OBS","X2","Y")
  return(params)
}

PredictMarginalisation = function(modelMarg, dTest){
  mu = modelMarg[["mu"]]
  sigma_XX = modelMarg[["sigma"]][c("X1OBS","X2"),c("X1OBS","X2")]
  sigma_YX = modelMarg[["sigma"]]["Y",c("X1OBS","X2")]
  
  X_centered = cbind(dTest[["X1OBS"]] - mu["X1OBS"],
                     dTest[["X2"]] - mu["X2"]) 
  
  predicted = as.vector(mu["Y"] + X_centered %*% solve(sigma_XX) %*% sigma_YX)
  
  sigma_YX2 = modelMarg[["sigma"]]["Y","X2"]
  sigma_X2 = modelMarg[["sigma"]]["X2","X2"]
  
  predictedY_X2 = mu[["Y"]] + sigma_YX2 / sigma_X2 * (dTest[["X2"]] - mu[["X2"]])
  predicted[dTest[["M1"]] == 1] = predictedY_X2[dTest[["M1"]] == 1]
  return(predicted)
}

FitMarginalisationMI = function(dTrain){
  if(sum(dTrain[["M1"]]) == 0){
    result = FitMarginalisation(dTrain)
    result[["NO_MISSING_AT_ESTIMATION"]] = T
    return(result)
  }else{
    # Note: need to consider interaction parameters
    require(norm)
    dMat = as.matrix(dTrain[,c("X1OBS","X2","M1","Y")])
    s = prelim.norm(dMat)
    theta = em.norm(s,showits=F)
    params = getparam.norm(s,theta)
    names(params[["mu"]]) = c("X1OBS","X2","M1","Y")
    colnames(params[["sigma"]]) = c("X1OBS","X2","M1","Y")
    rownames(params[["sigma"]]) = c("X1OBS","X2","M1","Y")
    params[["NO_MISSING_AT_ESTIMATION"]] = F
    return(params)
  }
}

PredictMarginalisationMI = function(modelMargMI, dTest){
  if(modelMargMI[["NO_MISSING_AT_ESTIMATION"]]){
    if(sum(dTest[["M1"]]) > 0){
      warning("Marginalisation with missingness indicators:
              \nMissing values were absent at estimation but are present at deployment.
              \nPrediction will be performed without missingness indicator")
    }
    return(PredictMarginalisation(modelMargMI, dTest))
  }else{
    mu     = modelMargMI[["mu"]]
    Sigma  = modelMargMI[["sigma"]]
    
    #### -------------------------------------------------
    #### 1. CASE M1 = 0 → predict E[Y|X1OBS, X2, M1=0]
    #### -------------------------------------------------
    # Variables in order: X1OBS, X2, M1, Y
    
    idx_X1X2M1 = c("X1OBS", "X2", "M1")
    
    # Σ_XX = Cov([X1OBS,X2,M1], [X1OBS,X2,M1])
    sigma_XX = Sigma[idx_X1X2M1, idx_X1X2M1]
    
    # Σ_YX = Cov(Y, [X1OBS,X2,M1])
    sigma_YX = Sigma["Y", idx_X1X2M1]
    
    # Center predictors
    X_centered_full = cbind(
      dTest[["X1OBS"]] - mu["X1OBS"],
      dTest[["X2"]]    - mu["X2"],
      dTest[["M1"]]    - mu["M1"]
    )
    
    # Full conditional prediction
    pred_full = as.vector(mu["Y"] + X_centered_full %*% solve(sigma_XX) %*% sigma_YX)
    
    
    #### -------------------------------------------------
    #### 2. CASE M1 = 1 → predict E[Y|X2,M1=1]
    #### -------------------------------------------------
    # This is marginalised over X1OBS.
    
    idx_X2M1 = c("X2", "M1")
    
    sigma_X2M1 = Sigma[idx_X2M1, idx_X2M1]      # Cov([X2,M1])
    sigma_YX2M1 = Sigma["Y", idx_X2M1]          # Cov(Y,[X2,M1])
    
    X_centered_X2M1 = cbind(
      dTest[["X2"]] - mu["X2"],
      dTest[["M1"]] - mu["M1"]
    )
    
    pred_marg = as.vector(mu["Y"] + 
                            X_centered_X2M1 %*% solve(sigma_X2M1) %*% sigma_YX2M1)
    
    
    #### -------------------------------------------------
    #### 3. Combine predictions
    #### -------------------------------------------------
    predicted = pred_full
    predicted[dTest$M1 == 1] = pred_marg[dTest$M1 == 1]
    
    return(predicted)
  }
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

FitSingleConditionalImputation = function(dTrain, m = 5){
  imp = mice(dTrain[,c("X1OBS","X2","Y")],
             m = m, method = "norm", printFlag = F)
  poolImp = pool(with(data=imp,exp=lm(formula("X1OBS ~ X2"))))
  impModel = setNames(poolImp[,"pooled"][,"estimate"],
                      poolImp[,"pooled"][,"term"])
  
  poolPred = pool(with(data=imp,exp=lm(formula("Y ~ X1OBS*X2"))))
  predModel = setNames(poolPred[,"pooled"][,"estimate"],
                      poolPred[,"pooled"][,"term"])
  
  return(list("impModel" = impModel,
              "predModel" = predModel))
}

PredictSingleConditionalImputation = function(modelSCI, dTest){
  M1 = dTest[["M1"]]
  X2 = dTest[["X2"]]
  X1OBS = dTest[["X1OBS"]]
  X1IMP = ifelse(M1 == 1,
                 cbind(1,X2) %*% modelSCI[["impModel"]],
                 X1OBS)
  pred = as.numeric(cbind(1,X1IMP,X2,X1IMP*X2) %*% modelSCI[["predModel"]])
  return(pred)
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
                impModelWithoutMI = list(), # Model used to impute X1OBS if no
                predModel = list()) 
  for(impSet in 1:m){
    if(any(dTrain[["M1"]] == 1)){
      result[["impModel"]][[impSet]] = lm(formula("X1OBS ~ X2"),
                                          data = complete(imp,impSet)[dTrain[["M1"]] == 1,])
    }else{
      result[["impModel"]][[impSet]] = NULL
    }
    result[["impModelWithoutMI"]][[impSet]] = lm(formula("X1OBS~X2"), data = complete(imp,impSet))
    result[["predModel"]][[impSet]] = lm(formula("Y ~ X1OBS*X2*M1"),
                                         data = complete(imp,impSet))
  }
  return(result)
}


PredictMultipleImputationMI = function(modelMIMI, dTest){
  result = list()
  for(m in 1:modelMIMI[["m"]]){
    if(any(dTest[["M1"]] == 1)){
      if(length(modelMIMI[["impModel"]]) > 0){
        dTest[dTest[["M1"]] == 1,"X1OBS"] = predict(modelMIMI[["impModel"]][[m]],
                                                    newdata = dTest[dTest[["M1"]] == 1,])
      }else{
        warning("Missing values are present at testing but not at estimation.")
        dTest[dTest[["M1"]] == 1,"X1OBS"] = predict(modelMIMI[["impModelWithoutMI"]][[m]],
                                                    newdata = dTest[dTest[["M1"]] == 1,])
      }
    }
    result[[m]] = predict(modelMIMI[["predModel"]][[m]],
                          newdata = dTest[,c("X1OBS","X2","M1")])
  }
  return(Reduce(`+`, result) / length(result))
}
