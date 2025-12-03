source("R/01_data.R")
source("R/02_methods.R")
source("R/03_metrics.R")
source("R/04_plotting.R")

n = 20

theta = list(X1 = list(mu = 0, sigma = 1),
             X2 = list(mu = 0, sigma = 1),
             Y = list(beta = c(0,1,1), sigma = .5))

phi = list("M1" = c(0,0,0,0),
           "M2" = c(0,0,1,0),
           "M3" = c(0,1,0,0),
           "M4" = c(0,1,0,1),
           "M5" = c(0,0,0,1))

missCoefs = seq(0,0.7,by = 0.1)

methodsList = list("PS" = list("fit" = FitPatternSubmodels,
                               "predict" = PredictPatternSubmodels),
                   "MARG" = list("fit" = fitMarg,
                                 "predict" = PredictMarg),
                   "MARGMI" = list("fit" = FitMargMI,
                                   "predict" = PredictMargMI),
                   "UI" = list("fit" = FitUnconditionalImputation,
                               "predict" = PredictUnconditionalImputation),
                   "SCI" = list("fit" = FitSingleConditionalImputation,
                                "predict" = PredictSingleConditionalImputation),
                   "SCIMI" = list("fit" = FitSingleConditionalImputationMI,
                                  "predict" = PredictSingleConditionalImputationMI),
                   "MI" = list("fit" = FitMultipleImputation,
                               "predict" = PredictMultipleImputation),
                   "MIMI" = list("fit" = FitMultipleImputationMI,
                                 "predict" = PredictMultipleImputationMI))

methodsKeys = c("PS","MARG","MARGMI","UI","SCI","SCIMI","MI","MIMI")
metrics = c("MSPE-OMU","MSPE-OMC","MSE")

datasetList = list()
predictions = list()
performanceMetrics = list()

for(model in names(phi)){
  message(model)
  datasetList[[model]] = list("TRAIN" = list(),
                              "TEST" = list())
  
  predictions[[model]] = list()
  performanceMetrics[[model]] = list("MISSPROP" = numeric())
  for(metric in metrics){
    performanceMetrics[[model]][[metric]] = list()
    for(methodKey in methodsKeys){
      performanceMetrics[[model]][[metric]][[methodKey]] = numeric()
    }
  }
  
  
  for(i in 1:length(missCoefs)){
    predictions[[model]][[i]] = list()
    missCoef = missCoefs[i]
    message(missCoef)
    
    # --- Generate the training and testing datasets ---
    dTrain = SimulateDataContinuous(n,theta,phi[[model]],missCoef, seed = 42, n_mc = 1e5, type = "train")
    dTest = SimulateDataContinuous(n,theta,phi[[model]],missCoef, seed = 42, n_mc = 1e5, type = "test")
    
    # --- Store the datasets ---
    datasetList[[model]][["TRAIN"]][[i]] = dTrain
    datasetList[[model]][["TEST"]][[i]] = dTest
    
    # --- Store the realised proportion of missingness in the testing set ---
    performanceMetrics[[model]][["MISSPROP"]][[i]] = mean(dTest[["M1"]])
    
    # --- store the Oracle and Pragmatic probabilities --- #
    predictions[[model]][[i]][["ORACLE_MU"]] = dTest[["ORACLE_MU"]]
    predictions[[model]][[i]][["ORACLE_MC"]] = dTest[["ORACLE_MC"]]
    predictions[[model]][[i]][["PRAGMATIC_MU"]] = dTest[["PRAGMATIC_MU"]]
    predictions[[model]][[i]][["PRAGMATIC_MC"]] = dTest[["PRAGMATIC_MC"]]
    predictions[[model]][[i]][["Y"]] = dTest[["Y"]]
      
    # --- Apply each prediction algorithm ---
    for(methodKey in methodsKeys){
      # --- Fitting the model on dTrain
      mod = methodsList[[methodKey]][["fit"]](dTrain)
      
      # --- Predict on dTest ---
      predicted = methodsList[[methodKey]][["predict"]](mod,
                                                        dTest[,c("X1OBS","X2",'M1')])
      
      predictions[[model]][[i]][[methodKey]] = predicted
      performanceMetrics[[model]][["MSPE-OMU"]][[methodKey]][[i]] = ms(predicted,dTest[["ORACLE_MU"]])
      performanceMetrics[[model]][["MSPE-OMC"]][[methodKey]][[i]] = ms(predicted,dTest[["ORACLE_MC"]])
      performanceMetrics[[model]][["MSE"]][[methodKey]][[i]] = ms(predicted,dTest[["Y"]])
    }
  }
}

GeneratePerformanceMetricsTables(performanceMetrics)
simulated_datasets = MergeDatasets(datasetList,missCoefs,predictions)
write.csv(simulated_datasets, file = "results/datasets/simulated_datasets_continuous.csv")
