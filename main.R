source("01_data.R")
source("02_methods.R")
source("03_metrics.R")
source("04_plotting.R")

library(parallel)

set.seed(42)

message("Number of cores: ", detectCores())

n = 1000

theta = list(X1 = list(mu = 0, sigma = 1),
             X2 = list(mu = 0, sigma = 1),
             Y = list(beta = c(0,1,-1), sigma = 1.3))

phi = list("M1" = c(0,0,0,0),
           "M2" = c(0,0,-1,0),
           "M3" = c(0,-2,0,0),
           "M4" = c(0,-1,0,1),
           "M5" = c(0,0,0,1))

missCoefs = seq(0,0.7,by = 0.01)

methodsList = list("PS" = list("fit" = FitPatternSubmodels,
                               "predict" = PredictPatternSubmodels),
                   "CCS" = list("fit" = FitCompleteCasesSubmodels,
                               "predict" = PredictCompleteCasesSubmodels),
                   "MARG" = list("fit" = FitMarginalisation,
                                 "predict" = PredictMarginalisation),
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

methodsKeys = c("PS","CCS","MARG","UI","SCI","SCIMI","MI","MIMI")
metrics = c("MSPE_OMU","MSPE_OMC","MSE")

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
    for(methodKey in c(methodsKeys,"PRAGMATIC_MU","PRAGMATIC_MC")){
      performanceMetrics[[model]][[metric]][[methodKey]] = numeric()
    }
  }
  
  
  for(i in 1:length(missCoefs)){
    predictions[[model]][[i]] = list()
    missCoef = missCoefs[i]
    message(missCoef)
    start = Sys.time()
    # --- Generate the training and testing datasets ---
    dTrain = SimulateDataContinuous(n,theta,phi[[model]],missCoef, n_mc = 1e4, type = "train")
    dTest = SimulateDataContinuous(n,theta,phi[[model]],missCoef, n_mc = 1e4, type = "test")
    
    # --- Store the realised proportion of missingness in the testing set ---
    performanceMetrics[[model]][["MISSPROP"]][[i]] = mean(dTest[["M1"]])
    
    # --- store the Oracle and Pragmatic probabilities --- #
    predictions[[model]][[i]][["ORACLE_MU"]] = dTest[["ORACLE_MU"]]
    predictions[[model]][[i]][["ORACLE_MC"]] = dTest[["ORACLE_MC"]]
    predictions[[model]][[i]][["PRAGMATIC_MU"]] = dTest[["PRAGMATIC_MU"]]
    predictions[[model]][[i]][["PRAGMATIC_MC"]] = dTest[["PRAGMATIC_MC"]]
    predictions[[model]][[i]][["Y"]] = dTest[["Y"]]
    
    # --- Store the datasets ---
    datasetList[[model]][["TRAIN"]][[i]] = dTrain
    datasetList[[model]][["TEST"]][[i]] = dTest
      
    # --- Apply each prediction algorithm ---
    for(methodKey in methodsKeys){
      # --- Fitting the model on dTrain
      mod = methodsList[[methodKey]][["fit"]](dTrain)
      
      # --- Predict on dTest ---
      predicted = methodsList[[methodKey]][["predict"]](mod,dTest[,c("X1OBS","X2",'M1')])
      
      predictions[[model]][[i]][[methodKey]] = predicted
      performanceMetrics[[model]][["MSPE_OMU"]][[methodKey]][[i]] = ms(predicted,dTest[["ORACLE_MU"]])
      performanceMetrics[[model]][["MSPE_OMC"]][[methodKey]][[i]] = ms(predicted,dTest[["ORACLE_MC"]])
      performanceMetrics[[model]][["MSE"]][[methodKey]][[i]] = ms(predicted,dTest[["Y"]])
    }
    
    for(reference in c("PRAGMATIC_MU","PRAGMATIC_MC")){
      performanceMetrics[[model]][["MSPE_OMU"]][[reference]][[i]] = ms(dTest[[reference]],dTest[["ORACLE_MU"]])
      performanceMetrics[[model]][["MSPE_OMC"]][[reference]][[i]] = ms(dTest[[reference]],dTest[["ORACLE_MC"]])
      performanceMetrics[[model]][["MSE"]][[reference]][[i]] = ms(dTest[[reference]],dTest[["Y"]])
    }
    
   message("Set time: ", Sys.time()-start) 
  }
  message("Model time: ",Sys.time()-start) 
}

save.image("save.RData")

#performanceMetricsTableList = GeneratePerformanceMetricsTables(performanceMetrics)
#simulated_datasets = MergeDatasets(datasetList,missCoefs,predictions)
#write.csv(simulated_datasets, file = "results/datasets/simulated_datasets_continuous.csv")
#plotPerformanceMetrics(performanceMetricsTableList,
#                       methodsKeys = c("SCI","MI","MARG","SCIMI","MIMI","UI","PS","CCS"), opacity = .5)

