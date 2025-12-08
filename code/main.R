source("scripts/R/01_data.R")
source("scripts/R/02_methods.R")
source("scripts/R/03_metrics.R")
source("scripts/R/04_plotting.R")

library(parallel)

set.seed(42)

message("Number of cores: ", detectCores())

# --- Create the directories for the outputs ---

dir.create("outputs/continuous/datasets", recursive = T)
dir.create("outputs/continuous/tables", recursive = T)
dir.create("outputs/continuous/plots", recursive = T)
dir.create("outputs/continuous/Rimages", recursive = T)

# --- Simulation parameters ---

n = 1000

theta = list(X1 = list(mu = 0, sigma = 1),
             X2 = list(mu = 0, sigma = 1),
             Y = list(beta = c(0,1,-1), sigma = 1.3))

phi = list("M1" = c(0,0,0,0),
           "M2" = c(0,0,-1,0),
           "M3" = c(0,-2,0,0),
           "M4" = c(0,2,0,-2),
           "M5" = c(0,0,0,1))

missCoefs = seq(0,0.7,by = 0.001)

methodsList = list("PS" = list("fit" = FitPatternSubmodels,
                               "predict" = PredictPatternSubmodels),
                   "CCS" = list("fit" = FitCompleteCasesSubmodels,
                                "predict" = PredictCompleteCasesSubmodels),
                   "MARG" = list("fit" = FitMarginalisation,
                                 "predict" = PredictMarginalisation),
                   "MARGMI" = list("fit" = FitMarginalisationMI,
                                   "predict" = PredictMarginalisationMI),
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

metrics = c("MSPE_OMU","MSPE_OMC","MSE")

datasetList = list()
predictions = list()

for(model in names(phi)){
  startModel = Sys.time()
  datasetList[[model]] = list("TRAIN" = list(),
                              "TEST" = list())
  
  predictions[[model]] = list()
  
  for(i in 1:length(missCoefs)){
    predictions[[model]][[i]] = list()
    missCoef = missCoefs[i]
    startSet = Sys.time()
    # --- Generate the training and testing datasets ---
    dTrain = SimulateDataContinuous(n,theta,phi[[model]],missCoef, n_mc = 1e4, type = "train")
    dTest = SimulateDataContinuous(n,theta,phi[[model]],missCoef, n_mc = 1e4, type = "test")
    
    # --- store the Reference Probabilities, the outcome and the missingness indicator --- #
    for(var in c("ORACLE_MU","ORACLE_MC","PRAGMATIC_MU","PRAGMATIC_MC","Y","M1")){
      predictions[[model]][[i]][[var]] = dTest[[var]]
    }
    
    # --- Store the datasets ---
    datasetList[[model]][["TRAIN"]][[i]] = dTrain
    datasetList[[model]][["TEST"]][[i]] = dTest
    
    # --- Apply each prediction algorithm ---
    for(methodKey in names(methodsList)){
      # --- Fitting the model on dTrain
      mod = methodsList[[methodKey]][["fit"]](dTrain)
      
      # --- Predict on dTest ---
      predicted = methodsList[[methodKey]][["predict"]](mod,dTest[,c("X1OBS","X2",'M1')])
      predictions[[model]][[i]][[methodKey]] = predicted
    }
    message("Model: ", model, 
            ", \t MissCoef: ", format(round(missCoef, 3), nsmall = 3), 
            ", \t Set time: ", format(round(Sys.time()-startSet, 3), nsmall = 3))
  }
  message("Model: ", model, ", \t Model time: ",format(round(Sys.time()-startModel, 3), nsmall = 3)) 
}

# --- Save the datasets ---
simulated_datasets = MergeDatasets(datasetList,missCoefs)
write.csv(simulated_datasets, file = "outputs/continuous/datasets/simulated_datasets_continuous.csv")

# --- Compute the performance metrics ---
performanceMetrics = computePerformanceMetrics(predictions,
                                               writeTab =  F)

# --- plot the results ---
plotPerformanceMetrics(performanceMetrics,
                       methodsKeys = c("MI","MARG","MIMI","MARGMI","PS","CCS"),
                       opacity = .15,
                       savePDF = F)

save.image("outputs/continuous/Rimages/output_main_continuous.RData")