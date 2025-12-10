#-------------------------------------------------------#
#                                                       #
#                       Header                          #
#                                                       #
#-------------------------------------------------------#

# --- Load the scripts ---

source("scripts/R/01_data.R")
source("scripts/R/02_methods.R")
source("scripts/R/03_metrics.R")
source("scripts/R/04_plotting.R")

# --- Set the seed for reproducible results ---

set.seed(42)

# --- Get the number of cores ---
# (simulation of datasets for continuous variables are parallelised)

message("Number of cores: ", detectCores())

# --- Create the directories for the outputs ---

dir.create("outputs/continuous/datasets", recursive = T)
dir.create("outputs/continuous/tables", recursive = T)
dir.create("outputs/continuous/plots", recursive = T)
dir.create("outputs/discrete/datasets", recursive = T)
dir.create("outputs/discrete/tables", recursive = T)
dir.create("outputs/discrete/plots", recursive = T)
dir.create("outputs/Rimages", recursive = T)

#-------------------------------------------------------#
#                                                       #
#                Continuous variables                   #
#                                                       #
#-------------------------------------------------------#

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

# --- Methods tested ---

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

# --- Performance metrics measured ---
metrics = c("MSPE_OMU","MSPE_OMC","MSE")

# --- Initialisation ---

datasetList = list()
predictions = list()
iStorage = 1

# --- Loop ---

for(model in names(phi)){
  startModel = Sys.time()
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
    
    # --- Apply each prediction algorithm ---
    for(methodKey in names(methodsList)){
      # --- Fitting the model on dTrain
      mod = methodsList[[methodKey]][["fit"]](dTrain)
      
      # --- Predict on dTest ---
      predicted = methodsList[[methodKey]][["predict"]](mod,dTest[,c("X1OBS","X2",'M1')])
      predictions[[model]][[i]][[methodKey]] = predicted
      dTest[[methodKey]] = predicted
    }
    
    # --- Store the datasets ---
    dTrain[["SET"]] = "TRAIN"
    dTest[["SET"]] = "TEST"
    dSim = dplyr::bind_rows(dTrain, dTest)
    dSim[["MISSCOEF"]] = missCoef
    dSim[["MODEL"]] = model
    datasetList[[iStorage]] = dSim
    iStorage = iStorage + 1
    
    message("Model: ", model, 
            ", \t MissCoef: ", format(round(missCoef, 3), nsmall = 3), 
            ", \t Set time: ", format(round(Sys.time()-startSet, 3), nsmall = 3))
  }
  message("Model: ", model, ", \t Model time: ",format(round(Sys.time()-startModel, 3), nsmall = 3)) 
}

# --- Save the datasets ---
simulated_datasets = dplyr::bind_rows(datasetList)
# write.csv(simulated_datasets, file = "outputs/continuous/datasets/simulated_datasets_continuous.csv")

#-------------------------------------------------------#
#                                                       #
#                  Continuous variables                 #
#                                                       #
#-------------------------------------------------------#

# --- Compute the performance metrics ---
performanceMetrics = computePerformanceMetrics(predictions,
                                               writeTab =  F)

# --- Plot the results ---
plotPerformanceMetrics(performanceMetrics,
                       methodsKeys = c("MI","MARG","MIMI","MARGMI","PS","CCS"),
                       opacity = .15,
                       savePDF = F)

#-------------------------------------------------------#
#                                                       #
#         Discrete variables - Models 1 to 5            #
#                                                       #
#-------------------------------------------------------#

simulated_datasets_discrete = read.csv("code/outputs/discrete/datasets/saveDataDiscrete.csv")
performanceMetricsDiscrete = dataset2performanceMetrics(simulated_datasets_discrete,
                                                        methodsKeys = c("MARG","MARGMI","PS","MIA","CCS"),
                                                        identifier = "MISSCOEF", missInd = "M")

plotPerformanceMetrics(performanceMetricsDiscrete,
                       opacity = .5,
                       savePDF = F)
#-------------------------------------------------------#
#                                                       #
#         Discrete variables - Models 6 and 7           #
#                                                       #
#-------------------------------------------------------#

# --- extract the results of the simulations in python

simulated_datasets_discrete_M6M7 = read.csv("code/outputs/discrete/datasets/saveDataDiscreteM6M7.csv")
performanceMetricsDiscreteM6M7 = dataset2performanceMetrics(simulated_datasets_discrete_M6M7,
                                                            methodsKeys = c("MARG","MARGMI","PS","MIA","CCS"),
                                                            identifier = "ITER", missInd = "M")
plotPerformanceMetricsM6M7(performanceMetricsDiscreteM6M7,
                           opacity = .5,
                           savePDF = F)

#-------------------------------------------------------#
#                                                       #
#                       Footer                          #
#                                                       #
#-------------------------------------------------------#


# --- save the final images ---
# save.image("outputs/Rimages/output_main_continuous.RData")
