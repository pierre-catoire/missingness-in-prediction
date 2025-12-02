source("R/01_data.R")
source("R/02_methods.R")

# TO REMOVE
theta = list(X1 = list(mu = 0, sigma = 1),
             X2 = list(mu = 0, sigma = 1),
             Y = list(beta = c(0,1,1), sigma = .5))
phi = c(0,1,-1,1)
dTrain = SimulateDataContinuous(1000,theta,phi,missCoef = 0,
                              seed = 42, n_mc = 1e5)
dTest = SimulateDataContinuous(1000,theta,phi,missCoef = 0,
                               seed = 42, n_mc = 1e5)

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

for(methodKey in c("PS","MARG","MARGMI","UI","SCI","SCIMI","MI","MIMI")){
  message(methodKey)
  message("fit ...")
  mod = methodsList[[methodKey]][["fit"]](dTrain)
  message("pred ...")
  pred = methodsList[[methodKey]][["predict"]](mod,dTest)
}
