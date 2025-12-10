import pyagrum as gum
import pandas as pd
import numpy as np

def predictFromEvidence(model, df, target = "Y"):
    ie = gum.LazyPropagation(model)  # reuse inference engine
    cols = list(df.columns)
    data = df.values              # numpy array for fast access
    n = data.shape[0]
    out = np.zeros(n)
    
    for i in range(n):
        row = data[i]

        # Build evidence dictionary
        evidence = {}
        for j, col in enumerate(cols):

            val = row[j]

            # Skip missing values
            if pd.isna(val):   # pd.NA is not equal to itself
                continue

            evidence[col] = int(val)

        # Perform inference
        ie.setEvidence(evidence)
        ie.makeInference()
        out[i] = ie.posterior(target)[1]

    return out

def fitPatternSubmodelsM6M7(dTrain):
    dag0 = gum.fastBN("X1OBS->Y<-X2OBS")
    dag1 = gum.fastBN("X2OBS->Y")
    dag2 = gum.fastBN("X1OBS->Y")

    subModels = {}

    if dTrain["M"].eq("0").any():
        learner = gum.BNLearner(dTrain[dTrain["M"] == "0"], dag0)
        learner.useSmoothingPrior()
        subModels["mod0"] = learner.learnParameters(dag0.dag())

    if dTrain["M"].eq("1").any():
        learner = gum.BNLearner(dTrain[dTrain["M"] == "1"], dag1)
        learner.useSmoothingPrior()
        subModels["mod1"] = learner.learnParameters(dag1.dag())
    else:
        subModels["mod1"] = subModels["mod0"]

    if dTrain["M"].eq("2").any():
        learner = gum.BNLearner(dTrain[dTrain["M"] == "2"], dag2)
        learner.useSmoothingPrior()
        subModels["mod2"] = learner.learnParameters(dag2.dag())
    else:
        subModels["mod2"] = subModels["mod0"]  

    return subModels

def predictPatternSubmodelsM6M7(psModel, dTest):
    preds = np.empty(dTest.shape[0])

    preds[dTest["M"] == "0"] = predictFromEvidence(psModel["mod0"], dTest[dTest["M"] == "0"][["X1OBS","X2OBS"]], "Y")
    preds[dTest["M"] == "1"] = predictFromEvidence(psModel["mod1"], dTest[dTest["M"] == "1"][["X2OBS"]], "Y")
    preds[dTest["M"] == "2"] = predictFromEvidence(psModel["mod2"], dTest[dTest["M"] == "2"][["X1OBS"]], "Y")

    return preds

def fitCCS(dTrain):
    dag0 = gum.fastBN("X1OBS->Y<-X2OBS")
    dag1 = gum.fastBN("X2OBS->Y")
    dag2 = gum.fastBN("X1OBS->Y")

    subModels = {}

    if dTrain["M"].eq("0").any():
        learner = gum.BNLearner(dTrain[dTrain["M"] == "0"], dag0)
        learner.useSmoothingPrior()
        subModels["mod0"] = learner.learnParameters(dag0.dag())

    if dTrain["M"].eq("1").any():
        learner = gum.BNLearner(dTrain[dTrain["M"] != "2"], dag1)
        learner.useSmoothingPrior()
        subModels["mod1"] = learner.learnParameters(dag1.dag())
    else:
        subModels["mod1"] = subModels["mod0"]

    if dTrain["M"].eq("2").any():
        learner = gum.BNLearner(dTrain[dTrain["M"] != "1"], dag2)
        learner.useSmoothingPrior()
        subModels["mod2"] = learner.learnParameters(dag2.dag())
    else:
        subModels["mod2"] = subModels["mod0"]  

    return subModels

def predictCCS(ccsModel, dTest):
    preds = np.empty(dTest.shape[0])

    preds[dTest["M"] == "0"] = predictFromEvidence(ccsModel["mod0"], dTest[dTest["M"] == "0"][["X1OBS","X2OBS"]], "Y")
    preds[dTest["M"] == "1"] = predictFromEvidence(ccsModel["mod1"], dTest[dTest["M"] == "1"][["X2OBS"]], "Y")
    preds[dTest["M"] == "2"] = predictFromEvidence(ccsModel["mod2"], dTest[dTest["M"] == "2"][["X1OBS"]], "Y")

    return preds

def fitMIA(dTrain):
    dTrainMIA = dTrain.copy()
    dTrainMIA.loc[dTrainMIA["M"] == "1", "X1OBS"] = "2"
    dTrainMIA.loc[dTrainMIA["M"] == "2", "X2OBS"] = "2"
    dagMIA = gum.fastBN("X1OBS{0|1|2}->Y{0|1}<-X2OBS{0|1|2}")
    learner = gum.BNLearner(dTrainMIA, dagMIA)
    learner.useSmoothingPrior()
    mod = learner.learnParameters(dagMIA.dag())
    return mod

def predictMIA(MIAmodel, dTest):
    dTestMIA = dTest.copy()
    dTestMIA.loc[dTestMIA["M"] == "1", "X1OBS"] = "2"
    dTestMIA.loc[dTestMIA["M"] == "2", "X2OBS"] = "2"
    preds = predictFromEvidence(MIAmodel, dTestMIA[["X1OBS","X2OBS"]], "Y")
    return preds

def fitMarginalisation(dTrain):
    dagMarg = gum.fastBN("X1OBS->Y<-X2OBS")
    learner = gum.BNLearner(dTrain, dagMarg)
    learner.EMsetMaxTime(30000)  
    learner.EMenableMaxTime()
    learner.useSmoothingPrior()
    learner.useEM(1e-5)
    modMarg = learner.learnParameters(dagMarg.dag())
    return modMarg

def predictMarginalisation(margModel, dTest):
    return predictFromEvidence(margModel, dTest[["X1OBS","X2OBS"]], target = "Y")

def fitMarginalisationMI(dTrain):
    dagMargMI = gum.fastBN("X1OBS->Y<-X2OBS->M{0|1|2};X1OBS->M<-Y")
    learner = gum.BNLearner(dTrain, dagMargMI)
    learner.EMsetMaxTime(30000)  
    learner.EMenableMaxTime()
    learner.useSmoothingPrior()
    learner.useEM(1e-5)
    modMarg = learner.learnParameters(dagMargMI.dag())
    return modMarg

def predictMarginalisationMI(margModelMI, dTest):
    return predictFromEvidence(margModelMI, dTest[["X1OBS","X2OBS","M"]], target = "Y")

def predictMarginalisationMIMU(margModelMI, dTest):
    return predictFromEvidence(margModelMI, dTest[["X1OBS","X2OBS"]], target = "Y")