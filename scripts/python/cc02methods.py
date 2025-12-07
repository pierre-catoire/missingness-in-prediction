import pyagrum as gum
import numpy as np
from cc01Data import predictFromEvidence

def fitPatternSubmodels(dTrain):
    dag0 = gum.fastBN("X1OBS->Y<-X2")
    learner = gum.BNLearner(dTrain[dTrain["M1"] == "0"], dag0)
    learner.useSmoothingPrior()
    mod0 = learner.learnParameters(dag0.dag())

    dag1 = gum.fastBN("X2->Y")
    learner = gum.BNLearner(dTrain[dTrain["M1"] == "1"], dag1)
    learner.useSmoothingPrior()
    mod1 = learner.learnParameters(dag1.dag())

    return {"mod0": mod0,
            "mod1": mod1}

def predictPatternSubmodels(psModel, dTest):
    preds = np.empty(dTest.shape[0])

    preds[dTest["M1"] == "0"] = predictFromEvidence(psModel["mod0"], dTest[dTest["M1"] == "0"][["X1OBS","X2"]], "Y")
    preds[dTest["M1"] == "1"] = predictFromEvidence(psModel["mod1"], dTest[dTest["M1"] == "1"][["X2"]], "Y")

    return preds

def fitCCS(dTrain):
    dag0 = gum.fastBN("X1OBS->Y<-X2")
    learner = gum.BNLearner(dTrain[dTrain["M1"] == "0"], dag0)
    mod0 = learner.learnParameters(dag0.dag())

    dag1 = gum.fastBN("X2->Y")
    learner = gum.BNLearner(dTrain, dag1)
    learner.useSmoothingPrior()
    mod1 = learner.learnParameters(dag1.dag())

    return {"mod0": mod0,
            "mod1": mod1}

def predictCCS(ccsModel, dTest):
    preds = np.empty(dTest.shape[0])

    preds[dTest["M1"] == "0"] = predictFromEvidence(ccsModel["mod0"], dTest[dTest["M1"] == "0"][["X1OBS","X2"]], "Y")
    preds[dTest["M1"] == "1"] = predictFromEvidence(ccsModel["mod1"], dTest[dTest["M1"] == "1"][["X2"]], "Y")

    return preds

def fitMIA(dTrain):
    dTrainMIA = dTrain.copy()
    dTrainMIA.loc[dTrainMIA["M1"] == "1", "X1OBS"] = "2"
    dagMIA = gum.fastBN("X1OBS{0|1|2}->Y{0|1}<-X2{0|1}")
    learner = gum.BNLearner(dTrainMIA, dagMIA)
    learner.useSmoothingPrior()
    mod = learner.learnParameters(dagMIA.dag())
    return mod

def predictMIA(MIAmodel, dTest):
    dTestMIA = dTest.copy()
    dTestMIA.loc[dTestMIA["M1"] == "1", "X1OBS"] = "2"
    preds = predictFromEvidence(MIAmodel, dTestMIA[["X1OBS","X2"]], "Y")
    return preds

def fitMarginalisation(dTrain):
    dagMarg = gum.fastBN("X1OBS->Y<-X2")
    learner = gum.BNLearner(dTrain, dagMarg)
    learner.EMsetMaxTime(30000)  
    learner.EMenableMaxTime()
    learner.useSmoothingPrior()
    learner.useEM(1e-5)
    modMarg = learner.learnParameters(dagMarg.dag())
    return modMarg

def predictMarginalisation(margModel, dTest):
    return predictFromEvidence(margModel, dTest[["X1OBS","X2"]], target = "Y")

def fitMarginalisationMI(dTrain):
    dagMargMI = gum.fastBN("X1OBS->Y<-X2->M1;X1OBS->M1<-Y")
    learner = gum.BNLearner(dTrain, dagMargMI)
    learner.EMsetMaxTime(30000)  
    learner.EMenableMaxTime()
    learner.useSmoothingPrior()
    learner.useEM(1e-5)
    modMarg = learner.learnParameters(dagMargMI.dag())
    return modMarg

def predictMarginalisationMI(margModelMI, dTest):
    return predictFromEvidence(margModelMI, dTest[["X1OBS","X2","M1"]], target = "Y")

def predictMarginalisationMIMU(margModelMI, dTest):
    return predictFromEvidence(margModelMI, dTest[["X1OBS","X2"]], target = "Y")