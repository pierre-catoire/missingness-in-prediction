import pandas as pd # managing databases
import pyagrum as gum # bayesian networks
import time # measuring time of processes
import random as rd # setting seed
import numpy as np
import warnings
import scripts.python.cc01Data as ccData
import scripts.python.cc02Methods as ccMethods
from pathlib import Path
import os

# Creating subfolders for storing outputs
ROOT = Path(__file__).resolve().parent
base = ROOT / "outputs" / "discrete"
subfolders = ["datasets", "plots", "tables"]
for folder in subfolders:
    (base / folder).mkdir(parents=True, exist_ok=True)

rd.seed(42)
warnings.filterwarnings("ignore", category=UserWarning)

modelsDict = {"M1": gum.fastBN("X1{0|1}->Y{0|1}<-X2{0|1};M1{0|1}"),
              "M2": gum.fastBN("X1{0|1}->Y{0|1}<-X2{0|1};M1{0|1}<-X2"),
              "M3": gum.fastBN("X1{0|1}->Y{0|1}<-X2{0|1};M1{0|1}<-X1"),
              "M4": gum.fastBN("X1{0|1}->Y{0|1}<-X2{0|1};X1->M1{0|1}<-Y"),
              "M5": gum.fastBN("X1{0|1}->Y{0|1}<-X2{0|1};M1{0|1}<-Y")}

theta = {"x1":          .5 , # P(X1=1)
         "x2":          .5 , # P(X2=1)
         "y.x1_0.x2_0": .6, # P(Y=1|X1=0,X2=0)
         "y.x1_0.x2_1": .35, # P(Y=1|X1=0,X2=1)
         "y.x1_1.x2_0": .15, # P(Y=1|X1=1,X2=0)
         "y.x1_1.x2_1": .8} # P(Y=1|X1=1,X2=1)

methodsDict = {"PS": {"fit": ccMethods.fitPatternSubmodels, "predict": ccMethods.predictPatternSubmodels},
               "CCS": {"fit": ccMethods.fitCCS, "predict": ccMethods.predictCCS},
               "MIA": {"fit": ccMethods.fitMIA, "predict": ccMethods.predictMIA},
               "MARG": {"fit": ccMethods.fitMarginalisation, "predict": ccMethods.predictMarginalisation},
               "MARGMI": {"fit": ccMethods.fitMarginalisationMI, "predict": ccMethods.predictMarginalisationMI},
               "MARGMIMU": {"fit": ccMethods.fitMarginalisationMI, "predict": ccMethods.predictMarginalisationMIMU}}

dfList = []

for modelLabel in modelsDict.keys():
    startModelTime = time.time()

    for missCoef in np.linspace(0, 0.7, 1):
        startMissCoefTime = time.time()
        mod = modelsDict[modelLabel]
        ccData.setTheta(mod, theta)
        ccData.setPhi(mod, modelLabel, missCoef)

        # Simulate data
        N = 1000
        dTrain, _ = gum.generateSample(mod, N)
        dTrain["X1OBS"] = dTrain["X1"].where(dTrain["M1"] == "0", pd.NA)
        dTrain["SET"] = "TRAIN"
        dTrain["MISSCOEF"] = missCoef
        dTrain["MODEL"] = modelLabel

        dTest, _ = gum.generateSample(mod, N)
        dTest["X1OBS"] = dTest["X1"].where(dTest["M1"] == "0", pd.NA)
        dTest["SET"] = "TEST"
        dTest["MISSCOEF"] = missCoef
        dTest["MODEL"] = modelLabel

        # predict OMU, OMC, PMU, PMC
        probY_X1X2 = ccData.predictFromEvidence(mod, dTest[["X1","X2"]], target = "Y")
        probY_X1X2M1 = ccData.predictFromEvidence(mod, dTest[["X1","X2","M1"]], target = "Y")
        probY_X2 = ccData.predictFromEvidence(mod, dTest[["X2"]], target = "Y")
        probY_X2M1 = ccData.predictFromEvidence(mod, dTest[["X2","M1"]], target = "Y")

        dTest["ORACLE_MU"] = probY_X1X2
        dTest["PRAGMATIC_MU"] = np.where(dTest["M1"] == "0", probY_X1X2, probY_X2)
        dTest["ORACLE_MC"] = probY_X1X2M1
        dTest["PRAGMATIC_MC"] = np.where(dTest["M1"] == "0", probY_X1X2M1, probY_X2M1)

        for methodKey in methodsDict.keys():
            predModel = ccMethods.fitPatternSubmodels(dTrain)
            predModel = methodsDict[methodKey]["fit"](dTrain)
            dTest[methodKey] = methodsDict[methodKey]["predict"](predModel, dTest)
        
        dfList.append(dTrain)
        dfList.append(dTest)
        elapsedMissCoef = time.time() - startMissCoefTime
        print(f"Model {modelLabel}, missCoef {missCoef} time: {elapsedMissCoef:.4f} seconds")
        
    elapsedModel = time.time() - startModelTime
    print(f"Model {modelLabel} time: {elapsedModel:.4f} seconds")
    
dfSim = pd.concat(dfList, axis=0, ignore_index=True)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.join(BASE_DIR, "outputs", "discrete", "datasets")
save_path = os.path.join(save_dir, "saveDataDiscrete.csv")
dfSim.to_csv(save_path, index=False)