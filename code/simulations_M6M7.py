import pandas as pd # managing databases
import pyagrum as gum # bayesian networks
import time # measuring time of processes
import random as rd # setting seed
import numpy as np
import warnings
import scripts.python.cc01Data as ccData
import scripts.python.cc03MethodsM6M7 as ccMethodsM6M7
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

# --- Set the model M6 ---

m6 = gum.fastBN("X1{0|1}->X2{0|1}->Y{0|1}->M{0|1|2};X1->Y;X1->M<-X2")

m6.cpt("X1").fillWith([51/105,54/105])

m6.cpt("X2")[{"X1":0}] = [24/51,27/51]
m6.cpt("X2")[{"X1":1}] = [30/54,24/54]

m6.cpt("Y")[{"X1":0,"X2":0}] = [18/24,6/24]
m6.cpt("Y")[{"X1":0,"X2":1}] = [18/27,9/27]
m6.cpt("Y")[{"X1":1,"X2":0}] = [18/30,12/30]
m6.cpt("Y")[{"X1":1,"X2":1}] = [9/24,15/24]

m6.cpt("M")[{"X1":0,"X2":0,"Y":0}] = [3/18, 7/18, 8/18]
m6.cpt("M")[{"X1":0,"X2":0,"Y":1}] = [1/6, 3/6, 2/6]
m6.cpt("M")[{"X1":0,"X2":1,"Y":0}] = [6/18, 4/18, 8/18]
m6.cpt("M")[{"X1":0,"X2":1,"Y":1}] = [3/9, 3/9, 3/9]
m6.cpt("M")[{"X1":1,"X2":0,"Y":0}] = [3/18, 7/18, 8/18]
m6.cpt("M")[{"X1":1,"X2":0,"Y":1}] = [2/12, 6/12, 4/12]
m6.cpt("M")[{"X1":1,"X2":1,"Y":0}] = [3/9, 2/9, 4/9]
m6.cpt("M")[{"X1":1,"X2":1,"Y":1}] = [5/15, 5/15, 5/15]

# --- Set the model M7 ---

m7 = gum.fastBN("X1{0|1}->X2{0|1}->Y{0|1}->M{0|1|2};X1->Y;X1->M<-X2")

m7.cpt("X1").fillWith([55/105,50/105])

m7.cpt("X2")[{"X1":0}] = [28/55,27/55]
m7.cpt("X2")[{"X1":1}] = [20/50,30/50]

m7.cpt("Y")[{"X1":0,"X2":0}] = [14/28,14/28]
m7.cpt("Y")[{"X1":0,"X2":1}] = [6/27,21/27]
m7.cpt("Y")[{"X1":1,"X2":0}] = [4/20,16/20]
m7.cpt("Y")[{"X1":1,"X2":1}] = [21/30,9/30]

m7.cpt("M")[{"X1":0,"X2":0,"Y":0}] = [3/14,4/14,7/14]
m7.cpt("M")[{"X1":0,"X2":0,"Y":1}] = [3/14,3/14,8/14]
m7.cpt("M")[{"X1":0,"X2":1,"Y":0}] = [2/6,3/6,1/6]
m7.cpt("M")[{"X1":0,"X2":1,"Y":1}] = [7/21,8/21,6/21]
m7.cpt("M")[{"X1":1,"X2":0,"Y":0}] = [1/4,2/4,1/4]
m7.cpt("M")[{"X1":1,"X2":0,"Y":1}] = [4/16,7/16,5/16]
m7.cpt("M")[{"X1":1,"X2":1,"Y":0}] = [7/21,6/21,8/21]
m7.cpt("M")[{"X1":1,"X2":1,"Y":1}] = [3/9,2/9,4/9]

modelsDict = {"M6": m6,
              "M7": m7}

methodsDict = {"PS": {"fit": ccMethodsM6M7.fitPatternSubmodelsM6M7, "predict": ccMethodsM6M7.predictPatternSubmodelsM6M7},
               "CCS": {"fit": ccMethodsM6M7.fitCCS, "predict": ccMethodsM6M7.predictCCS},
               "MIA": {"fit": ccMethodsM6M7.fitMIA, "predict": ccMethodsM6M7.predictMIA},
               "MARG": {"fit": ccMethodsM6M7.fitMarginalisation, "predict": ccMethodsM6M7.predictMarginalisation},
               "MARGMI": {"fit": ccMethodsM6M7.fitMarginalisationMI, "predict": ccMethodsM6M7.predictMarginalisationMI},
               "MARGMIMU": {"fit": ccMethodsM6M7.fitMarginalisationMI, "predict": ccMethodsM6M7.predictMarginalisationMIMU}}

dfList = []

nSim = 1
N = 1000

for modelLabel in modelsDict.keys():
    startModelTime = time.time()

    for i in range(nSim):
        startMissCoefTime = time.time()
        mod = modelsDict[modelLabel]

        # Simulate data
        dTrain, _ = gum.generateSample(mod, N)
        dTrain["X1OBS"] = dTrain["X1"].where(dTrain["M"] != "1", pd.NA)
        dTrain["X2OBS"] = dTrain["X2"].where(dTrain["M"] != "2", pd.NA)
        dTrain["SET"] = "TRAIN"
        dTrain["ITER"] = i
        dTrain["MODEL"] = modelLabel

        dTest, _ = gum.generateSample(mod, N)
        dTest["X1OBS"] = dTest["X1"].where(dTest["M"] != "1", pd.NA)
        dTest["X2OBS"] = dTest["X2"].where(dTest["M"] != "2", pd.NA)
        dTest["SET"] = "TEST"
        dTest["ITER"] = i
        dTest["MODEL"] = modelLabel

        # predict OMU, OMC, PMU, PMC
        probY_X1X2 = ccMethodsM6M7.predictFromEvidence(mod, dTest[["X1","X2"]], target = "Y") # vector of 1000
        probY_X1X2M = ccMethodsM6M7.predictFromEvidence(mod, dTest[["X1","X2","M"]], target = "Y")
        probY_X1 = ccMethodsM6M7.predictFromEvidence(mod, dTest[["X1"]], target = "Y")
        probY_X1M = ccMethodsM6M7.predictFromEvidence(mod, dTest[["X1","M"]], target = "Y")
        probY_X2 = ccMethodsM6M7.predictFromEvidence(mod, dTest[["X2"]], target = "Y")
        probY_X2M = ccMethodsM6M7.predictFromEvidence(mod, dTest[["X2","M"]], target = "Y")
        conditions = [dTest["M"] == "0",
                      dTest["M"] == "1",
                      dTest["M"] == "2"]
        
        choicesMU = [probY_X1X2,
                     probY_X2,
                     probY_X1]
        
        choicesMC = [probY_X1X2M,
                     probY_X2M,
                     probY_X1M]
        
        dTest["ORACLE_MU"] = probY_X1X2
        dTest["ORACLE_MC"] = probY_X1X2M
        dTest["PRAGMATIC_MU"] = np.select(conditions, choicesMU)
        dTest["PRAGMATIC_MC"] = np.select(conditions, choicesMC)

        for methodKey in methodsDict.keys():
            predModel = methodsDict[methodKey]["fit"](dTrain)
            dTest[methodKey] = methodsDict[methodKey]["predict"](predModel, dTest)
        
        dfList.append(dTrain)
        dfList.append(dTest)
        elapsedMissCoef = time.time() - startMissCoefTime
        print(f"Model {modelLabel}, iter {i} time: {elapsedMissCoef:.4f} seconds")
        
    elapsedModel = time.time() - startModelTime
    print(f"Model {modelLabel} time: {elapsedModel:.4f} seconds")
    
dfSim = pd.concat(dfList, axis=0, ignore_index=True)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
save_dir = os.path.join(BASE_DIR, "outputs", "discrete", "datasets")
save_path = os.path.join(save_dir, "saveDataDiscreteM6M7.csv")
dfSim.to_csv(save_path, index=False)