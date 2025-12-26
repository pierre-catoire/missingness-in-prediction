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

result = [1, 1, 1, 3, 1, 5, 6, 1, 8, 2, 1, 2, 2, 2, 8, 1, 2, 3, 6, 8, 1, 6, 5, 4]

a = result[0]
b = result[1]
c = result[2]
d = result[3]
e = result[4]
f = result[5]
g = result[6]
h = result[7]
i = result[8]
j = result[9]
k = result[10]
l = result[11]
m = result[12]
n = result[13]
o = result[14]
p = result[15]
q = result[16]
r = result[17]
s = result[18]
t = result[19]
u = result[20]
v = result[21]
w = result[22]
x = result[23]

mod = gum.fastBN("X1{0|1}->Y{0|1}<-X2{0|1};X1->M{0|1|2}<-X2;X1->X2;Y->M")
pX1 = (a+b+c+d+e+f+g+h+i+j+k+l)/(a+b+c+d+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t+u+v+w+x)
mod.cpt("X1").fillWith([pX1,1-pX1])
mod.cpt("X2")[{"X1":0}] = [(a+b+c+d+e+f)/(a+b+c+d+e+f+g+h+i+j+k+l),1-(a+b+c+d+e+f)/(a+b+c+d+e+f+g+h+i+j+k+l)]
mod.cpt("X2")[{"X1":1}] = [(m+n+o+p+q+r)/(m+n+o+p+q+r+s+t+u+v+w+x),1-(m+n+o+p+q+r)/(m+n+o+p+q+r+s+t+u+v+w+x)]
mod.cpt("Y")[{"X1":0,"X2":0}] = [(a+b+c)/(a+b+c+d+e+f),1-(a+b+c)/(a+b+c+d+e+f)]
mod.cpt("Y")[{"X1":0,"X2":1}] = [(g+h+i)/(g+h+i+j+k+l),1-(g+h+i)/(g+h+i+j+k+l)]
mod.cpt("Y")[{"X1":1,"X2":0}] = [(m+n+o)/(m+n+o+p+q+r),1-(m+n+o)/(m+n+o+p+q+r)]
mod.cpt("Y")[{"X1":1,"X2":1}] = [(s+t+u)/(s+t+u+v+w+x),1-(s+t+u)/(s+t+u+v+w+x)]
mod.cpt("M")[{"X1":0,"X2":0, "Y":0}] = [a/(a+b+c), b/(a+b+c), c/(a+b+c)]
mod.cpt("M")[{"X1":0,"X2":0, "Y":1}] = [d/(d+e+f), e/(d+e+f), f/(d+e+f)]
mod.cpt("M")[{"X1":0,"X2":1, "Y":0}] = [g/(g+h+i), h/(g+h+i), i/(g+h+i)]
mod.cpt("M")[{"X1":0,"X2":1, "Y":1}] = [j/(j+k+l), k/(j+k+l), l/(j+k+l)]
mod.cpt("M")[{"X1":1,"X2":0, "Y":0}] = [m/(m+n+o), n/(m+n+o), o/(m+n+o)]
mod.cpt("M")[{"X1":1,"X2":0, "Y":1}] = [p/(p+q+r), q/(p+q+r), r/(p+q+r)]
mod.cpt("M")[{"X1":1,"X2":1, "Y":0}] = [s/(s+t+u), t/(s+t+u), u/(s+t+u)]
mod.cpt("M")[{"X1":1,"X2":1, "Y":1}] = [v/(v+w+x), w/(v+w+x), x/(v+w+x)]

methodsDict = {"PS": {"fit": ccMethodsM6M7.fitPatternSubmodelsM6M7, "predict": ccMethodsM6M7.predictPatternSubmodelsM6M7},
               "CCS": {"fit": ccMethodsM6M7.fitCCS, "predict": ccMethodsM6M7.predictCCS},
               "MIA": {"fit": ccMethodsM6M7.fitMIA, "predict": ccMethodsM6M7.predictMIA},
               "MARG": {"fit": ccMethodsM6M7.fitMarginalisation, "predict": ccMethodsM6M7.predictMarginalisation},
               "MARGMI": {"fit": ccMethodsM6M7.fitMarginalisationMI, "predict": ccMethodsM6M7.predictMarginalisationMI}}

dfList = []

nSim = 100
N = 100000

modelsDict = {"M7":mod}

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
save_dir = os.path.join(BASE_DIR, "outputs", "dxiscrete", "datasets")
save_path = os.path.join(save_dir, "saveDataDiscreteM7_new.csv")
dfSim.to_csv(save_path, index=False)