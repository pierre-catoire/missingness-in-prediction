import pyagrum as gum # bayesian networks
import numpy as np
import pandas as pd

# Feed the model with theta, the parameters of Pr(X1, X2, Y)
def setTheta(model, theta):
    model.cpt("X1").fillWith([1-theta["x1"], theta["x1"]])
    model.cpt("X2").fillWith([1-theta["x2"], theta["x2"]])
    model.cpt("Y")[{"X1":0,"X2":0}] = [1-theta["y.x1_0.x2_0"], theta["y.x1_0.x2_0"]]
    model.cpt("Y")[{"X1":0,"X2":1}] = [1-theta["y.x1_0.x2_1"], theta["y.x1_0.x2_1"]]
    model.cpt("Y")[{"X1":1,"X2":0}] = [1-theta["y.x1_1.x2_0"], theta["y.x1_1.x2_0"]]
    model.cpt("Y")[{"X1":1,"X2":1}] = [1-theta["y.x1_1.x2_1"], theta["y.x1_1.x2_1"]]

# Feed the model with phi, the parameters of P(M1|X1,X2,Y)
# missCoef is the target probability of missingness.
# the function tunes the CPT of M1 such as Pr(M1 = 1) approx. missCoef
def setPhi(model, modelLabel, missCoef):

    def set_M1_CPT_with_missCoef(model, modelLabel, missCoef):
        if modelLabel == "M1":
            model.cpt("M1").fillWith([1-missCoef,missCoef])
        elif modelLabel == "M2":
            model.cpt("M1")[{"X2":0}] = [1-(.5*missCoef),.5*missCoef]
            model.cpt("M1")[{"X2":1}] = [1-(1.2*missCoef),1.2*missCoef]
        elif modelLabel == "M3":
            model.cpt("M1")[{"X1":0}] = [1-(.5*missCoef),.5*missCoef]
            model.cpt("M1")[{"X1":1}] = [1-(1.2*missCoef),1.2*missCoef]
        elif modelLabel == "M4":
            model.cpt("M1")[{"X1":0, "Y":0}] = [1-(.5*missCoef),.5*missCoef]
            model.cpt("M1")[{"X1":0, "Y":1}] = [1-(1.2*missCoef),1.2*missCoef]
            model.cpt("M1")[{"X1":1, "Y":0}] = [1-(.2*missCoef),.2*missCoef]
            model.cpt("M1")[{"X1":1, "Y":1}] = [1-missCoef,missCoef]
        elif modelLabel == "M5":
            model.cpt("M1")[{"Y":0}] = [1-(.5*missCoef),.5*missCoef]
            model.cpt("M1")[{"Y":1}] = [1-(1.2*missCoef),1.2*missCoef]

    def PM1_given_coef(model, modelLabel, missCoef):
        set_M1_CPT_with_missCoef(model, modelLabel, missCoef)
        inference = gum.LazyPropagation(model)
        inference.makeInference()
        return inference.posterior("M1")[1]
    
    # Iterative search of the best tune for the CPT of M1
    low = 0
    high = 1
    maxIter = 100
    for _ in range(maxIter):
        mid = 0.5 * (low + high)
        p = PM1_given_coef(model, modelLabel, mid)

        if abs(p - missCoef) < .001:
            return mid

        if p < missCoef:
            low = mid
        else:
            high = mid

    set_M1_CPT_with_missCoef(model, modelLabel, mid)

def predictFromEvidence(model, df, target = "Y"):
    ie = gum.LazyPropagation(model)  # reuse inference engine
    cols = list(df.columns)
    data = df.values              # numpy array for fast access
    n = data.shape[0]
    out = np.zeros(n)

    idxM1 = cols.index("M1") if "M1" in cols else None
    
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