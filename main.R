source("R/01_data.R")

# TO REMOVE
theta = list(X1 = list(mu = 0, sigma = 1),
             X2 = list(mu = 0, sigma = 1),
             Y = list(beta = c(0,1,1), sigma = .1))
phi = c(0,1,-1,1)
df = SimulateDataContinuous(1000,theta,phi,missCoef = 0,
                              seed = 42, n_mc = 1e5)