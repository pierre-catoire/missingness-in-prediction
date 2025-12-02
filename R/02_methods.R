# 02_methods.R
# Description: Functions to train the models and perform the predictions
# Author: Pierre CATOIRE
# Date: 2025-12-02
# Purpose:
# - train each model on a training set
# - perform the prediction on a testing set

# Preamble
library(tidyverse)

# Functions
fitPatternSubmodels = function(dTrain){
  modelPS = list()
  modelPS[["M0"]] = lm()
}