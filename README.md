# RobustFLR

This repo contains code to run the simulations in the paper Robust Function-on-Function regression.

### Prerequisites

The simulation were built in R version 3.0.2

# Required packages 

library('MASS')
library('fda')
library('fda.usc')
library('rrcov')
library('robustbase')
library('ROCR')

## Authors

* Harjit Hullait - https://github.com/hullait/RobustFLR.git

## Procedure

* Run RobustFLR.R with a set value of alpha and a, to create simulations and output models selected by BIC and RBIC. Also gives the residual error using robust and classical FLR.
* Simulations_FLR.R creates the simulation for given alpha and a
* Model_select_FLR.R calculates BIC and Robust BIC 
* RobustFLR_functions.R contains functions which will be called
