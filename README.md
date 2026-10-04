# GMM Misspecification: Dynamic model

This repository contains the code for the dynamic model presented in the research paper:
**"The effect of misspecification in Gaussian mixture models"** (Kate Roberts, 2026).

## Overview
This script ( `dynamic_model.R` ) generates normal and skew-normal data (with labels). A Monte Carlo simulation is then used to fit the GMM on the data. For the normal data an incorrect number of clusters will also be chosen and the GMM will be fit. All of the parameter estimates and model evaluations are saved within the functions.

## Prerequisites

To run the code, ensure that **R** is installed with the following packages:

```R
install.packages(c("MASS", "sn", "mvtnorm", "mixsmsn"))
