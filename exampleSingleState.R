library(nimble)
library(tidyverse)
library(MCMCvis)

# Functions, samplers, distrib and model ----

source("./functions/LMM_Pt_DistribModel.R")
source("./functions/LMM_Pt_functions.R")
source("./functions/LMM_Pt_samplers.R")

# Parameters ----

S <- 9
N <- 500
alpha <- 0.9
set.seed(1235)
capture <- runif(S, min = 0.2, max = 0.8)


# Simulation ----

dta <- simulCMR(N, S, capture, alpha)


# Estimate parameters with LMM ----

samples <- nimEstimateLMM(dta)

MCMCsummary(samples , round = 3)
MCMCtrace(samples, pdf = F, params = c("D", "N", "alpha"))
MCMCtrace(samples, pdf = F, params = c("capture"))



