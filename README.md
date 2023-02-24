
# Table of contents

1. [Introduction](#introduction)
2. [Simulating data](#simulation)
3. [Running the LMM](#runLMM)
3. [Reproducing results of the paper](#reproduce)


# Introduction <a name="introduction"></a>

This folder contains all the necessary material to reproduce results presented in the paper "Population size estimation with capture-recapture in presence of individual misidentification and low recapture" (Fraysse *et al.* 2023).
There is also all base material needed to run the latent multinomial model in closed population with single or multi-state with nimble. 

Running the LMM might might need to modify some functions for all cases different from what the paper covers (*i.e.* for anything other than closed population with single or three-state). Basic instruction to do so are in section [Running the LMM](#runLMM). If you encounter a problem or find that the instructions are not sufficient or clear enough, please ask/report at ???

# Simulating data <a name="simulation"></a>

In order to simulate CMR data, use either function ??? for single state or function ??? for multistate. They take the following parameters:

* N: Population size
* S: number of capture occasions
* p: Capture vector or matrix with occasions in cols and states in rows
* $\psi$: Transition probability matrix

In order to simulation misidentification in a CMR dataset, use function ???. It takes the following parameters:

* dta: a matrix with CMR data, individuals in rows, capture occasions in cols
* $\alpha$: the identification probability


# Running the LMM  <a name="runLMM"></a>




# Reproducing results of the paper  <a name="reproduce"></a>


