
# Table of contents

1. [Introduction](#introduction)
2. [Reproducing results of the paper](#reproduce)
3. [Simulating data](#simulation)
4. [Running the LMM](#runLMM)


# Introduction <a name="introduction"></a>

This folder contains all the necessary material to reproduce results presented in the paper "Population size estimation with capture-recapture in presence of individual misidentification and low recapture" (Fraysse *et al.* 2023).
There is also all base material needed to run the latent multinomial model in closed population with single or multi-state with nimble. 

Running the LMM might might need to modify some functions for all cases different from what the paper covers (*i.e.* for anything other than closed population with single or three-state). Basic instruction to do so are in section [Running the LMM](#runLMM). If you encounter a problem or find that the instructions are not sufficient or clear enough, please ask/report at ???


# Reproducing results of the paper  <a name="reproduce"></a>

All the script that generated the results are in the "Reproducible" folder. In order to reproduce it, one just has to run all of them. There are several folder inside:

* data: it contains several script that were used to simulate the data.
* functions: contains all the functions, nimble code and distributions used to run the model. 
* singleStateResults: contains all the main scripts used to run the model on single state simulations, one file to extract the estimates from all the saved markov chains in a single table and one file to calculate the overlaps between priors and posteriors of the parameter $\alpha$.
* multiStateResults: contains all the main scripts used to run the model on multi state simulations and an additional file that extract the estimates from all the saved markov chains in a single table. Some files for S=7 are duplicated and renamed with an ending "_p02", they are only to deal with simulations for which $p_{capture} = 0.2$ with more iterations than the other scripts.
* paperPlots: scripts used to generate the figures from the paper.


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



