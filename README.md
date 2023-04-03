
# Table of contents

1. [Introduction](#introduction)
2. [Reproducing results of the paper](#reproduce)


# Introduction <a name="introduction"></a>

This folder contains all the necessary material to reproduce results presented in the paper "Population size estimation with capture-recapture in presence of individual misidentification and low recapture" (Fraysse *et al.* 2023).


# Reproducing results of the paper  <a name="reproduce"></a>

In order to reproduce the results, one just has to run the scripts in the data folder to get the simulations. Then the scripts in the results folder to get the results. 

* data: it contains several script that were used to simulate the data.
* functions: contains all the functions, nimble code and distributions used to run the model. 
* singleStateResults: contains all the main scripts used to run the model on single state simulations, one file to extract the estimates from all the saved markov chains in a single table and one file to calculate the overlaps between priors and posteriors of the parameter $\alpha$.
* multiStateResults: contains all the main scripts used to run the model on multi state simulations and an additional file that extract the estimates from all the saved markov chains in a single table. Some files for S=7 are duplicated and renamed with an ending "_p02", they are only to deal with simulations for which $p_{capture} = 0.2$ with more iterations than the other scripts.
* paperPlots: scripts used to generate the figures from the paper.




