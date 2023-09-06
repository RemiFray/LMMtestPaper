library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/LMM_Pt_functions.R")

# ----- Parameters ----- ----

S <- 7

p_test <- c(0.4, 0.3, 0.2, 0.1)
N_test <- c(500, 1000)
alpha_test <- c(0.8, 0.9, 0.95)

bias_plan <- data.frame(N = rep(N_test, each = length(p_test)*length(alpha_test)),
                        a = rep(rep(alpha_test, each = length(p_test)), length(N_test)),
                        p = rep(rep(p_test, length(alpha_test)), length(N_test)))


# ------- Data folder ------- ----

data_dir <- paste0("data/Simul_S", S)
if(!file.exists(data_dir)) dir.create(data_dir)


# ------- Big loop ------ ----

for(line in 1:nrow(bias_plan)){
  set.seed(1234)
  
  capture <- rep(bias_plan$p[line], S)
  alpha <- bias_plan$a[line]
  N <- bias_plan$N[line]
  
  for(iter in 1:10){
    
    # ------- Preparing data ------ ----
    
    dta <- simulCMR(N, S, capture, alpha)
    summaryDta <- getSummaryCMR(dta)
    
    if(!any(summaryDta$index == 1)) 
      summaryDta <- rbind(c(deparse(rep(0,S)), 1, 0), summaryDta)
    
    save(summaryDta,
         file = paste("./data/Simul_S",S, "/simul_S", S, "_N", N, "_a", alpha, 
                      "_p", capture[1], "_iter", iter, ".Rdata", sep = ""))
  }
}

