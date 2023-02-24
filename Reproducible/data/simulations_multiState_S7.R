library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----


# Simulation of data
simulCMRStades <- function(N, S, pcapture, identification, 
                           pStade, ptransition){
  
  cr <- matrix(ncol = S, nrow = N, data = 0) # Matrix of capture NbIndiv * NbSample
  
  for(i in 1:N){
    stade <- sample(1:3, 1, prob = pStade)
    for(t in 1:S){
      if(runif(1) < pcapture[stade, t]) { #capture
        if (runif(1) < identification[stade])
          cr[i,t] <- stade # Adding a succesful capture to cr
        else{
          ccr <- cr
          cr <- matrix(nrow = nrow(ccr)+1, ncol = S, data = 0)
          cr[1:nrow(ccr),] <- ccr[,]
          cr[nrow(cr),t] <- stade # add one hist with unique observation
        }
      }
      stade <- sample(1:3, 1, prob = ptransition[stade,])
    }
  }
  
  return(cr)
  
}

# Counts of histories in sorted vector
getSummaryCMR <- function(dta){
  summaryDta <- as.data.frame(table(apply(dta, 1, deparse)), 
                              stringsAsFactors = F)
  summaryDta$index <- sapply(summaryDta$Var1, 
                             function(h) nimGetLatentIndex(eval(parse(text=h))))
  summaryDta <- summaryDta %>% 
    arrange(index) %>% 
    select(history=Var1, index, Freq)
  summaryDta
}

# calculate index of an history
# to be stored in a vector
nimGetLatentIndex <- nimbleFunction(
  run = function(h = double(1)){
    s <- 1 
    for(t in 1:length(h)) s <- s + h[t]*7^(t-1)
    return(s)
    returnType(double(0))
  })

# ----- Parameters ----- ----

S <- 7

simulType <- "fullms_meanTrans"

p_test <- c(0.4, 0.3, 0.2, 0.1)
N_test <- c(500, 1000)
alpha_test <- c(0.8, 0.9, 0.95)

pStade <- c(0.33, 0.4, 0.27)
transition <- t(matrix(ncol = 3, data = c(0.76, 0.12, 0.12,
                                          0.1, 0.8, 0.1,
                                          0.15, 0.15, 0.7)))


bias_plan <- data.frame(N = rep(N_test, each = length(p_test)*length(alpha_test)),
                        a = rep(rep(alpha_test, each = length(p_test)), length(N_test)),
                        p = rep(rep(p_test, length(alpha_test)), length(N_test)))

# ------- Big loop ------ ----


for(line in 1:nrow(bias_plan)){
  set.seed(1234)
  
  capture <- matrix(nrow = 3, ncol = S, data = bias_plan$p[line])
  alpha <- rep(bias_plan$a[line], 3)
  N <- bias_plan$N[line]
  
  for(iter in 1:10){
  
  dta <- simulCMRStades(N, S, capture, alpha, pStade, transition)
  summaryDta <- getSummaryCMR(dta)
  
  if(!any(summaryDta$index == 1)) 
    summaryDta <- rbind(c(deparse(rep(0,S)), 1, 0), summaryDta)
  
  
  save(summaryDta,
       file = paste("./data/", simulType, "_S",S, "/simul_S", S, "_N", N, "_a", alpha[1], 
                    "_p", capture[1], "_iter", iter, ".Rdata", sep = ""))
  }
}
