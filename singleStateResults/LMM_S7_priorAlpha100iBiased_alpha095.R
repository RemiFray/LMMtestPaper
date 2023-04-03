library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/LMM_Pt_DistribModel.R")
source("./functions/LMM_Pt_functions.R")
source("./functions/LMM_Pt_samplers.R")

# ----- Parameters ----- ----

S <- 7
nbHlatent <- 3^S

alpha <- 0.95
p_test <- c(0.4, 0.3, 0.2, 0.1)
N_test <- c(500, 1000)

priorType <- "priorAlpha100iBiased"
priorAlpha <- c(91, 9)

bias_plan <- data.frame(N = rep(N_test, each = length(p_test)*10),
                        p = rep(rep(p_test, each = 10), length(N_test)),
                        iter = rep(rep(1:10, length(p_test)), length(N_test)))

burnin <- 30000
nthin <- 200
niter <- 5000*nthin+burnin

# ------- Big loop ------ ----

for(line in 1:nrow(bias_plan)){
  
  N <- bias_plan$N[line]
  
  capture <- rep(bias_plan$p[line], S)
  
  iter <- bias_plan$iter[line]
  
  # ------- Data ------ ----
  
  load(file = paste0("./data/Simul_S",S, "/simul_S", S, "_N", N, "_a", alpha, 
                     "_p", capture[1], "_iter", iter, ".Rdata"))
  
  # MATRIX OF ORDERED OBSERVED HISTORIES
  observation <- as.data.frame(t(sapply(summaryDta$history, 
                                        function(h) eval(parse(text=h)))))
  rownames(observation) <- NULL
  observation$index <- apply(observation, 1, nimGetLatentIndex)
  observation <- observation %>% 
    arrange(index) %>% select(-index) %>% as.matrix()
  colnames(observation) <- NULL
  
  # MATRIX OF HISTORIES AND SPACE FOR OTHERS
  n <- nrow(observation)
  M <- min(sum(summaryDta$Freq[-1]), 3^S-n)
  latentObservation <- matrix(nrow=n+M, ncol = S, data=0)
  latentObservation[1:n,] <- observation
  latentIndex <- numeric(n+M)
  latentIndex[1:n] <- summaryDta$index
  
  xInit <- numeric(n+M)
  xInit[1:n] <- summaryDta$Freq
  
  # MATRIX A (just to check if xInit if valid)
  A <- getA(latentObservation, observation)
  
  # INITIAL X TO START MCMC
  set.seed(1248)
  xInit2 <- xInit
  latObs2 <- latentObservation
  latIdx2 <- latentIndex
  for(j in 1:40){
    tmp <- addError(xInit2, latObs2, latIdx2, S, n, M, n+M)
    xInit2 <- tmp$x
    latObs2 <- tmp$latObs
    latIdx2 <- tmp$latIdx
    
    A2 <- getA(latObs2, observation)
    if(!(all(xInit2 >= 0) & all(summaryDta$Freq[-1] == A2 %*% xInit2)))
      print(j)
  }
  print(paste0("N", N, ", p", capture[1], ", iter", iter))
  cat(all(summaryDta$Freq[-1] == A %*% xInit), "  ")
  A2 <- getA(latObs2, observation)
  cat(all(xInit2 >= 0) & all(summaryDta$Freq[-1] == A2 %*% xInit2), "  ")
  print(sum(xInit-xInit2))
  
  # ------- Constructing nimble model ------ ----
  
  
  LMMPtConsts <- list(S = S, nbLatentObs=n+M)
  
  # Initialisation
  LMMPtInits <- function(i) list(capture = rep(0.5, S),
                                 alpha = 0.5, 
                                 N = c(sum(xInit),sum(xInit2))[[i]],
                                 x = list(xInit, xInit2)[[i]],
                                 latentObservation = list(latentObservation, latObs2)[[i]],
                                 latentIndex = list(latentIndex, latIdx2)[[i]])
  
  suppressMessages(
    LMMPtModel <- nimbleModel(code = LMMPtCode, name = "Pt", 
                              constants = LMMPtConsts,
                              data = list(), inits = LMMPtInits(1))
  )
  
  # Compilation modèle OK
  suppressMessages(
    CLMMPt <- compileNimble(LMMPtModel, showCompilerOutput = FALSE)
  )
  
  # ------- MCMC configuration ------ ----
  
  # Configuration of MCMC
  LMMPtConf <- configureMCMC(LMMPtModel,
                             monitors = c("N", "D", "capture", "alpha"), 
                             thin = 1,
                             # enableWAIC = T,
                             print = F)
  
  LMMPtConf$removeSampler("x", "capture", "alpha", "N")
  LMMPtConf$addSampler(target = c("alpha"), type = nimUpdateA, 
                       control = list(prior = priorAlpha))
  LMMPtConf$addSampler(target = paste0("capture[1:",S,"]"),
                       type = nimUpdatePt)
  LMMPtConf$addSampler(target = c("x"), type = nimSamplerXmove,
                       control = list(D = 1, latentIndex=as.double(latentIndex),
                                      n = n, M = M))
  LMMPtConf$addSampler(target = c("x"), type = nimSamplerX0,
                       control = list(D = 5))
  # LMMPtConf$printSamplers()
  
  # Compilation MCMC
  LMMPtMCMC <- buildMCMC(LMMPtConf)
  
  suppressMessages(
    CLMMPtMCMC <- compileNimble(LMMPtMCMC, project = LMMPtModel,
                                showCompilerOutput = F, resetFunctions = T)
  )
  
  # ------- run MCMC ------ ----
  
  inits <- list(LMMPtInits(1), LMMPtInits(2))
  
  system.time(
    samples <- runMCMC(CLMMPtMCMC, niter = niter, nburnin = burnin, 
                       thin = nthin, nchains = 2,inits = inits, setSeed = 1234,
                       progressBar = F)
  )
  
  save(samples,
       file = paste0("./results/LMM_S",S, "_", priorType, "/LMM_S", S, 
                     "_", priorType, "_N", N, "_a", alpha, 
                     "_p", capture[1], "_iter", iter, ".Rdata"))
  
}



