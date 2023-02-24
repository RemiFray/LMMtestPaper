library(nimble)
library(tidyverse)


# ----- Functions, samplers, distrib and model ----- ----

source("./functions/fullmsPst_functions.R")
source("./functions/fullmsPst_samplers.R")
source("./functions/fullmsPst_DistribModel.R")

# ----- Parameters ----- ----

S <- 9
nbHlatent <- 3^S

alpha <- 0.95
p_test <- c(0.4, 0.3, 0.2)
N_test <- c(500, 1000)

priorType <- "noPrior"
priorAlpha <- c(1, 1)

Gamma <- t(matrix(ncol = 3, data = c(0.76, 0.12, 0.12,
                                          0.1, 0.8, 0.1,
                                          0.15, 0.15, 0.7)))

modelType <- "fullmsPts"
simulType <- "fullms_meanTrans"

bias_plan <- data.frame(N = rep(N_test, each = length(p_test)*10),
                        p = rep(rep(p_test, each = 10), length(N_test)),
                        iter = rep(rep(1:10, length(p_test)), length(N_test)))

burnin <- 60000
nthin <- 100
niter <- 5000*nthin+burnin

# ------- Big loop ------ ----

for(line in 1:nrow(bias_plan)){
  
N <- bias_plan$N[line]

capture <- matrix(nrow = 3, ncol = S, data = bias_plan$p[line])

iter <- bias_plan$iter[line]

# ------- Data ------ ----
  
load(file = paste0("./data/", simulType, "_S",S, "/simul_S", S, "_N", N, "_a", alpha, 
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
m <- sum(summaryDta$Freq[-1])
latentObservation <- matrix(nrow=n+m, ncol = S, data=0)
latentObservation[1:n,] <- observation
latentIndex <- numeric(n+m)
latentIndex[1:n] <- summaryDta$index

xInit <- numeric(n+m)
xInit[1:n] <- summaryDta$Freq

# MATRIX A (just to check if xInit if valid)
A <- getA(latentObservation, observation)

# INITIAL X TO START MCMC
set.seed(1248)
xInit2 <- xInit
latObs2 <- latentObservation
latIdx2 <- latentIndex
for(j in 1:50){
  tmp <- addError(xInit2, latObs2, latIdx2, S, n, m, n+m)
  xInit2 <- tmp$x
  latObs2 <- tmp$latObs
  latIdx2 <- tmp$latIdx
  
  A2 <- getA(latObs2, observation)
  if(!(all(xInit2 >= 0) & all(summaryDta$Freq[-1] == A2 %*% xInit2)))
    print(j)
  # if(j==6) debugonce(addError)
}
cat(all(summaryDta$Freq[-1] == A %*% xInit), " ")
A2 <- getA(latObs2, observation)
cat(all(xInit2 >= 0) & all(summaryDta$Freq[-1] == A2 %*% xInit2), " ")
cat(sum(xInit - xInit2))


# ------- Constructing nimble model ------ ----

fullmsPtsConsts <- list(S = S, nbEvents = 4, nbStades = 3,
                        nbLatentObs=n+m)

# Initialisation
fullmsPtsInits <- function(i) list(capture = matrix(nrow = 3, ncol = S, 
                                                    data = 0.5),
                                   delta = c(0.4, 0.3, 0.3),
                                   Gamma = Gamma,
                                   alpha = 0.5, 
                                   N = c(sum(xInit),sum(xInit2))[[i]],
                                   x = list(xInit, xInit2)[[i]],
                                   latentObservation = list(latentObservation, latObs2)[[i]],
                                   latentIndex = list(latentIndex, latIdx2)[[i]])

fullmsPtsModel <- nimbleModel(code = LMstadesCode, name = "fullmsPts", 
                              constants = fullmsPtsConsts,
                              data = list(), inits = fullmsPtsInits(1), 
                              calculate = T)

# Compilation modèle OK
CfullmsPts <- compileNimble(fullmsPtsModel, showCompilerOutput = FALSE)


# ------- MCMC configuration ------ ----

# Configuration of MCMC
fullmsPtsConf <- configureMCMC(fullmsPtsModel,
                               monitors = c("N", "D", "capture", "alpha",
                                            "delta", "Gamma"), 
                               thin = 1,
                               # enableWAIC = T,
                               print = TRUE)

fullmsPtsConf$removeSampler("x", "capture", "delta", "Gamma", "alpha", "N")
fullmsPtsConf$addSampler(target = c("alpha"), type = nimSamplerAlpha, 
                         control = list(prior = priorAlpha))
fullmsPtsConf$addSampler(target = paste0("capture[1:3,1:",S,"]"),
                         type = nimSamplerCapture)
fullmsPtsConf$addSampler(target = "Gamma[1:3, 1:3]",
                         type = nimSamplerGamma)
fullmsPtsConf$addSampler(target = "delta[1:3]",
                         type = nimSamplerDelta)
fullmsPtsConf$addSampler(target = c("x"), type = nimSamplerXmove,
                         control = list(D = 1, n = n, m = m))
fullmsPtsConf$addSampler(target = c("x"), type = nimSamplerX0,
                         control = list(D = 5))
fullmsPtsConf$printSamplers()

# Compilation MCMC
fullmsPtsMCMC <- buildMCMC(fullmsPtsConf)


CfullmsPtsMCMC <- compileNimble(fullmsPtsMCMC, project = fullmsPtsModel,
                                showCompilerOutput = F, resetFunctions = T)

# ------- run MCMC ------ ----

inits <- list(fullmsPtsInits(1), fullmsPtsInits(2))

system.time(
  samples <- runMCMC(CfullmsPtsMCMC, niter = niter, nburnin = burnin, 
                     thin = nthin, nchains = 2,inits = inits, setSeed = c(1234, 777))
)

save(samples,
     file = paste0("./results/", modelType, "_S",S, "_", priorType, "/", 
                   modelType, "_S", S, "_", priorType, "_N", N, "_a", alpha, 
                   "_p", capture[1], "_iter", iter, ".Rdata"))

}



