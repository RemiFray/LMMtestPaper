
# gives a vector of 0 and 1 indicating which are the 
#    observation generated for a given latent history 
latentObsGenerates <- function(histo, vectHists){
  index2 <- which(histo == 2)
  hist2 <- histo
  hist2[index2] <- 0
  
  a <- integer(length(vectHists))
  
  if(any(hist2 > 0)){
    a[which(vectHists == deparse(hist2))] <- 1
  }
  
  if(length(index2) > 0){
    for(i in index2){
      histErr <- numeric(length(histo))
      histErr[i] <- 1
      a[which(vectHists == deparse(histErr))] <- 1
    }
  }
  
  a
}

# get matrix A such as y = Ax
getA <- function(latentObservation, observation){
  vectHists <- apply(observation, 1, deparse)
  A <- t(apply(latentObservation, 1, latentObsGenerates, 
               vectHists = vectHists))
  t(A)[-which(vectHists==deparse(numeric(S))),]
}

# Simulation of data
simulCMR <- function(N, S, pcapture, identification){
  
  cr <- matrix(ncol = S, nrow = N, data = 0) # Matrix of capture NbIndiv * NbSample
  
  for(i in 1:N){
    for(t in 1:S){
      if(runif(1) < pcapture[t]) { #capture
        if (runif(1) < identification)
          cr[i,t] <- 1 # Adding a succesful capture to cr
        else{
          ccr <- cr
          cr <- matrix(nrow = nrow(ccr)+1, ncol = S, data = 0)
          cr[1:nrow(ccr),] <- ccr[,]
          cr[nrow(cr),t] <- 1 # add one hist with unique observation
        }
      }
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

# get intial state for the MCMC (an x without errors)
getInit <- function(dta){
  summaryDta <- getSummaryCMR(dta)
  
  if(!any(summaryDta$index == 1)) 
    summaryDta <- rbind(c(deparse(rep(0,S)), 1, 0), summaryDta)
  
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
  
  # VECTOR WITH INDEX OF THE HISTORIES IN LATENTOBSERVATION
  latentIndex <- numeric(n+M)
  latentIndex[1:n] <- summaryDta$index
  
  # XINIT
  xInit <- numeric(n+M)
  xInit[1:n] <- summaryDta$Freq
  
  list(dta = dta,
       observation = observation, 
       latentObservation = latentObservation,
       latentIndex = latentIndex,
       xInit = xInit)
}

# Adds error to generate different xInit
addError <- function(x, latentObservation, latentIndex, S, n, M, nM){
  
  c <- 1
  xProp <- x
  
  latentObsProp <- latentObservation
  latentIdxProp <- latentIndex
  
  nu0 <- nimNumeric(S, value = 0.0)
  nu2 <- nimNumeric(S, value = 0.0)
  
  X0 <- xProp > 0
  # def nu0
  which_nu0 <- nimNumeric()
  which_nu0 <- which(X0)
  nu0i <- nimSample(which_nu0)
  nu0 <- latentObsProp[nu0i, ]
  idxt <- nimNumeric() ; idxt <- which(nu0 == 0)
  if(length(idxt) > 0){
    t <- nimSample(idxt)
    # def nu1 and nu2
    nu1iIdx <- which(latentIdxProp == 1+3^(t-1))
    if(length(nu1iIdx) > 0){
      nu1i <- nu1iIdx[1]
      nu2 <- nu0 ; nu2[t] <- 2
      idxNu2 <- nimGetLatentIndex(nu2)
      # nu2 already exists
      if(any(latentIdxProp==idxNu2))
        nu2i <- which(latentIdxProp==idxNu2)[1]
      # nu2 needs to be added
      else {
        nu2i <- which(xProp[(n+1):nM] == 0)[1] + n # first free spot
        latentObsProp[nu2i, ] <- nu2
        latentIdxProp[nu2i] <- idxNu2
      }
      b <- numeric(nM)
      b[nu0i] <- - c
      b[nu1i] <- - c
      b[nu2i] <- + c
      xProp <- xProp + b
    }
  }
  
  
  if(all(xProp >= 0)) 
    return(list(x=xProp, latObs=latentObsProp, latIdx=latentIdxProp))
  else return(list(x=x, latObs=latentObservation, latIdx=latentIndex))
}

# Create alternative initialization with random errors added
getAlternatInit <- function(init, nbErr, seed = 1234){
  
  A <- getA(init$latentObservation, init$observation)
  summaryDta <- getSummaryCMR(dta)
  
  if(any(init$xInit < 0)) 
    stop("There is a negative value in the xInit.")
  if(!all(summaryDta$Freq[-1] == A %*% init$xInit)) 
    stop("The initialization given does not validate y=Ax.")
  
  S <- ncol(init$observation)
  maxErr <- sum(init$xInit[which(init$latentIndex %in% 
                                   sapply(1:S, function(t) 1+3^(t-1)))])
  if(nbErr > maxErr)
    warning(paste0("It is not possible there are ", nbErr,
                   " errors given the data."))
  
  xInit2 <- init$xInit
  latObs2 <- init$latentObservation
  latIdx2 <- init$latentIndex
  n <- nrow(init$observation)
  M <- nrow(init$latentObservation) - n
  
  set.seed(seed)
  for(j in 1:nbErr){
    tmp <- addError(xInit2, latObs2, latIdx2, S, n, M, n+M)
    xInit2 <- tmp$x
    latObs2 <- tmp$latObs
    latIdx2 <- tmp$latIdx
    
    A2 <- getA(latObs2, init$observation)
    if(!(all(xInit2 >= 0) & all(summaryDta$Freq[-1] == A2 %*% xInit2)))
      stop(paste0("Something went wrong when adding the error n°", j, 
                  ". Try enter the the function in debug mode to see what happened."))
  }
  
  message(paste0(sum(init$xInit-xInit2),
                 " ghosts have been added to the second initialization."))
  
  list(dta = init$dta,
       observation = init$observation, 
       latentObservation = latObs2,
       latentIndex = latIdx2,
       xInit = xInit2)
}

# nimble sample function (only one element)
nimSample <-  nimbleFunction(
  run = function(x = double(1)){
    l <- length(x)
    ind <- rcat(1, nimRep(1/l, l))
    nsamp <- x[ind]
    
    returnType(double(0))
    return(nsamp)
  })
# CnimSample <- compileNimble(nimSample)

# get the index of an history
nimGetLatentIndex <- nimbleFunction(
  run = function(h = double(1)){
    s <- 1 
    for(t in 1:length(h)) s <- s + h[t]*3^(t-1)
    return(s)
    returnType(double(0))
  })
# CnimGetLatentIndex <- compileNimble(nimGetLatentIndex)



nimEstimateLMM <- function(dta, nbErrSecInit = 50, priorAlpha=c(1,1),
                           niter=130000, burnin=30000, nthin=1,
                           monitors = c("N", "D", "capture", "alpha"),
                           seed = 1234){
  
  seeds <- sample(1:10000, 2)
  init1 <- getInit(dta)
  init2 <- getAlternatInit(init1, nbErrSecInit, seeds[1])
  
  n <- nrow(init1$observation)
  M <- nrow(init1$latentObservation) -  n
  
  # Constants
  LMMPtConsts <- list(S = S, nbLatentObs=n+M)
  
  # Initialisation
  LMMPtInits <- function(i) list(capture = rep(0.5, S),
                                 alpha = 0.5, 
                                 N = c(sum(init1$xInit),sum(init2$xInit))[[i]],
                                 x = list(init1$xInit, init2$xInit)[[i]],
                                 latentObservation = list(init1$latentObservation,
                                                          init2$latentObservation)[[i]],
                                 latentIndex = list(init1$latentIndex, 
                                                    init2$latentIndex)[[i]])
  
  
  LMMPtModel <- suppressMessages(
                  nimbleModel(code = LMMPtCode, name = "LMMPtModel", 
                              constants = LMMPtConsts,
                              data = list(), inits = LMMPtInits(1))
                            )
  
  message("Compiling model...")
  # Compilation modèle
  CLMMPt <- suppressMessages(
              compileNimble(LMMPtModel, showCompilerOutput = FALSE)
              )
  
  
  # Configuration of MCMC
  LMMPtConf <- configureMCMC(LMMPtModel,
                             monitors = c("N", "D", "capture", "alpha"), 
                             thin = 1, print = FALSE)
  
  LMMPtConf$removeSampler("x", "capture", "alpha", "N")
  LMMPtConf$addSampler(target = c("alpha"), type = nimUpdateA, 
                       control = list(prior = priorAlpha))
  LMMPtConf$addSampler(target = paste0("capture[1:",S,"]"),
                       type = nimUpdatePt)
  LMMPtConf$addSampler(target = c("x"), type = nimSamplerXmove,
                       control = list(D = 1, n = n, M = M))
  LMMPtConf$addSampler(target = c("x"), type = nimSamplerX0,
                       control = list(D = 5))
  
  # Compilation MCMC
  LMMPtMCMC <- buildMCMC(LMMPtConf)
  
  message("Compiling MCMC...")
  
  CLMMPtMCMC <- suppressMessages(
                  compileNimble(LMMPtMCMC, project = LMMPtModel,
                                showCompilerOutput = F, resetFunctions = T)
                  )
  
  
  # ------- run MCMC ------ ----
  
  inits <- list(LMMPtInits(1), LMMPtInits(2))
  
  runMCMC(CLMMPtMCMC, 
          niter = niter, nburnin = burnin, thin = nthin, 
          nchains = 2, inits = inits, setSeed = seeds)

}



