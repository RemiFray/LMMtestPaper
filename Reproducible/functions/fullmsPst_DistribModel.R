

# ------- Dummy distribution ------ ----


dLMultiStates <- nimbleFunction(
  run = function(x = double(1), 
                 capture = double(2),
                 identification = double(0),
                 delta = double(1),
                 Gamma = double(2),
                 B=double(3),
                 latentObservation = double(2),
                 latentIndex=double(1),
                 nbLatentObs=double(0),
                 S = double(0), nbStades = double(0),
                 log = integer(0, default = 1) 
  ) {
    
    N <- sum(x[1:nbLatentObs])
    
    logProbData <- lfactorial(N)
    
    indexs <- which(x > 0)
    
    for(i in 1:length(indexs)){
      h <- indexs[i]
      logProbData <- logProbData - lfactorial(x[h])
      
      hist <- latentObservation[h,]
      
      hist2 <- hist+1
      hist2[which(hist2 > 3)] <- hist2[which(hist2 > 3)] - 3
      
      alphas <- nimMatrix(nrow = nbStades, ncol = S)
      for(s in 1:nbStades){
        alphas[s, 1] <- delta[s] * B[hist2[1], s, 1]
      }
      for(t in 2:S){
        for(s in 1:nbStades){
          p <- 0
          for(r in 1:nbStades){
            p <- p + alphas[r, t-1] * Gamma[r, s]
          }
          alphas[s, t] <-  p * B[hist2[t], s, t]
        }
      }
      
      probThisObs <- sum(alphas[,S]) * 
        identification^(sum(hist == 1) + sum(hist == 2)) *
        (1-identification)^sum(hist > 3)
      
      logProbData <- logProbData + x[h]*log(probThisObs)
    }
    
    if(log) return(logProbData)
    return(exp(logProbData))
    returnType(double(0))
  })

rLMultiStates <- nimbleFunction(
  run = function(n = integer(), capture = double(2), 
                 identification = double(0), delta = double(1),
                 Gamma = double(2), B=double(3),
                 latentObservation = double(2),
                 latentIndex = double(1), 
                 nbLatentObs = double(0),
                 S = double(0), nbStades = double(0)) {
    returnType(double(1))
    return(delta) ## dummy behavior
  })

deregisterDistributions("dLMultiStates")
registerDistributions(
  list(dLMultiStates = list(BUGSdist = paste0("dLMultiStates(capture, identification,",
                                              "delta, Gamma, B, latentObservation, ",
                                              "latentIndex, nbLatentObs,", 
                                              "S, nbStades)"),
                            types =c("value = double(1)", "capture = double(2)", 
                                     "delta = double(1)", "Gamma = double(2)",
                                     "B=double(3)",
                                     "latentObservation = double(2)", 
                                     "latentIndex = double(1)"))))

# ------- Model Code ------ ----

LMstadesCode <- nimbleCode({
  
  # priors
  for(stade in 1:3){
    for(t in 1:S){
      capture[stade, t] ~ dbeta(1.0, 1.0)
    }
  }
  
  for(t in 1:S){
    B[1, 1:nbStades, t] <- 1 - capture[,t]
    B[2:nbEvents, 1:nbStades, t] <- diag(capture[,t])
  }
  
  paramDirich2[1:2] <- c(1,1)
  paramDirich3[1:3] <- c(1,1,1)
  
  delta[1:3] ~ ddirch(paramDirich3[1:3])
  for(st1 in 1:3){
    Gamma[st1, 1:3] ~ ddirch(paramDirich3[1:3])
  }
  
  alpha ~ dbeta(1.0, 1.0)
  
  N ~ dunif(1, 10e4)
  D <- x[1]
  
  x[1:nbLatentObs] ~ dLMultiStates(capture[1:3, 1:S], alpha, 
                                   delta[1:3], Gamma[1:3,1:3], 
                                   B[1:nbEvents, 1:nbStades, 1:S],
                                   latentObservation[1:nbLatentObs, 1:S],  
                                   latentIndex[1:nbLatentObs],
                                   nbLatentObs, S, nbStades)
  
})

