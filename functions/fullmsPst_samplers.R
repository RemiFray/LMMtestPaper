
# Sampler error of identifications (x)
nimSamplerXmove <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    D <- if(!is.null(control$D))         control$D         else 1
    S <- dim(model$latentObservation)[2]
    n <- control$n
    m <- control$m
    nm <- n+m
    nbStades <- dim(model$B)[2]
    forward <- forwardGen(model, nbStades)
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    
    c <- nimSample(c(1:D))
    xProp <- model$x
    
    latentObsProp <- model$latentObservation
    latentIdxProp <- model$latentIndex
    
    log_qratio <- 0
    skip<-TRUE
    
    nu0 <- nimNumeric(S, value = 0.0)
    nu2 <- nimNumeric(S, value = 0.0)
    
    choice <- runif(1, 0, 1)
    u <- nimSample(1:3)
    
    ## ----------------------------------------------- add error
    if(choice <= 0.5){
      
      # assuming its better to check if there is a 0 in the 
      #   one history sampled than sample one in which we know there is one
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
        nu1iIdx <- which(latentIdxProp == 1+u*7^(t-1))
        if(length(nu1iIdx) > 0){
          nu1i <- nu1iIdx[1]
          nu2 <- nu0 ; nu2[t] <- u+3
          idxNu2 <- nimGetLatentIndex(nu2)
          # nu2 already exists
          if(any(latentIdxProp==idxNu2))
            nu2i <- which(latentIdxProp==idxNu2)[1]
          # nu2 needs to be added
          else {
            nu2i <- which(xProp[(n+1):nm] == 0)[1] + n # first free spot
            latentObsProp[nu2i, ] <- nu2
            latentIdxProp[nu2i] <- idxNu2
          }
          b <- numeric(nm)
          b[nu0i] <- - c
          b[nu1i] <- - c
          b[nu2i] <- + c
          xProp <- xProp + b
          # # calculate q(x|xprop) / x(xprop|x)
          Xu <- getXu(xProp[(n+1):nm], u+3, latentObsProp)
          qx <- 1/(sum(Xu)*sum(nu2 == (u+3)))
          qxprop <- 1/(sum(X0)*sum(nu0 == 0))
          log_qratio <- log(qx / qxprop)
          
          skip <- FALSE
        }
      }
    }
    
    # ----------------------------------------------- remove error 
    else {
      
      Xu <- getXu(xProp[(n+1):nm], u+3, latentObsProp)
      if(sum(Xu) != 0){
        # def nu2
        which_nu2 <- nimNumeric()
        which_nu2 <- which(Xu == 1)
        nu2i <- nimSample(which_nu2)+n
        nu2 <- latentObsProp[nu2i, ]
        # def t
        idxt <- nimNumeric() ; idxt <- which(nu2 == u+3)
        t <- nimSample(idxt)
        # def nu1 and nu2
        nu1i <- which(latentIdxProp == 1+u*7^(t-1))[1]
        nu0 <- nu2 ; nu0[t] <- 0
        idxNu0i <- nimGetLatentIndex(nu0)
        if(any(latentIdxProp == idxNu0i))
          nu0i <- which(latentIdxProp == idxNu0i)[1]
        else{
          nu0i <- which(xProp[(n+1):nm] == 0)[1] + n # first free spot
          latentObsProp[nu0i, ] <- nu0
          latentIdxProp[nu0i] <- idxNu0i
        }
        # def xProp
        b <- numeric(nm)
        b[nu0i] <- + c
        b[nu1i] <- + c
        b[nu2i] <- - c
        xProp <- xProp + b
        if(xProp[nu2i] == 0)
          latentIdxProp[nu2i] <- 0
        # calculate q(x|xprop) / x(xprop|x)
        X0 <- xProp > 0
        qx <- 0.25/(sum(X0)*sum(nu0 == 0))
        qxprop <- 0.25/(sum(Xu)*sum(nu2 == (u+3)))
        log_qratio <- log(qx / qxprop)
        skip <- FALSE
      }
    }
    ## ------------------------------ metrolis ratio and acceptance
    if(all(xProp >= 0) & !skip){
      
      idxNu <- c(nu0i, nu1i, nu2i)
      model_lprop <- logRatio(idxNu, b, xProp, latentObsProp)
      
      log_MH_ratio <- model_lprop + log_qratio
      
      u <- runif(1, 0, 1)
      if(u < exp(log_MH_ratio)){
        model[[target]] <<- xProp
        model[["N"]] <<- sum(xProp)
        model[["latentObservation"]] <<- latentObsProp
        model[["latentIndex"]] <<- latentIdxProp
        
        copy(from = model, to = mvSaved, row = 1,
             nodes = c(calcNodes, "N"), logProb = T)
      }
    }
  },
  methods = list(
    getXu = function(y2 = double(1), u = double(0), 
                     latentObservation=double(2)){
      Xu <- nimNumeric(length = m)
      indexs <- which(y2 > 0)
      if(length(indexs) > 0){
        for(i in 1:length(indexs)){
          I <- indexs[i]+n
          hist <- latentObservation[I,]
          if(any(hist == u)) Xu[indexs[i]] <- 1
        }
      }
      return(Xu)
      returnType(double(1))},
    logRatio = function (idxNu=double(1), b=double(1), xProp = double(1), 
                         histories=double(2)) {
      N <- sum(model$x)
      Nprop <- sum(xProp)
      
      logRes <- lfactorial(Nprop)-lfactorial(N)
      for(i in 1:3){
        I <- idxNu[i]
        hist <- histories[I,]
        
        als <- forward$run(hist, S)
        prob <- sum(als[,S]) *
          model$alpha[1]^( sum(hist == 1)+sum(hist == 2)+sum(hist == 3)) *
          (1-model$alpha[1])^sum(hist > 3)
        
        logRes <- logRes + lfactorial(model$x[I]) - lfactorial(xProp[I]) +
          b[I] * log(prob)
      }
      return(logRes)
      returnType(double(0))
    },
    reset = function () {})
)

# Sampler unseen individuals (x[1])
nimSamplerX0 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    D <- if(!is.null(control$D))         control$D         else 10
    S <- dim(model$latentObservation)[2]
    nbStades <- dim(model$B)[2]
    forward <- forwardGen(model, nbStades)
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    
    x0move <- nimSample(-D:D)
    xProp <- model$x
    if(xProp[1] + x0move >= 0){
      
      log_MH_ratio <- logRatio(x0move)
      
      u <- runif(1, 0, 1)
      if(u < exp(log_MH_ratio)){
        xProp[1] <- xProp[1] + x0move
        model[[target]] <<- xProp
        model[["D"]] <<- xProp[1]
        model[["N"]] <<- sum(xProp)
        
        copy(from = model, to = mvSaved, row = 1,
             nodes = c(calcNodes, "N", "D"), logProb = T)
      }
    }
  },
  methods = list(
    logRatio = function (x0move=double(0)) {
      N <- sum(model$x)
      Nprop <- N+x0move
      
      logRes <- lfactorial(Nprop)-lfactorial(N)
      
      hist <- nimNumeric(S)
      als <- forward$run(hist, S)
      logRes <- logRes + lfactorial(model$x[1]) - lfactorial(model$x[1]+x0move) +
        x0move * log(sum(als[,S]))
      
      return(logRes)
      returnType(double(0))
    },
    reset = function () {} )
)


# Sampler proba capture (p_capture[stades, time])
nimSamplerCapture <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    S <- dim(model$latentObservation)[2]
    nbStades <- dim(model$B)[2]
    forward <- forwardGen(model, nbStades)
    backward <- backwardGen(model, nbStades)
  },
  run = function() {
    
    ast <- nimMatrix(nrow = 3, ncol = S, value = 1)
    bst <- nimMatrix(nrow = 3, ncol = S, value = 1)
    
    indexs <- which(model$x > 0)
    
    for(i in 1:length(indexs)){
      I <- indexs[i]
      hist <- model$latentObservation[I,]
      hist2 <- hist
      hist2[hist > 3] <- hist2[hist > 3] - 3
      
      als <- forward$run(hist, S)
      bes <- backward$run(hist, S)
      
      posts <- als * bes / sum(als[,S])
      
      for(t in 1:S){
        if(hist[t] > 0)
          ast[hist2[t], t] <- ast[hist2[t], t] + model$x[I]
        else{
          postStades <- rmulti(1, model$x[I], posts[, t])
          bst[, t] <- bst[, t] + postStades
        }
      } # end for t
    } # end for i
    
    new_p <- nimMatrix(nrow = 3, ncol = S)
    for (s in 1:3){
      for(t in 1:S)
        new_p[s, t] <- rbeta(n=1, shape1=ast[s, t], shape2=bst[s, t])
    }
    model[[target]] <<- new_p
    model$simulate("B")
    copy(from = model, to = mvSaved, row = 1,
         nodes = c(target, "B"), logProb = T)
  },
  methods = list( reset = function () {} )
)

# Sampler proba transition (Gamma)
nimSamplerGamma <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    S <- dim(model$latentObservation)[2]
    nbStades <- dim(model$B)[2]
    forward <- forwardGen(model, nbStades)
    backward <- backwardGen(model, nbStades)
  },
  run = function() {
    
    alphaDirich <- nimMatrix(nrow = 3, ncol = 3, value = 1)
    
    indexs <- which(model$x > 0)
    
    for(i in 1:length(indexs)){
      I <- indexs[i]
      hist <- model$latentObservation[I,]
      hist2 <- hist+1
      hist2[which(hist2>4)] <- hist2[which(hist2>4)] - 3
      
      als<- forward$run(hist, S)
      bes<- backward$run(hist, S)
      
      p <- nimMatrix(ncol = 3, nrow = 3)
      for(t in 1:(S-1)){
        for(s in 1:3){
          for(r in 1:3){
            p[s,r]<-als[s, t]*model$Gamma[s, r]*model$B[hist2[t+1], r, t+1]*bes[r, t+1]
          }
        }
        p <- p / sum(p)
        
        sampT <- rmulti(1, model$x[I], c(p[1,], p[2,], p[3,]))
        alphaDirich[1, ] <- alphaDirich[1, ] + sampT[1:3]
        alphaDirich[2, ] <- alphaDirich[2, ] + sampT[4:6]
        alphaDirich[3, ] <- alphaDirich[3, ] + sampT[7:9]
      }
    } # end for i
    
    new_transition <- nimMatrix(nrow = 3, ncol = 3)
    new_transition[1,] <- rdirch(n=1, alpha=alphaDirich[1,])
    new_transition[2,] <- rdirch(n=1, alpha=alphaDirich[2,])
    new_transition[3,] <- rdirch(n=1, alpha=alphaDirich[3,])
    model[[target]] <<- new_transition
    copy(from = model, to = mvSaved, row = 1,
         nodes = target, logProb = T)
  },
  methods = list(
    reset = function() {}
  )
)

# Sampler proba initial stades (delta)
nimSamplerDelta <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    S <- dim(model$latentObservation)[2]
    nbStades <- dim(model$B)[2]
    forward <- forwardGen(model, nbStades)
    backward <- backwardGen(model, nbStades)
  },
  run = function() {
    
    alphaDirich <- nimNumeric(length = 3, value = 1)
    
    indexs <- which(model$x > 0)
    
    for(i in 1:length(indexs)){
      I <- indexs[i]
      hist <- model$latentObservation[I,]
      
      if(hist[1] > 0){
        if(hist[1]> 3) st1 <- hist[1] - 3
        else st1 <- hist[1]
        alphaDirich[st1] <- alphaDirich[st1] + model$x[I]
      }
      
      else{
        if(any(hist > 0)) hist2 <- hist[1:(which(hist > 0)[1])]
        else hist2 <- hist
        tmax <- length(hist2)
        
        als<- forward$run(hist2, tmax)
        bes<- backward$run(hist2, tmax)
        
        posts <- als[,1] * bes[,1] / sum(als[,tmax])
        
        sampStds <- rmulti(1, model$x[I], posts)
        for(s in 1:nbStades){
          alphaDirich[s] <- alphaDirich[s] + sampStds[s]
        }
      }
      
    } # end for i
    
    new_stade <- rdirch(n=1, alpha=alphaDirich)
    model[[target]] <<- new_stade
    copy(from = model, to = mvSaved, row = 1,
         nodes = target, logProb = T)
  },
  methods = list(
    reset = function() {}
  )
)

# Sampler proba identification (alpha)
nimSamplerAlpha <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    nbLatentObs <- dim(model$latentObservation)[1]
    prior <- control$prior
  },
  run = function() {
    aa = prior[1]
    ba = prior[2]
    
    indexs <- which(model$x > 0)
    
    for (i in 1:length(indexs)){
      I <- indexs[i]
      h <- model$latentObservation[I,]
      aa <- aa + model$x[I]*(sum(h[] == 1) + sum(h[] == 2) + sum(h[] == 3))
      ba <- ba + model$x[I]*(sum(h[] == 4) + sum(h[] == 5) + sum(h[] == 6))
    }
    new_a <- rbeta(1, aa, ba)
    
    model[[target]] <<- new_a
    copy(from = model, to = mvSaved, row = 1,
         nodes = target, logProb = T)
  },
  methods = list( reset = function () {} )
)