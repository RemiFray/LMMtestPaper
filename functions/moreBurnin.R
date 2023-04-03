
moreBurnin <- function(model, samples){
  
  if( model == "fullmsPts"){
    # S == 5
    if(s == 5 & pr == "noPrior" & N == 1000){
      if(p == 0.4 & a == 0.8)
        samples <- list(chain1 = samples$chain1[1001:5000,],
                        chain2 = samples$chain2[1001:5000,])
      else if(p == 0.3 & a == 0.8)
        samples <- list(chain1 = samples$chain1[1501:5000,],
                        chain2 = samples$chain2[1501:5000,])
      else if(p == 0.2 & a == 0.8)
        samples <- list(chain1 = samples$chain1[3001:5000,],
                        chain2 = samples$chain2[3001:5000,])
      else if(p == 0.1 & a == 0.8)
        samples <- list(chain1 = samples$chain1[1001:5000,],
                        chain2 = samples$chain2[1001:5000,])
      else if(p == 0.1 & a == 0.9)
        samples <- list(chain1 = samples$chain1[501:5000,],
                        chain2 = samples$chain2[501:5000,])
      # S == 7
    } else if(s == 7 & pr == "noPrior" & N == 1000){
      if(p == 0.3 & a == 0.8)
        samples <- list(chain1 = samples$chain1[501:5000,],
                        chain2 = samples$chain2[501:5000,])
      else if(p == 0.2 & a == 0.8)
        samples <- list(chain1 = samples$chain1[1501:5000,],
                        chain2 = samples$chain2[1501:5000,])
      else if(p == 0.1)
        samples <- list(chain1 = samples$chain1[1001:5000,],
                        chain2 = samples$chain2[1001:5000,])
      # S == 9
    } else if(s == 9 & N == 1000){
      if(pr == "noPrior" & p == 0.2 & a == 0.8)
        samples <- list(chain1 = samples$chain1[1501:5000,],
                        chain2 = samples$chain2[1501:5000,])
      else if(pr == "noPrior" & p == 0.1 & a == 0.8)
        samples <- list(chain1 = samples$chain1[3501:5000,],
                        chain2 = samples$chain2[3501:5000,])
      else if(pr == "noPrior" & p == 0.1 & a == 0.9)
        samples <- list(chain1 = samples$chain1[1501:5000,],
                        chain2 = samples$chain2[1501:5000,])
      else if(pr == "noPrior" & p == 0.1 & a == 0.95)
        samples <- list(chain1 = samples$chain1[1001:5000,],
                        chain2 = samples$chain2[1001:5000,])
      else if(pr == "priorAlpha100i" & p == 0.2 & a == 0.8)
        samples <- list(chain1 = samples$chain1[501:5000,],
                        chain2 = samples$chain2[501:5000,])
      else if(pr == "priorAlpha100i" & p == 0.1 & a == 0.8)
        samples <- list(chain1 = samples$chain1[1001:5000,],
                        chain2 = samples$chain2[1001:5000,])
    }
  } else{
    warning("Model name not found, no burnin was done.")
  }
  
  
  samples
}
