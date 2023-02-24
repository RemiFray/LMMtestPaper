
moreBurnin <- function(model, samples){
  
  if( model == "fullmsPts"){
    if(s == 5 & pr == "noPrior" & N == 1000 & a == 0.8){
      if(p == 0.4)
        samples <- list(chain1 = samples$chain1[1001:5000,],
                        chain2 = samples$chain2[1001:5000,])
      else if(p == 0.3){
        samples <- list(chain1 = samples$chain1[1501:5000,],
                        chain2 = samples$chain2[1501:5000,])
      }
      else if(p == 0.2){
        samples <- list(chain1 = samples$chain1[3001:5000,],
                        chain2 = samples$chain2[3001:5000,])
      }
    } else if(s == 7 & pr == "noPrior" & N == 1000 & a == 0.8){
      if(p == 0.3){
        samples <- list(chain1 = samples$chain1[501:5000,],
                        chain2 = samples$chain2[501:5000,])
      } else if(p == 0.2){
        samples <- list(chain1 = samples$chain1[1501:5000,],
                        chain2 = samples$chain2[1501:5000,])
      }
    } else if(s == 9 & N == 1000 & p == 0.2 & a == 0.8){
      if(pr == "noPrior")
        samples <- list(chain1 = samples$chain1[1501:5000,],
                        chain2 = samples$chain2[1501:5000,])
      else if(pr == "priorAlpha100i")
        samples <- list(chain1 = samples$chain1[501:5000,],
                        chain2 = samples$chain2[501:5000,])
    }
  } else{
    warning("Model name not found, no burnin was done.")
  }
  
  
  samples
}
