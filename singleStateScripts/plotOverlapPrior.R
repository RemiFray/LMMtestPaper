
library(MCMCvis)
library(tidyverse)
library(ggplot2)
library(gridExtra)

# ----- Parameters ----- ----CC


S <- c(5, 7, 9)

alpha_test <- c(0.8, 0.9, 0.95)
p_test <- c(0.4, 0.3, 0.2, 0.1)
N_test <- c(500, 1000)
model <- "LMM_Pt"
prior <- c("noPrior", "priorAlpha100i", "priorAlpha100iBiased")

aPrior <- list(priorAlpha100i = list(a0.8 = c(80, 20),
                                     a0.9 = c(90, 10),
                                     a0.95 = c(95, 5)),
               priorAlpha100iBiased = list(a0.8 = c(74, 26),
                                     a0.9 = c(85, 15),
                                     a0.95 = c(91, 9)))

# empty table of overlaps and empty list of plots
allOvlaps <- as.data.frame(matrix(ncol= 15, nrow=0))


for(pr in prior){
  print(pr)
  for(s in S){
    cat("\t", s, "\n")
    for(N in N_test){
      cat("\t\t", N, "\n")
      for(p in p_test){
        cat("\t\t\t", p, "\n")
        for(a in alpha_test){
          # preparing to load chain files
          params <- paste0("N", N, " p", p, " a", a)
          # plts[[pr]][[paste0("S", s)]][[params]] <- list()
          directory <- paste0("./results/LMM_S", s, "_", pr)
          files <- list.files(directory)
          
          # Prior on alpha
          if(grepl("priorAlpha100i", pr)){
            apr <- paste0("a", a)
            aPr <- aPrior[[pr]][[apr]]
            alphaPrior <- rbeta(1000000, aPr[1], aPr[2])
          } else alphaPrior <- rbeta(1000000, 1, 1)
          
          # data.frames of the overlap
          ovlapTmp <- as.data.frame(matrix(ncol= 15, nrow=1))
          ovlapTmp[1,1] <- pr
          ovlapTmp[1,2:5] <- c(s, N, p, a)
          
          for(iter in 1:10){
            # read chain files
            fileExp <- paste0("N", N, "_a", a, "_p", p, "_iter", iter, ".Rdata")
            filePath <- files[which(grepl(fileExp, files))]
            load(file.path(directory, filePath))
            
            # Some chain hadn't really reached convergence
            if((pr == "noPrior" & N == 1000 & p == 0.1) |
               (pr == "noPrior" & N == 1000 & p == 0.2 & a == 0.8))
              samples <- list(chain1 = samples$chain1[1001:5000,],
                                  chain2 = samples$chain2[1001:5000,])
            else if(pr == "noPrior" & N == 1000 & p == 0.3 & a == 0.8)
              samples <- list(chain1 = samples$chain1[501:5000,],
                                  chain2 = samples$chain2[501:5000,])
            else
              samples <- samples
            
            # calculation of the posterior and overlap
            allAlpha <- c(samples[[1]][,"alpha"], samples[[2]][,"alpha"])
            alphaPost <- density(allAlpha)
            x <- alphaPost$x
            if(grepl("priorAlpha100i", pr)){
              alphaPriorLocal <- dbeta(x, aPr[1], aPr[2])
            } else alphaPriorLocal <- dbeta(x, 1, 1)
            z <- sapply(1:length(alphaPriorLocal), 
                        function(j) min(alphaPost$y[j], alphaPriorLocal[j]))
            ovlapValue <- round(sum((x[2]-x[1])*z),2)
            
            ovlapTmp[, iter+5] <- ovlapValue
          }
          # bind overlap in the global table
          allOvlaps <- allOvlaps %>% rbind(ovlapTmp)
        }
      }
    }
  }
}

# save table with overlaps values
colnames(allOvlaps) <- c("Prior", "S", "N", "p", "a", paste0("iter", 1:10))
allOvlaps <- allOvlaps %>% 
  pivot_longer(starts_with("iter"), 
               names_to = "iter", names_prefix = "iter",
               values_to = "overlap")
write.csv(allOvlaps, file="./figures/CorrAndOverlaps/allOverlaps.csv", row.names = F)
