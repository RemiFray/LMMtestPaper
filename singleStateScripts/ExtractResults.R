
library(dplyr)
library(MCMCvis)

wideSummary <- function(tmp){
  if("D" %in% rownames(tmp)) idx <- c(2:4,6,8)
  else idx <- c(1:3,5,7)
  
  rnames <- gsub("\\[|\\]", "", rownames(tmp)[idx], perl=TRUE)
  rnames <- gsub("capture", "p", rnames, perl=TRUE)
  rnames <- gsub("alpha", "a", rnames, perl=TRUE)
  cnames <- gsub("%", "", colnames(tmp), perl=TRUE)
  allCols <- as.character(sapply(rnames, function(x) paste(x, cnames, sep="_")))

  df <- as.matrix(tmp[idx,]) %>% 
    t() %>% as.numeric() %>% 
    matrix(nrow=1) %>% as.data.frame()
  colnames(df) <- allCols
  df
}

S <- c(5, 7, 9)

alpha_test <- c(0.8, 0.9, 0.95)
p_test <- c(0.4, 0.3, 0.2, 0.1)
N_test <- c(500, 1000)
model <- "LMM_Pt"
prior <- c("noPrior", "priorAlpha100i", "priorAlpha100iBiased")

testSummaries <- as.data.frame(matrix(ncol = 5*7, nrow=0))

for(pr in prior){
  print(pr)
  for(s in S){
    cat("\t", s, "\n")
    for(N in N_test){
      cat("\t\t", N, "\n")
      for(p in p_test){
        cat("\t\t\t", p, "\n")
        for(a in alpha_test){
          directory <- paste0("./results/LMM_S", s, "_", pr)
          # if(!file.exists(directory)) cat(directory, "\n")
          files <- list.files(directory)
          for(iter in 1:10){
            fileExp <- paste0("N", N, "_a", a, "_p", p, "_iter", iter, ".Rdata")
            filePath <- files[which(grepl(fileExp, files))]
            load(file.path(directory, filePath))
            
            if((pr == "noPrior" & N == 1000 & p == 0.1) |
               (pr == "noPrior" & N == 1000 & p == 0.2 & a == 0.8))
              reducedSamp <- list(chain1 = samples$chain1[1001:5000,],
                                  chain2 = samples$chain2[1001:5000,])
            else if(pr == "noPrior" & N == 1000 & p == 0.3 & a == 0.8)
              reducedSamp <- list(chain1 = samples$chain1[501:5000,],
                                  chain2 = samples$chain2[501:5000,])
            else
              reducedSamp <- samples
            
            tmp <- data.frame(Prior=pr, S=s, N=N, p=p, 
                       alpha=a, iter=iter) %>% 
              cbind(wideSummary(MCMCsummary(reducedSamp)))
            testSummaries <- rbind(testSummaries, tmp)
            
          }
        }
      }
    }
  }
}

write.csv(testSummaries, file="./results/summaries/AllResults.csv", row.names = F)
  
  
  
