library(tidyverse)
library(scales)
library(ggplot2)
library(gridExtra)
library(cowplot)

# Data
biasSummary <- read.csv("./results/summaries/AllResults.csv") %>% 
  filter(N == 500) %>% 
  group_by(p, alpha, S, Prior) %>% 
  summarise(N_estim = mean(N_mean),
            N_2.5 = mean(N_2.5),
            N_97.5 = mean(N_97.5)) %>% 
  ungroup() %>% 
  mutate(p = case_when(alpha == 0.8 ~ p-0.01,
                       alpha == 0.9 ~ p,
                       alpha == 0.95 ~ p+0.01),
         alpha = as.factor(alpha),
         Prior = fct_recode(Prior, "Uninformative" = "noPrior",
                            "Unbiased" = "priorAlpha100i",
                            "Biased" = "priorAlpha100iBiased"))
minN <- floor(min(biasSummary$N_2.5)/10)*10
maxN <- max(500, ceiling(max(biasSummary$N_97.5)/10)*10)

# Grid of plot coordinates
S <- c(5, 7, 9)
priors <- unique(biasSummary$Prior)

# Themes and graphic parameters
pointsSize <- 2
linesWidth <- 0.8
textSize <- 20

myThemePdf <- 
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', 
                                 size=linesWidth, linetype='solid'),
        axis.ticks = element_line(size = linesWidth),
        axis.title = element_blank(),
        strip.text = element_text(size = textSize),
        panel.grid.minor.x = element_line(size = linesWidth),
        panel.grid.minor.y = element_line(size = linesWidth),
        strip.background = element_rect(size = linesWidth),
        legend.position = "none")


# empty list of plots
plts <- list("5" = list(), "7" = list(), "9" = list())
for(s in unique(biasSummary$S)){
  for(pr in unique(priors)){
    
    # subset of the data
    biasSummaryTmp <- biasSummary %>%
      filter(S == s, Prior == pr)
    
    # boxplots of estimations of N 
    plts[[as.character(s)]][[pr]] <- biasSummaryTmp %>% 
      ggplot(aes(x=p, y=N_estim))+
      geom_line(aes(group = alpha)) +
      geom_point(aes(shape = alpha)) +
      geom_errorbar(aes(ymin = N_2.5, ymax = N_97.5), width = 0.01) +
      geom_hline(aes(yintercept = 500), lty = "dashed", col = "grey") +
      scale_shape_manual(values=c(16, 0, 17))+
      ylim(minN, maxN) +
      myThemePdf
    
    
  }
}


#legend
p <- biasSummaryTmp %>% 
  ggplot(aes(x=p, y=N_estim))+
  geom_point(aes(shape = alpha)) +
  scale_shape_manual(values=c(16, 0, 17))+
  theme_bw()
legend <- get_legend(p)

# labels of the prior used
prior1 <- ggdraw() + draw_label(unique(biasSummary$Prior)[1])
prior2 <- ggdraw() + draw_label(unique(biasSummary$Prior)[2])
prior3 <- ggdraw() + draw_label(unique(biasSummary$Prior)[3])
priorGrid <- plot_grid(prior1, prior2, prior3, ncol = 3)

# label of the number of occasion
s1 <- ggdraw() + draw_label("5")
s2 <- ggdraw() + draw_label("7")
s3 <- ggdraw() + draw_label("9")
sGrid <- plot_grid(s1, s2, s3, ncol = 1)

# grids of plots
plotGrid5 <- plot_grid(plotlist = plts[["5"]], ncol = 3)
plotGrid7 <- plot_grid(plotlist = plts[["7"]], ncol = 3)
plotGrid9 <- plot_grid(plotlist = plts[["9"]], ncol = 3)
plotGrid <- plot_grid(plotGrid5, plotGrid7, plotGrid9, ncol = 1)

# grid with labels and boxplots
pGrid <- plot_grid(NULL, priorGrid, NULL,
                   sGrid, plotGrid, legend,
                   rel_widths = c(0.05, 1, 0.1), 
                   rel_heights = c(0.05, 1))
  
  


# Save plots
pdf(file=paste("./figures/", "EstimPlots.pdf", sep=""),
    onefile = TRUE)
plot(pGrid)
dev.off()

