
# Data
biasSummary <- read.csv("./data/fullmsPts_AllResults.csv") %>% 
  filter(N == 500) %>% 
  mutate(p = case_when(alpha == 0.8 ~ p-0.01,
                       alpha == 0.9 ~ p,
                       alpha == 0.95 ~ p+0.01),
         alpha = as.factor(alpha),
         Prior = fct_recode(Prior, "Uninformative" = "noPrior",
                            "Unbiased" = "priorAlpha100i"))

biasSummaryMeans <- biasSummary %>% 
  group_by(p, alpha, S, Prior) %>% 
  summarise(N_estim = mean(N_mean),
            N_2.5 = mean(N_2.5),
            N_97.5 = mean(N_97.5)) %>% 
  ungroup()

minN <- min(floor(min(biasSummaryMeans$N_2.5)/10)*10,
            min(biasSummary$N_mean))
maxN <- max(500, 
            ceiling(max(biasSummaryMeans$N_97.5)/10)*10, 
            max(biasSummary$N_mean))


# Grid of plot coordinates
S <- c(5, 7, 9)
priors <- unique(biasSummary$Prior)

# Themes and graphic parameters
pointsSize <- 2
boxLinesWidth <- 0.6
plotLinesWidth <- 0.2
labelSize <- 11

myThemePdf <- 
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', 
                                 size=boxLinesWidth, linetype='solid'),
        axis.ticks = element_line(size = boxLinesWidth),
        axis.title = element_blank(),
        panel.grid.minor.x = element_line(size = boxLinesWidth*0.8),
        panel.grid.minor.y = element_line(size = boxLinesWidth*0.8),
        legend.position = "none")


# empty list of plots
plts <- list("5" = list(), "7" = list(), "9" = list())
for(s in unique(biasSummary$S)){
  for(pr in unique(priors)){
    
    # subset of the data
    biasSummaryTmp <- biasSummary %>%
      filter(S == s, Prior == pr)
    biasSummaryMeansTmp <- biasSummaryMeans %>% 
      filter(S == s, Prior == pr)
    
    # boxplots of estimations of N 
    plts[[as.character(s)]][[pr]] <- ggplot()+
      geom_hline(aes(yintercept = 500), lty = "dashed", col = "grey") +
      # all points
      geom_point(data = biasSummaryTmp,
                 aes(x=p, y=N_mean, shape = alpha),
                 col="grey") +
      # means
      geom_line(data = biasSummaryMeansTmp,
                aes(x=p, y=N_estim, group = alpha), 
                size = plotLinesWidth) +
      geom_point(data = biasSummaryMeansTmp,
                 aes(x=p, y=N_estim, shape = alpha)) +
      # means of 95% interval limits
      geom_point(data = biasSummaryMeansTmp,
                 aes(x=p, y=N_2.5), pch = 25) +
      geom_point(data = biasSummaryMeansTmp,
                 aes(x=p, y=N_97.5), pch = 24) +
      # theme
      scale_shape_manual(values=c(16, 0, 18))+
      scale_x_continuous(breaks = seq(0.2, 0.4, by = 0.1), 
                         limits = c(0.18, 0.42)) +
      ylim(minN, maxN) +
      myThemePdf
    
    
  }
}


#legend
p <- biasSummaryMeansTmp %>% 
  ggplot(aes(x=p, y=N_estim))+
  geom_point(aes(shape = alpha)) +
  scale_shape_manual(values=c(16, 0, 17))+
  theme_bw()
legAlpha <- get_legend(p)
p2 <- biasSummaryMeansTmp %>% 
  select(S, Prior, p, alpha, N_2.5, N_97.5) %>% 
  rename("q 2.5" = N_2.5,
         "q 97.5" = N_97.5) %>% 
  pivot_longer(cols = c("q 97.5", "q 2.5"), 
               names_to = "limit", values_to = "estim") %>% 
  mutate(limit = factor(limit, 
                        levels = c("q 97.5", "q 2.5"))) %>% 
  ggplot(aes(x=p, y=estim))+
  geom_point(aes(shape = limit)) +
  scale_shape_manual(values=c(24, 25))+
  guides(shape = guide_legend(title = "Average\nestimated\nquantile")) +
  theme_bw()
leg95 <- get_legend(p2)
legend <- plot_grid(NULL, legAlpha, leg95, NULL, ncol = 1, 
                    rel_heights = c(1, 1, 1, 1))

# labels of the prior used
prior1 <- ggdraw() + draw_label(paste(priors[1], "prior"), size = labelSize)
prior2 <- ggdraw() + draw_label(paste(priors[2], "prior"), size = labelSize)
priorGrid <- plot_grid(prior1, prior2, ncol = 2)

# label of the number of occasion
s1 <- ggdraw() + draw_label("5", size = labelSize)
s2 <- ggdraw() + draw_label("7", size = labelSize)
s3 <- ggdraw() + draw_label("9", size = labelSize)
sGrid <- plot_grid(s1, s2, s3, ncol = 1)

# grids of plots
plotGrid5 <- plot_grid(plotlist = plts[["5"]], ncol = 2)
plotGrid7 <- plot_grid(plotlist = plts[["7"]], ncol = 2)
plotGrid9 <- plot_grid(plotlist = plts[["9"]], ncol = 2)
plotGrid <- plot_grid(plotGrid5, plotGrid7, plotGrid9, ncol = 1)

# grid with labels and boxplots
pGrid <- plot_grid(NULL, priorGrid, NULL,
                   sGrid, plotGrid, legend,
                   rel_widths = c(0.05, 1, 0.15), 
                   rel_heights = c(0.05, 1))
  
plot(pGrid)

rm(plts, plotGrid, plotGrid5, plotGrid7, plotGrid9, pGrid, sGrid, p, p2, legend, priorGrid)

