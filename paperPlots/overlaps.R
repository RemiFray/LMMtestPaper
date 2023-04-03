allOvlaps <- read.csv("./data/LMMPt_allOverlaps.csv")


plt <- allOvlaps %>% 
  mutate(ap = paste0("a", a, " p", p),
         Prior = fct_recode(Prior, "Uninformative" = "noPrior",
                            "Unbiased" = "priorAlpha100i",
                            "Biased" = "priorAlpha100iBiased")) %>%
  ggplot() +
  geom_boxplot(aes(x = ap, y=overlap), outlier.alpha = 0) +
  geom_hline(aes(yintercept = 0.35)) +
  facet_wrap(vars(Prior), labeller = label_wrap_gen(multi_line=FALSE))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        axis.title = element_blank())

plot(plt)
