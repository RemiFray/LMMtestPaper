
a <- c(74, 85, 91)
b <- c(26, 15,9)

plts <- list()
i <- 1
x <- seq(0.5,1, 0.001)
dta <- data.frame(x = x, y = dbeta(x, a[i], b[i]))
plts[[1]] <- ggplot(dta, aes(x, y)) + 
  geom_line() +
  geom_vline(aes(xintercept = 0.8), lty = "dotted") +
  geom_vline(aes(xintercept = qbeta(0.5, a[1], b[1]))) +
  geom_vline(aes(xintercept = qbeta(0.95, a[1], b[1])), lty= "longdash") +
  xlab(paste0("Beta(", a[i], ", ", b[i], ") density")) +
  theme_bw() +
  theme(axis.title.y = element_blank())
i <- 2
x <- seq(0.65,1, 0.01)
dta <- data.frame(x = x, y = dbeta(x, a[i], b[i]))
plts[[2]] <- ggplot(dta, aes(x, y)) + 
  geom_line() +
  geom_vline(aes(xintercept = 0.9), lty = "dotted") +
  geom_vline(aes(xintercept = qbeta(0.5, a[2], b[2]))) +
  geom_vline(aes(xintercept = qbeta(0.95, a[2], b[2])), lty= "longdash") +
  xlab(paste0("Beta(", a[i], ", ", b[i], ") density")) +
  theme_bw() +
  theme(axis.title.y = element_blank())
i <- 3
x <- seq(0.75,1, 0.01)
dta <- data.frame(x = x, y = dbeta(x, a[i], b[i]))
plts[[3]] <- ggplot(dta, aes(x, y)) + 
  geom_line() +
  geom_vline(aes(xintercept = 0.95), lty = "dotted") +
  geom_vline(aes(xintercept = qbeta(0.5, a[3], b[3]))) +
  geom_vline(aes(xintercept = qbeta(0.95, a[3], b[3])), lty= "longdash") +
  xlab(paste0("Beta(", a[i], ", ", b[i], ") density")) +
  theme_bw() +
  theme(axis.title.y = element_blank())

plot(plot_grid(plotlist = plts, ncol = 3))
