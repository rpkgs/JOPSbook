# Binomial P-spline smoothing (Kyphosis data)
# A graph in the book "Practical Smoothing. The Joys of P-splines"
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)
library(rpart)

# Extract data
Kyphosis <- kyphosis$Kyphosis
Age <- kyphosis$Age
Number <- kyphosis$Number
Start <- kyphosis$Start

y <- 1 * (Kyphosis == "present") # make y 0/1

nseg <- 20
bdeg <- 3
fit1 <- psBinomial(Age, y, nseg = 20, bdeg = 3, pord = 2, lambda = 1, show = F)
fit2 <- psBinomial(Age, y, nseg = 20, bdeg = 3, pord = 2, lambda = 100, show = F)

xl <- min(Age)
xr <- max(Age)
n <- length(fit1$pcoef)
knots <- ((1:n) - 2) / nseg
knots <- xl + knots * (xr - xl)

F1 <- data.frame(Age, y)
F2 <- data.frame(
  xg1 = fit1$xgrid, xg2 = fit2$xgrid, yg1 = fit1$pgrid, yg2 = fit2$pgrid,
  eta1 = fit1$ygrid, eta2 = fit2$ygrid
)
F3 <- data.frame(knots = knots, pcoef1 = fit1$pcoef, pcoef2 = fit2$pcoef)
plt1 <- ggplot(F1) +
  geom_point(data = F1, aes(x = Age, y = y), size = 1.5, color = "darkgray") +
  geom_line(aes(x = xg1, y = yg1), data = F2, size = 1, colour = I("blue"), lty = 2) +
  geom_line(aes(x = xg2, y = yg2), data = F2, size = 1, colour = I("red"), lty = 1) +
  xlab("Age (months)") + ylab("P(presence)") +
  ggtitle("Binomial P-spline fits") + # ylim(c(0, 15)) +
  JOPS_theme()

plt2 <- ggplot(F2, aes(x = Age, y = eta1)) +
  geom_point(data = F3, aes(x = knots, y = pcoef1), size = 1.5, pch = 25, colour = I("blue")) +
  geom_point(data = F3, aes(x = knots, y = pcoef2), size = 1.5, pch = 24, colour = I("red")) +
  geom_line(aes(x = xg1, y = eta1), data = F2, size = 1, colour = I("blue"), lty = 2) +
  geom_line(aes(x = xg2, y = eta2), data = F2, size = 1, colour = I("red"), lty = 1) +
  xlab("Age (months)") + ylab("logit(p)") + ggtitle("Linear predictor") + # ylim(c(0, 15)) +
  JOPS_theme()

# Save graph
grid.arrange(plt1, plt2, ncol = 2, nrow = 1)
