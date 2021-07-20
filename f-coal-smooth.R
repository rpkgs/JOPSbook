# Poisson smoothing (Coal mining data)
# A graph in the book "Practical Smoothing. The Joys of P-splines"
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)
library(boot)

# Extract the data
Count <- hist(coal$date, breaks = c(1851:1963), plot = F)$counts
Year <- c(1851:1962)
xl <- min(Year)
xr <- max(Year)

# Poisson smoothing
nseg <- 20
bdeg <- 3
fit1 <- psPoisson(Year, Count, xl, xr, nseg, bdeg, pord = 2, lambda = 1, show = F)
fit2 <- psPoisson(Year, Count, xl, xr, nseg, bdeg, pord = 2, lambda = 100, show = F)

n <- length(fit1$pcoef)
knots <- ((1:n) - 2) / nseg
knots <- xl + knots * (xr - xl)

# Plotting on fine grid
F1 = data.frame(Year, Count)
F2 = data.frame(xg1 = fit1$xgrid,xg2 = fit2$xgrid,yg1 = fit1$mugrid, yg2 = fit2$mugrid,
                eta1 = fit1$ygrid, eta2 = fit2$ygrid)
F3=data.frame(knots = knots, pcoef1 = fit1$pcoef, pcoef2 = fit2$pcoef)
plt1 = ggplot(F1, aes(x = Year, y = Count)) +
  geom_point(data=F1,size = 1.5,,color = "darkgray") +
  geom_line(aes(x = xg1, y = yg1), data = F2, size = 1, colour = I("blue"), lty = 2) +
  geom_line(aes(x = xg2, y = yg2), data = F2, size = 1, colour = I("red"), lty = 1) +
  xlab(" ") + ylab("Accident count") +
  ggtitle("Poisson P-spline fits") +
  JOPS_theme()

plt2 = ggplot(F1, aes(x = Year, y = Count)) +
  geom_point(data = F3,aes(x = knots, y = pcoef1), size = 1.5, pch = 25, colour = I("blue")) +
  geom_point(data = F3,aes(x = knots, y = pcoef2), size = 1.5, pch = 24, colour = I("red")) +
  geom_line(aes(x = xg1, y = eta1), data = F2, size = 1, colour = I("blue"), lty = 2) +
  geom_line(aes(x = xg2, y = eta2), data = F2, size = 1, colour = I("red"), lty = 1) +
  xlab(" ") + ylab(" ") +
  ggtitle("Linear predictor") +
  JOPS_theme()

grid.arrange(plt1, plt2, ncol = 2, nrow = 1)
