# Two-dimensional isotropic density estimation (Old Faithful data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(reshape2)
library(ucminf)
library(JOPS)
library(ggplot2)

# Get the data
data(faithful)
u = faithful[, 1]
v = faithful[, 2]

# Create 2-d histogram
h = hist2d(u, v, xlim = c(1, 6), ylim = c(40, 100))
Y = h$H

# Smooth histogram and return AIC
aicfun = function(lla) {
  lambda = 10^lla
  fit = hist2dsm(Y, lambdax = lambda, lambday = lambda, dx = 2)
  Mu = fit$Mu
  ok = Y > 0
  dev = 2 * sum(Y[ok] * log(Y[ok]/Mu[ok]))
  ed = fit$ed
  aic = dev + 2 * ed
}

# Search for best (log) lambdas
op = ucminf(0, aicfun)
lambda = 10^op$par
cat("log10(lambdas), AIC:", op$par, op$value, "\n")
fit = hist2dsm(Y, lambdax = lambda, lambday = lambda, dx = 2)
Mu = fit$Mu

# Turn matrix into a 'long' data frame
row.names(Mu) = h$xgrid
names(Mu) = h$ygrid
dens <- melt(sqrt(Mu))
names(dens) = c("x", "y", "z")

# Plot density with contours
plt = ggplot(dens,  aes(x, y, fill = z)) +
  geom_raster(show.legend = F) +
  scale_fill_gradient(high = 'darkgreen', low = 'white') +
  xlab('Eruption length (min)') + ylab('Waiting time (min)') +
  ggtitle('Old Faithful, isotropic smoothing (square root of density)') +
  geom_contour(aes(z = z), color = "steelblue", show.legend = T) +
  JOPS_theme()

# Make and save plot
plot(plt)
save_PDF("f-geyser-iso")


