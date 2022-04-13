# Histogram smoothing with INLA (Old Faithful geyser data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(INLA)
library(JOPS)
library(ggplot2)
library(gridExtra)

# First plot ======

# Compute the histogram
data(faithful)
h = hist(faithful$eruptions, breaks = seq(1, 6, by = 0.05), plot = F)
x = h$mids
y = h$counts
m = length(x)

# Prepare the B-spline basis
xlo = min(x)
xhi = max(x)
B = bbase(x, xlo, xhi, nseg = 20)
n = ncol(B)

# Compute fit with INLA
data.inla = list(y = y, x = 1:n)
formula = y ~ -1 + f(x, model = "rw2", constr = F)
result = inla(formula, family = "poisson", data = data.inla, control.predictor = list(A = B,
    compute = T))

# Extract coefficient statistics
S = result$summary.linear.predictor[m + (1:n), ]
ahat = S$mean

# Compute fit on fine grid
xg = seq(xlo, xhi, length = 200)
Bg = bbase(xg, xlo, xhi, nseg = 20)
mu = exp(Bg %*% ahat)

Data = data.frame(x = x, y = y)
Fit = data.frame(x = xg, y = mu)

plt1 = ggplot(aes(x = x, y = y,  fill = I("wheat3")), data = Data) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  #ggtitle('Bayesian fit to Old Faithful histogram, using INLA') +
  xlab("Eruption length (min)") + ylab("Frequency") +
  geom_line(data = Fit, col = I("steelblue"), size = 1) +
  JOPS_theme()

# Second plot=====

# Compute the histogram
data(faithful)
h = hist(faithful$eruptions, breaks = seq(1,6, by =0.1), plot = F)
x = h$mids
y = h$counts
m = length(x)

# Prepare the B-spline basis
xlo = min(x)
xhi = max(x)
B = bbase(x, xlo, xhi, nseg = 20)
n = ncol(B)

# Compute fit with INLA
data.inla = list(y = y, x = 1:n)
formula = y ~ -1 + f(x, model = 'rw2', constr = F)
result = inla(formula, family = 'poisson', data = data.inla, control.predictor = list(A = B, compute = T))

# Extract coefficient statistics
S = result$summary.linear.predictor[m + (1:n), ]
ahat = S$mean

# Compute fit on fine grid
xg = seq(xlo, xhi, length = 200)
Bg = bbase(xg, xlo, xhi, nseg = 20)
mu = exp(Bg %*% ahat)

Data = data.frame(x = x, y = y)
Fit = data.frame(x = xg, y = mu)

plt2 = ggplot(aes(x = x, y = y,  fill = I("wheat3")), data = Data) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  #ggtitle('Bayesian fit to Old Faithful histogram, using INLA') +
  xlab("Eruption length (min)") + ylab("Frequency") +
  geom_line(data = Fit, col = I("steelblue"), size = 1) +
  JOPS_theme()

# Save graphs
grid.arrange(plt1, plt2, ncol = 2, nrow = 1)
