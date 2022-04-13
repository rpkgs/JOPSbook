# Smoothing the sum of absolute differences in penalty (Chromosome data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(Vega)
library(spam)
library(JOPS)
library(ggplot2)
library(gridExtra)

# Get the data
data(G519)
chr = G519[, 1]
sel = chr == '18'
x = as.numeric(G519[sel, 2])
y = as.numeric(G519[sel, 4])

# Prepare for smoothing
m = length(y)
E = diag.spam(m)
D = diff(E)
v = rep(1, m - 1)
beta = 0.001
nit = 200
Z = matrix(0, m, nit)
z = 0

lambdas = c(0.1, 0.2, 1)
Z = NULL

# Do the smoothing
for (lambda in lambdas) {
  v = rep(1, m - 1)
  for (it in 1:nit) {
    V = diag.spam(as.vector(v))
    P = lambda * t(D) %*% V %*% D
    znew = solve(E + P, y)
    dz = max(abs(znew - z))
    z = znew
    if (dz < 1e-4) break
    g = diff(z)
    v = 1 / (g ^ 2 + beta ^ 2)
  }
  Z = cbind(Z, z)
  cat(lambda, '\n')
}

# Plot results
pp = list()
Dxy = data.frame(x = x, y = y)

# Make and collect the graphs
for (k in 1:3) {
  Dz = data.frame(x = x, z = Z[, k])
  plt = ggplot(Dxy, aes (x = x, y = y)) +
    geom_point(col = 'darkgrey', size = 0.6) +
    geom_line(col = 'darkgrey') +
    geom_line(data = Dz, aes(x = x, y = z), color = "blue") +
    xlab('Position (Mb)') + ylab('LRR') +
    ggtitle(bquote("G519, chromosome 18," ~ lambda == .(lambdas[k]))) +
    JOPS_theme()
    pp[[k]] = plt
}

# Save pdf
grid.arrange(grobs = pp, nrow = 3, ncol = 1)

