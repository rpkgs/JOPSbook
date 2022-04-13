# Impulse response of the Whittaker smoother with harmonic penalty
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(JOPS)
library(gridExtra)


# Create the impulse
n = 200
y = rep(0, n)
y[n/2] = 1

psis = c(1, 0.98, 0.95, 0.9)

# Apply the Whittaker smoother with a harmonic penalty
E = diag(n)
plts = list()
for (jp in 1:4) {
  psi = psis[jp]
dp = c(1, -2 * psi, 1)
D = matrix(0, n - 2, n)
for (j in 1:(n - 2)) D[j, j:(j + 2)] = dp
Ph = t(D) %*% D
lambda1 = 1000
lambda2 = 5000
mu1 = solve(E + lambda1 * Ph, y)
mu2 = solve(E + lambda2 * Ph, y)

# Data frames for plotting
x = seq(-1, 1, length = n)
F1 = data.frame(x = x, mu1 = mu1, mu2 = mu2)

ttl = bquote('Harmonic penalty,' ~ psi == .(psi))
if (psi ==1) ttl = bquote('Standard penalty,' ~ psi == .(psi))

# Make the plot
plt = ggplot(F1)  +
  geom_line(aes(x = x, y = mu1), color = "red", size = 0.7) +
  geom_line(aes(x = x, y = mu2), color = "blue", size = 0.7) +
  ggtitle(ttl) +
  xlab('') + ylab('') +
  JOPS_theme() +
  theme(plot.title = element_text(size = 11))

  plts[[jp]] = plt
}

grid.arrange(grobs = plts, nrow = 2, ncol = 2)
