# Double penalization for impulse response without negative side lobes
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(JOPS)

# Simulate data
m0 = 100
x = ((1:m0) - 0.5)/m0
y = 0 * x
y[51] = 1
m = length(x)

# Prepare basis and penalty matrices
nseg = 100
B = bbase(x, xl = 0, xr = 1, nseg = nseg)
n = ncol(B)
D1 = diff(diag(n))
D2 = diff(D1)

# Do the smoothing
lambda = 1000
gams = c(0, 0.2, 0.5, 1, 2)
nc = length(gams)
Z = NULL
for (gam in gams) {
  P = lambda * t(D2) %*% D2 + gam * sqrt(lambda) * t(D1) %*% D1
  a = solve(t(B) %*% B + P, t(B) %*% y)
  z = B %*% a
  Z = cbind(Z, z)
}

# Make data frames for ggplot
Fz = data.frame(x = rep(x, nc), y = as.vector(Z), gamma = as.factor(rep(gams,
    each = m)))

# Make the plot
plt = ggplot(data = Fz) +
      geom_line(aes(x = x, y = y, color = gamma), size = 1) +
      ggtitle("Impulse response with double penalty") +
      scale_color_discrete(name = bquote(~ gamma)) +
      JOPS_theme()

# Save graph
print(plt)

