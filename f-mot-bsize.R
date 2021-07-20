# Compare a small and large number of B-spline segments (Motorcyle data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(colorspace)
library(MASS)
library(JOPS)

# Get the data
data(mcycle)
x = mcycle$times
y = mcycle$accel
Data = data.frame(x, y)

# Make a matrix containing the small B-spline basis
deg = 3
xlo = min(x)
xhi = max(x)
xlo = 0
xhi = 60

# Basis for fit on grid
ng = 1000
xg = seq(min(x), max(x), length = ng)

# Estimate the coefficients and compute fit on the grid Small basis
ndx = 10
B = bbase(x, xlo, xhi, nseg = ndx, bdeg = deg)
Bg = bbase(xg, xlo, xhi, nseg = ndx, bdeg = deg)
a = solve(t(B) %*% B, t(B) %*% y)
z1 = Bg %*% a

# Large basis
ndx = 20
B = bbase(x, xlo, xhi, nseg = ndx, bdeg = deg)
Bg = bbase(xg, xlo, xhi, nseg = ndx, bdeg = deg)
a = solve(t(B) %*% B, t(B) %*% y)
z2 = Bg %*% a

# Create data frames for ggplot Can be combined into one data.frame
# (tbd)
Z1 = data.frame(x = xg, y = z1, id = as.factor(1))
Z2 = data.frame(x = xg, y = z2, id = as.factor(1))

# Build the graph
plt1 = ggplot(Data, aes(x = x, y = y)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_line(aes(colour = I("red")), data = Z1, size = 1.2, linetype = 2 ) +
  geom_line(aes(colour = I("red")), data = Z1, size = 0.6, linetype = 1 ) +
  geom_line(aes(colour = I("blue")), data = Z2, size = 1., linetype = 1) +
  xlab("Time (ms)") + ylab("Acceleration (g)") +
  ggtitle("Motorcycle helmet impact data") +
  ylim(c(-200, 100)) +
  JOPS_theme()

# Make and save graph
print(plt1)
