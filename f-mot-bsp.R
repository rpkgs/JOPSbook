# B-spline fits having differing support (Motorcycle data)
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

# Boundary for the subdomain
thr = 5
sel = x > thr
xsel = x[sel]
ysel = y[sel]

# Compute the B-spline basis
deg = 3
xlo = min(x)
xhi = max(x)
ndx = 5
B = bbase(x, xlo, xhi, nseg = ndx, bdeg = deg)

# Basis for fit on grid
ng = 1000
xg = seq(min(x), max(x), length = ng)
Bg = bbase(xg, xlo, xhi, nseg = ndx, bdeg = deg)

# Use 0/1 weight to select the subdomain
W = 1 * diag(x > thr)

# Estimate the coefficients and compute fit on the grid
a = solve(t(B) %*% B, t(B) %*% y)
z1 = Bg %*% a
asel = solve(t(B) %*% W %*% B, t(B) %*% W %*% y)
zsel = Bg %*% asel

# Create data frames for ggplot
Zf1 = data.frame(x = xg, y = z1, id = as.factor(1))
Zf2 = data.frame(x = xg[xg > thr], y = zsel[xg > thr], id = as.factor(1))

# Build the graph
plt1 = ggplot(Data, aes(x = x, y = y)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_line(data = Zf1, size = 1, colour = I("blue")) +
  geom_line(data = Zf2, size = 1, colour = I("red"), lty = 6) +
  xlab("Time (ms)") + ylab("Acceleration (g)") +
  ggtitle("Motorcycle helmet impact data") +
  ylim(c(-150, 100)) +
  JOPS_theme()

print(plt1)

