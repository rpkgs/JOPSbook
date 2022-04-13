# Illustration of interpolation and extrapolation by penalty order
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)

# Simulate data
m = 50
set.seed(123)
x = runif(m)
y = sin(2.5 * x) + rnorm(m) * 0.05 + 0.2
f = which(0.2 < x & x < 0.4 | 0.6 < x & x < 0.8)
x = x[f]
y = y[f]
Data = data.frame(x, y)

# Make a matrix containing the B-spline basis
ndx = 25
deg = 3
B = bbase(x, 0, 1, nseg = ndx, bdeg = deg)
nb = ncol(B)

# Basis for fit on grid
ng = 500
xg = seq(0, 1, length = ng)
Bg = bbase(xg, 0, 1, nseg = ndx, bdeg = deg)

# Fit
D = diff(diag(nb), diff = 1)
P = t(D) %*% D
lambda = 1
a = solve(t(B) %*% B + lambda * P, t(B) %*% y)
z = Bg %*% a

knots = ((1:nb) - 2)/ndx
Fa = data.frame(x = knots, y = a)

# Create data frame for ggplot
Zf = data.frame(x = xg, y = z)

plt1 = ggplot() +
  geom_line(data = Zf, aes(x = x, y = y), size = 0.6, col = 'blue') +
  ggtitle("First differences") +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(data = Data, aes(x = x, y = y), color = "darkgrey") +
  geom_point(data = Fa, aes(x = x, y = y), color = "red", shape = 1) +
  xlab("") + ylab("") +
  JOPS_theme() +
  theme(legend.position = "none")

# Fit
D = diff(diag(nb), diff = 2)
P = t(D) %*% D
lambda = 1;
a = solve(t(B) %*% B + lambda * P , t(B) %*% y)
z = Bg %*%  a

# Create data frame for ggplot
Zf = data.frame(x = xg, y = z, id  = as.factor(1))
Fa = data.frame(x = knots, y = a, id = 1)

plt2 =   ggplot(Zf) +
  geom_line(aes(x = xg, y = y), size = 0.6, color = 'blue') +
  ggtitle("Second differences") +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(data = Data, aes(x =x, y = y), color = "darkgrey") +
  geom_point(data = Fa, aes(x = x, y = y), color = "red", shape = 1) +
  xlab("") + ylab("") +
  JOPS_theme() +
  theme(legend.position = "none")


# Show the graphs
grid.arrange(plt1, plt2, ncol = 2, nrow = 1)
