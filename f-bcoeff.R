# View of B-spline coefficients
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(colorspace)
library(JOPS)

# Simulate data
n = 50
set.seed(2016)
x = runif(n)
y0 = 0.3 + sin(1.2 * x + 0.3)
noise = rnorm(n)

# Make a matrix containing the small B-spline basis
ndx = 10
deg = 3
B = bbase(x, 0, 1, nseg = ndx, bdeg = deg)
nb = ncol(B)
knots = ((1:nb) - 2)/ndx
D = diff(diag(nb), diff = 2)
lambda1 = 0.1

# Basis for fit on grid
ng = 500
xg = seq(0, 1, length = ng)
Bg = bbase(xg, 0, 1, nseg = ndx, bdeg = deg)

# Estimate the coefficients and compute fit on the grid
y = y0 + noise * 0.1
a1 = solve(t(B) %*% B + lambda1 * t(D) %*% D, t(B) %*% y)
z1 = Bg %*% a1
A1 = data.frame(x = knots, y = a1, id = as.factor(1))

# Make a matrix with B-splines scaled by coefficients
Bsc1 = Bg %*% diag(c(a1))

# Create data frames for ggplot
Zf1 = data.frame(x = xg, y = z1, id = as.factor(1))
Bf1 = data.frame(x = rep(xg, nb), y = as.vector(Bsc1), id = as.factor(rep(1:nb,
    each = ng)))
Bf1$y[abs(Bf1$y) < 1e-04] = NA
Bf1 = na.omit(Bf1)

# Build the graphs
plt1 = ggplot(Bf1, aes(x = x, y = y, group = id, colour = id)) +
    geom_line(size = 0.7) +
    ggtitle("Wiggly curve") +
    geom_hline(yintercept = 0, size = 0.3) +
    geom_line(data = Zf1, size = 1, colour = "blue") +
    geom_point(data = A1,  color = "red", size = 1.5) +
    xlab("") + ylab("") +
    JOPS_theme() +
    theme(legend.position = "none") +
    scale_color_manual(values = rainbow_hcl(nb + 1, start = 10, end = 350))

lambda2 = 3
a2 = solve(t(B) %*% B + lambda2 * t(D) %*% D, t(B) %*% y)
z2 = Bg %*% a2
Bsc2 = Bg %*% diag(c(a2))

# Create data frames for ggplot
Zf2 = data.frame(x = xg, y = z2, id  = as.factor(1))
Bf2 = data.frame(x = rep(xg, nb), y = as.vector(Bsc2),
                 id = as.factor(rep(1:nb, each = ng)))
Bf2$y[abs(Bf2$y) < 0.0001] = NA
Bf2 = na.omit(Bf2)
A2 = data.frame(x = knots, y = a2, id = as.factor(1))

# Build the graphs
plt2 = ggplot(Bf2, aes(x = x, y = y, group = id, colour = id)) +
  geom_line( size = 0.7) +
  ggtitle("Smooth curve") +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_line(data = Zf2, size = 1, color = 'blue') +
  geom_point(data = A2, color = "red", size = 1.5) +
  xlab("") + ylab("") +
  JOPS_theme() +
  theme(legend.position = "none") +
  scale_color_manual(values = rainbow_hcl(nb, start = 10, end = 350))

# Show the graphs
grid.arrange(plt1, plt2, ncol = 2, nrow = 1)
