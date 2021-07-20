# B-splines in perspective
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(colorspace)
library(JOPS)

# Basis on grid
ndx = 4
deg1 = 3
ng = 500
xmin = 0
xmax = 4
xg = seq(xmin, xmax, length = ng)
Bg = bbase(xg, xmin, xmax, nseg = ndx, bdeg = deg1)
nb1 = ncol(Bg)

# Make a matrix with B-splines scaled by coefficients
Bsc1 = Bg + outer(rep(1, ng), 1:ncol(Bg))

# Create data frames for ggplot
Bf1 = data.frame(x = rep(xg, nb1), y = as.vector(Bsc1), id = as.factor(rep(1:nb1,
    each = ng)))
Bf1$y[abs(Bf1$y) < 1e-04] = NA
Bf1 = na.omit(Bf1)

# Select one row, for visualization
k = 160
xk = xg[k]
bk2 = Bg[k, ]
bk1 = bk2 + 1:nb1
bk2[bk2 < 0.001] = NaN
Fk1 = data.frame(x = rep(xk, nb1), y = bk1, id = as.factor(1:nb1))
Fk2 = data.frame(x = rep(xk, nb1), y = bk2, id = as.factor(1:nb1))

# Build the graphs
plt1 = ggplot(Bf1,aes(x = x, y = y, group = id, colour = id)) +
  geom_line(size = 1.2) +
  ggtitle("Perspective view") +
  geom_hline(yintercept = 0, size = 0.3) +
  xlab("") + ylab("") +
  JOPS_theme() +
  geom_vline(xintercept = xk, lty = 2) +
  geom_point(aes(colour = id), size = 3, data = Fk1) +
  theme(legend.position = "none") +
  scale_color_manual(values = rainbow_hcl(nb1 + 1, start = 10, end = 350))

# Create data frames for ggplot
Bf2 = data.frame(x = rep(xg, nb1), y = as.vector(Bg),
                 id = as.factor(rep(1:nb1, each = ng)))
Bf2$y[abs(Bf2$y) < 0.0001] = NA
Bf2 = na.omit(Bf2)

# Build the graphs
plt2 = ggplot(Bf2,aes(x = x, y = y, group = id, colour = id)) +
  geom_line(size = 1.2) +
  ggtitle("Columns of a B-spline basis") +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_vline(xintercept = xk, lty = 2) +
  geom_point(aes(colour = id), size = 3,  data = Fk2) +
  xlab("") + ylab("") +
  JOPS_theme() +
  theme(legend.position = "none") +
  scale_color_manual(values = rainbow_hcl(nb1 + 1, start = 10, end = 350))

grid.arrange(plt1, plt2, ncol = 2, nrow = 1)
