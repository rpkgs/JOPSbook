# Illustration of linear and quadratic B-spline bases
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(colorspace)
library(JOPS)

# Basis on grid
ndx = 4
deg1 = 1
deg2 = 2
ng = 500
xmin = 0
xmax = 4
xg = seq(xmin, xmax, length = ng)
Bg1 = bbase(xg, xmin, xmax, nseg = ndx, bdeg = deg1)
Bg2 = bbase(xg, xmin, xmax, nseg = ndx, bdeg = deg2)
nb1 = ncol(Bg1)
nb2 = ncol(Bg2)

# Create data frames for ggplot
Bf1 = data.frame(x = rep(xg, nb1), y = as.vector(Bg1),
                 id = as.factor(rep(1:nb1, each = ng)))
Bf1$y[abs(Bf1$y) < 0.0001] = NA
Bf1 = na.omit(Bf1)
Bf2 = data.frame(x = rep(xg, nb2), y = as.vector(Bg2),
                 id = as.factor(rep(1:nb2, each = ng)))
Bf2$y[abs(Bf2$y) < 0.0001] = NA
Bf2 = na.omit(Bf2)
Fk = data.frame(x = seq(0, 1, length = 5), y = rep(0, 5), id = 1)

# Build the graphs
plt1 = ggplot(Bf1,aes(x = x, y = y, group = id, colour = id)) +
  geom_line(size = 1.2) +
  ggtitle("Linear B-splines") +
  geom_hline(yintercept = 0, size = 0.3) +
  xlab("") + ylab("") +
  JOPS_theme() +
  theme(legend.position = "none") +
  scale_color_manual(values = rainbow_hcl(nb1 + 1, start = 10, end = 350))

# Build the graphs
plt2 = ggplot(Bf2,aes(x = x, y = y, group = id, colour = id)) +
  geom_line(size = 1.2) +
  ggtitle("Quadratic B-splines") +
  geom_hline(yintercept = 0, size = 0.3) +
  xlab("") + ylab("") +
  JOPS_theme() +
  theme(legend.position = "none") +
  scale_color_manual(values = rainbow_hcl(nb1 + 1, start = 10, end = 350))

# Show the graphs
grid.arrange(plt1, plt2, ncol = 2, nrow = 1)
