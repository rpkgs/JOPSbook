# Circular B-spline basis
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(colorspace)
library(JOPS)

# Basis on grid
ndx = 6
deg = 3
ng = 500
xg = seq(0, 1, length = ng)
Bg = cbase(xg, 0, 1, nseg = ndx, bdeg = deg)
nb1 = ncol(Bg)
Bshift = Bg + outer(rep(1, ng), (1:nb1))

# Create data frame for ggplot
Bf1 = data.frame(x = rep(xg, nb1), y = as.vector(Bshift),
                 id = as.factor(rep(1:nb1, each = ng)))
Bf1$y[abs(Bf1$y) < 0.0001] = NA
Bf1 = na.omit(Bf1)

# Create data frame for ggplot
Bf2 = data.frame(x = rep(xg, nb1), y = as.vector(Bg),
                 id = as.factor(rep(1:nb1, each = ng)))
Bf2$y[abs(Bf2$y) < 1e-8] = NA
Bf2 = na.omit(Bf2)

# Build the graphs
plt1 = ggplot(Bf1,aes(x = x, y = y, colour = id)) +
  geom_line(size = 1.2) +
  ggtitle("Perspective view") +
  geom_hline(yintercept = 0, size = 0.3) +
  xlab("") + ylab("") +
  JOPS_theme() +
  theme(legend.position = "none") +
  scale_color_manual(values = rainbow_hcl(nb1 + 1, start = 10, end = 350))

plt2 = ggplot(Bf2,aes(x = x, y = y, colour = id)) +
  geom_line(size = 1.2) +
  ggtitle("Standard view") +
  geom_hline(yintercept = 0, size = 0.3) +
  xlab("") + ylab("") +
  JOPS_theme() +
  theme(legend.position = "none") +
  scale_color_manual(values = rainbow_hcl(nb1 + 1, start = 10, end = 350))

# Show and save graph
grid.arrange(plt1, plt2, ncol = 2, nrow = 1)
