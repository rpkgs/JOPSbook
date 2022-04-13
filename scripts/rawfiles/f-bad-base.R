# Compare proper and improper B-spline basis
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(splines)
library(JOPS)

# Set domain
xmax = 5
xmin = 0
m = 500
x = seq(xmin, xmax, length = m)

# Compute a proper basis
knots = (xmin - 3):(xmax + 3)
B <- splineDesign(knots, x = x, outer.ok = T)
B[B < 1e-05] = NaN
n = ncol(B)

# Compute an improper basis
knots1 = 0:(xmax - 1)
B1 <- bs(x = x, knots = knots1)
B1[B1 < 1e-05] = NaN

# Make data frames for ggplot
Bf = data.frame(x = rep(x, n), y = as.vector(B),
                id = as.factor(rep(1:n, each = m)))
Bf1 = data.frame(x = rep(x, n), y = as.vector(B1),
                 id = as.factor(rep(1:n, each = m)))

# Produce ggplot objects
lim = c(0,1)
plt1 = ggplot(Bf, aes(x = x, y = y, color = id)) +
  ylim(lim) +
  geom_line(size = 1)  +
  geom_hline(yintercept = 0, color = 'grey') +
  JOPS_theme() +
  xlab("") + ylab("") +
  ggtitle('Proper B-spline basis') +
  theme(legend.position = "none")

plt2 = ggplot(Bf1, aes(x = x, y = y, color = id), ylim = lim) +
  ylim(lim) +
  geom_line(size = 1)  +
  geom_hline(yintercept = 0, color = 'grey') +
  xlab("") + ylab("") +
  ggtitle('Improper B-spline basis') +
  JOPS_theme()  +
  theme(legend.position = "none")

# Make and save plot
grid.arrange(plt1, plt2, nrow = 1, ncol = 2)
