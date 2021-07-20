# B-spline fits with same basis having varying roughness
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)
library(colorspace)

# Make basis
nu = 200
u = seq(0, 1, length = nu)
nseg = 10
B = bbase(u, nseg = nseg)
n = ncol(B)
set.seed(123)

# Make coefficients
a1 = runif(n)
a2 = 0.8 * sin(2 * (1:n)/n) + runif(n) * 0.2
a3 = (1:n)/n
a4 = rep(1, n)
A = cbind(a1, a2, a3, a4)
Z = B %*% A

# Generate the plots
plts = list()
for (j in 1:4) {

  # Compute roughness
  aj = c(A[, j])
  da = diff(aj)
  r = sqrt(sum(da ^ 2)/(n - 1))
  sr = sprintf("r = %3.2f", r)

  # Scaled basis
  Bsc = B %*% diag(aj)

  # Create data frames for ggplot
  id = as.factor(rep(1:n, each = nu))
  Zf = data.frame(x = u, y = Z[, j], id  = as.factor(1))
  Bf = data.frame(x = rep(u, n), y = as.vector(Bsc), id = id)

  # Remove zero entries
  Bf$y[Bf$y < 0.0001] = NA
  Bf = na.omit(Bf)

  knots = ((1:n) - 2) / nseg
  Fa = data.frame(x = knots, y = aj, id = 1)

  # Build the graphs
  plt1 = ggplot(Bf, aes(x = x, y = y, group = id, colour = id)) +
    geom_line( size = 0.7) +
    ggtitle(sr) +
    geom_hline(yintercept = 0, size = 0.3) +
    geom_line(data = Zf, size = 1, col = "blue") +
    geom_point(data = Fa, color = "red", size = 2, shape = 1) +
    xlab("") + ylab("") +
    JOPS_theme() +
    theme(legend.position = "none") +
    scale_color_manual(values=rainbow_hcl(n, start = 10, end = 350))

  # Add to list of plots
  plts[[j]] = plt1
}

# Plot and save
grid.arrange(grobs = plts, nrow = 2, ncol = 2)
