# Show equivalent kernels for the harmonic Whittaker smoother
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(JOPS)
library(gridExtra)

# Apply the Whittaker smoother with a harmonic penalty
n <- 201
E <- diag(n)
plts <- list()
lambda = 300
psis = c(0.95, 0.9, 0.8, 0.7)

for (jp in 1:4) {
  psi = psis[jp]
  dp = c(1, -2 * psi, 1)
  D = matrix(0, n - 2, n)
  for (j in 1:(n - 2)) D[j, j:(j + 2)] = dp
  P = lambda * t(D) %*% D
  H = solve(E + P)

  # Prepare for ggplot
  h = as.vector(H[, seq(1, n, length = 3)])
  x <- 1:n
  F1 <- data.frame(x = rep(x, 3), y = h, id = as.factor(rep(1:3, each = n)))

  ttl = bquote(psi == .(psi) ~  lambda == .(lambda))
  cols = rainbow(5)

  # Make the plot
  plt <- ggplot(F1, aes(x = x, y = y, color = id)) +
    geom_hline(yintercept = 0, size = 0.3) +
    geom_line() +
    ggtitle(ttl) +
    xlab("") + ylab("") +
    JOPS_theme() +
    theme(plot.title = element_text(size = 11)) +
    theme(legend.position = "none")
  plts[[jp]] <- plt
}

grid.arrange(grobs = plts, nrow = 2, ncol = 2)
