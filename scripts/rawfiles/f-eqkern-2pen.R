# Show equivalent kernels for Whittaker smoother with double penalty
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(JOPS)
library(gridExtra)

# Create the impulse
n <- 201
y <- rep(0, n)
y[n / 2] <- 1
lambda= 100
gammas = c( 0, 1, 2, 4)

# Apply the Whittaker smoother with a harmonic penalty
E <- diag(n)
plts <- list()
cols = rainbow(5)
for (jp in 1:4) {
  gamma = gammas[jp]
  D <- diff(E, diff = 2)
  D1 = diff(E)
  P <- lambda * t(D) %*% D + gamma * sqrt(lambda) * t(D1) %*% D1
  H = solve(E + P)

  # Data frames for plotting
  h = as.vector(H[, seq(1, n, length = 5)])
  x <- 1:n
  F1 <- data.frame(x = rep(x, 5), y = h, id = as.factor(rep(1:5, each = n)))

  ttl = bquote(gamma==.(gamma))

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
