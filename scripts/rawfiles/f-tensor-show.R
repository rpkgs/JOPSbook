# Impression of one B-spline tensor product and a tensor basis
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(JOPS)

# Create one-dimensional basis
nseg = 6
x = seq(0, 1, length = 40)
B = bbase(x, nseg = nseg)

# Combine individual tensor products
S = outer(B[, 5], B[, 5])

# Plot one tensor product
par(mfrow = c(1, 2), mar = c(0, 0, 1, 0))
zmax = 1.5 * max(S)
thecol = 'blue'
thelwd = 0.7
persp(x, x, S, theta = 40, phi = 30, d = 0.9, zlim = c(0, zmax), axes = F,
    box = F, border = thecol, lwd = thelwd, r = 10)
title("One tensor product")

# Create one-dimensional basis
nseg = 8
x = seq(0, 1, length = 60)
B = bbase(x, nseg = nseg)

# Combine individual tensor products
S = 0 * outer(x, x)
step = 2
for (j in seq(4, 9, by = step)) {
  for (k in seq(4, 9, by = step)) {
    P = outer(B[, j], B[, k])
    S = pmax(S, P)
  }
}

# Make the plot
persp(x, x, S, theta = 40, phi = 40, d = 0.9, zlim = c(0, zmax), axes = F,
    box = F, border = thecol, lwd = thelwd, r = 10)
title("A tensor product basis")
