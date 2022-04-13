# Tensor product surfaces, with various row and column tuning
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(JOPS)

# Add noise
nx = ny = 30
x = 1:nx
y = 1:ny
Z = matrix(runif(nx * ny), nx, ny)

sm2d = function(Z, W, x, y, lamdas) {
    # Prepare bases
    Bx = bbase(x, nseg = 10)
    By = bbase(y, nseg = 10)
    nbx = ncol(Bx)
    nby = ncol(By)
    W = 0 * Z + 1

    # Prpare the penalty matrices
    Dx = diff(diag(nbx), diff = 2)
    Dy = diff(diag(nby), diff = 2)
    lambdax = lambday = 1
    Px = lamdas[1] * t(Dx) %*% Dx
    Py = lamdas[2] * t(Dy) %*% Dy
    P = kronecker(Py, diag(nbx)) + kronecker(diag(nby), Px)

    # Do the smoothing, using the array algorithm
    W = 0 * Z + 1
    Tx = rowtens(Bx)
    Ty = rowtens(By)
    Q = t(Tx) %*% W %*% Ty
    dim(Q) = c(nbx, nbx, nby, nby)
    Q = aperm(Q, c(1, 3, 2, 4))
    dim(Q) = c(nbx * nby, nbx * nby)
    r = t(Bx) %*% (Z * W) %*% By
    dim(r) = c(nbx * nby, 1)
    A = solve(Q + P, r)
    dim(A) = c(nbx, nby)
    Zhat = Bx %*% A %*% t(By)
}

# Setting various tuning to show flexiblity
par(mfrow = c(2, 2), mar = c(0, 0, 0, 0))
Zhat1 = sm2d(Z, W, x, y, c(0.001, 0.001))
Zhat2 = sm2d(Z, W, x, y, c(1000, 1e-04))
Zhat3 = sm2d(Z, W, x, y, c(1e-04, 1000))
Zhat4 = sm2d(Z, W, x, y, c(1000, 1000))

thecol = 'blue'
thelwd = 0.7
xlim = ylim = c(0.15, 0.85) * nx
zlim = c(0, 1) * max(Z)
persp(x, y, Zhat1, theta = 30, phi = 50, r = 10, d = 2, axes = F, box = F,
    zlim = zlim, xlim = xlim, ylim = ylim, border = thecol, lwd = thelwd)
persp(x, y, Zhat2, theta = 30, phi = 50, r = 10, d = 2, axes = F, box = F,
    zlim = zlim, xlim = xlim, ylim = ylim, border = thecol, lwd = thelwd)
persp(x, y, Zhat3, theta = 30, phi = 50, r = 10, d = 2, axes = F, box = F,
    zlim = zlim, xlim = xlim, ylim = ylim, border = thecol, lwd = thelwd)
persp(x, y, Zhat4, theta = 30, phi = 50, r = 10, d = 2, axes = F, box = F,
    zlim = zlim, xlim = xlim, ylim = ylim, border = thecol, lwd = thelwd)





