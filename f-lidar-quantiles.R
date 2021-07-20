# Quantile smoothing (Lidar data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(colorspace)
library(SemiPar)
library(JOPS)

# Get the data
data(lidar)
m = length(lidar$range)
sel = seq(1, m, by = 1)
x = lidar$range[sel]
y = lidar$logratio[sel]

# Compute the B-spline basis
deg = 3
xlo = min(x)
xhi = max(x)
ndx = 20
B = bbase(x, xlo, xhi, ndx, deg)
n = ncol(B)

d = 2
D = diff(diag(n), diff = d)
lambda = 10

# Fit and plot ALS curves
pp <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.98, 0.99)
pp <- c(0.02, 0.1, 0.5, 0.9, 0.98)
# pp = c(0.2, 0.8)
np <- length(pp)
A <- NULL
beta = 0.01 * (max(y) - min(y))
beta = 0.001

# Iterate to solution
P = lambda * t(D) %*% D
W = diag(length(y))
for (p in pp) {
    a = 0
    for (it in 1:50) {
        anew = solve(t(B) %*% W %*% B + P, t(B) %*% W %*% y)
        da = max(abs(a - anew))
        a = anew
        if (max(da) < 1e-04)
            break
        yhat = B %*% a
        r = y - yhat
        w = 1/sqrt(r^2 + beta^2)
        w = w * ifelse(r > 0, p, 1 - p)
        W = diag(c(w))
        # cat(p, it, da, max(abs(a)), da/ max(abs(a)), '\n')
    }
    cat(p, "\n")
    A = cbind(A, a)
}

# Compute fitted curves on grid
ng = 200
xg <- seq(from = min(x), to = max(x), length = ng)
Bg <- bbase(xg, xlo, xhi, ndx, deg)
Z = Bg %*% A

# Create data frames for ggplot
Data = data.frame(x, y)
Zf = data.frame(x = rep(xg, np), y = as.vector(Z), tau = as.factor(rep(pp,
    each = ng)))

# Build the graph
plt1 = ggplot(Data, aes(x = x, y = y)) +
       geom_point(size = 1) +
       geom_line(data = Zf, aes(x = x, y = y, group = tau, color = tau), size = 1, lty = 1) +
       ggtitle("LIDAR data") +
       xlab("Distance") +
       ylab("log(Ratio)") +
       JOPS_theme() +
       scale_color_manual(values = rainbow_hcl(np, start = 10,  end = 350))

# Make graph and pdf
print(plt1)
