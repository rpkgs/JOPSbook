# Convergence of the fast Harville algorithm (Motorcycle data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)
library(MASS)

# Get the data
data(mcycle)
x = mcycle$times
y = w = mcycle$accel

# Design parameters
pord = 2
bdeg = 3
nseg = 20
nit = 8

# Construct the basis matrix
B = bbase(x, nseg = 20, bdeg = 3)
n = ncol(B)
m = length(y)

# Construct the penalty matrix
D = diff(diag(n), diff = 2)
P = t(D) %*% D

# Precompute the inner products
BtB = t(B) %*% B
Bty = t(B) %*% y
lambda = 1
sigs = sig2 = 1000
taus = tau2 = 1000
lambda = sig2/tau2

# Iterate updates
for (it in 1:8) {
    G = BtB + lambda * P
    a = solve(G, Bty)
    mu = B %*% a
    r = y - mu
    H = solve(G, BtB)
    ed = sum(diag(H))
    sig2 = sum(r^2)/(m - ed - pord)
    tau2 = sum((D %*% a)^2)/(ed - pord)
    lanew = sig2/tau2
    dla = (lanew - lambda)/lambda
    lambda = lanew
    cat(it, ed, dla, "\n")
    sigs = c(sigs, sig2)
    taus = c(taus, tau2)
}

lsig = log(sqrt(sigs))
ltau = log(sqrt(taus))
xg = seq(min(x), max(x), length = 200)
Bg = bbase(xg, bdeg = 3, nseg = nseg)
yg = Bg %*% a

cols = c('red', 'blue')
F1 = data.frame(x = rep(1:nit, 2),
                y = log10(abs(c(diff(lsig), diff(ltau)))),
                Var = rep(c('sigma', 'tau'), each = nit))
plt1 = ggplot(F1, aes(x = x, y = y, group = Var, colour = Var, shape = Var)) +
  geom_point(data = F1, size = 2.5) +
  xlab("Iteration") + ylab("log10(relative change)") +
  ggtitle('Convergence of the HFS algorithm') +
  JOPS_theme()

# Plot and save
grid.arrange(plt1, ncol = 1, nrow = 1)

