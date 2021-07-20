# Penalized composite link model to estimate density of log transformed data (RDW data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(JOPS)
library(ggplot2)
library(gridExtra)

# Read the data
data(rdw)
set.seed(123)
y = sample(rdw, 2000)

# Compute histogram ----
ylo = 10
dy = 0.1
brk = seq(ylo, 20, by = dy)
h0 = hist(y, breaks = brk, plot = F)
u = h0$counts
nu = length(u)

# Set up C matrix ----
dx = 0.004
xlo = 1
xg = seq(xlo, 1.3, by = dx)
nx = length(xg)
C = matrix(0, nu, nx)
for (j in 1:nx) {
    xr = xlo + j * dx
    xl = xr - dx
    ul = 10^xl
    ur = 10^xr
    for (i in 1:nu) {
        yr = ylo + i * dy
        yl = yr - dy
        bl = max(ul, yl)
        br = min(ur, yr)
        if (br > bl) {
            v = (log10(br) - log10(bl))/(xr - xl)
            C[i, j] = v
        }
    }
}

# Prepare basis (the identity matrix here) and penalty
E = diag(nx)
D = diff(E, diff = 3)

# Places to store AIC and fits
aics = NULL
Gam = NULL

# Run over a range of (log-)lambdas
llas = seq(3, 5, by = 0.1)
eta = rep(log(sum(u)/nx), nx)
for (lla in llas) {
    lambda = 10^lla
    P = lambda * t(D) %*% D
    # Do the fiting
    for (it in 1:20) {
        gam = exp(eta)
        mu = c(C %*% gam)
        V = outer(1/mu, gam) * C
        M = diag(mu)
        G = t(V) %*% M %*% V
        enew = solve(G + P, t(V) %*% (u - mu) + G %*% eta)
        # Check convergence
        de = max(abs(enew - eta))
        # cat(it, de, '\n')
        eta = c(enew)
        if (de < 1e-04)
            break
    }
    # Compute AIC
    K = solve(G + P, G)
    ed = sum(diag(K))
    ok = u > 0
    dev = 2 * sum(u[ok] * log(u[ok]/mu[ok]))
    aic = dev + 2 * ed
    # Save results
    aics = c(aics, aic)
    Gam = cbind(Gam, gam)
    cat(lla, it, aic, "\n")
}

# Pick best lambda
k = which.min(aics)
gam = Gam[, k]

# Make the graphs
DF1 = data.frame(y = y)
DF2 = data.frame(x = h0$mids, y = mu)
DF3 = data.frame(x = xg, y = gam)
cfill = "wheat3"
cline = "steelblue"

plt1 = ggplot(DF1, aes(y)) +
  geom_histogram(binwidth = 0.1, fill = cfill, col = 'white') +
  ggtitle('Histogram and smooth estimate') +
  xlab('RDW') + ylab('Count') +
  geom_line(data = DF2, aes(x = x, y = y), col = cline, size = 1) +
  JOPS_theme()

plt2 = ggplot(DF3, aes(x = x, y = y)) +
  ggtitle('Smooth estimate on log scale') +
  geom_col(aes(x = x, y  = y), fill = cfill, col = 'white')+
#  geom_line(data = DF3, aes(x = x, y = y), col = cline, size = 1) +
  geom_hline(yintercept = 0) +
  xlim(1, 1.25)+
  xlab('log10(RDW)') + ylab('Count') +
  JOPS_theme()

# Show and save
grid.arrange(plt1, plt2, nrow = 2, ncol = 1)
