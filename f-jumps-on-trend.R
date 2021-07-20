# Jumps on trend, smooth and piecewise-constant (simulated data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(spam)
library(ggplot2)
library(JOPS)
library(gridExtra)

# Make steps
m = 300
set.seed(2013)
u = rep(0, m)
v = cumsum(rpois(m, 10) + 10)
v = v[v < m]
n = length(v)
u[v] = rnorm(n)
u = cumsum(u)
u = u - mean(u) - 0.5

# Make trend and add up with noise
x = 1:m
g = 1 * sin(2 * pi* x / m)
y = g + u + rnorm(m) * 0.2

# Basis matrix and inproducts
E = diag.spam(m)
B = cbind(E, E)
BtB = t(B) %*% B
Bty = t(B) %*% y

# Penalty matrices
D1 = diff(E)
D2 = diff(D1)
lambda_1 = 10000
kappa = 1e-6
P = diag.spam(2 *m) * kappa
P2 = lambda_1 * t(D2) %*% D2
ip = (1:m) + m
P[ip, ip] = P2 + P[ip, ip]
lambda_2 = 0.1
beta = 0.001

# Run the iterations
s = 0
W = diag.spam(m - 1)
for (it in 1: 20) {

  # Update penalty and solve
  P0 = lambda_2 * t(D1) %*% W %*% D1
  P[1:m, 1:m] = P0
  snew = solve(BtB + P, Bty)
  ds = max(abs(s - snew))
  cat(it, ds, '\n')
  s = snew
  if (ds < 1e-4) break

  # New weights from jumps
  q = s[1:m]
  dq = diff(q)
  w = c(1 / (dq ^ 2 + beta ^ 2))
  W = diag.spam(w)
}
z = s[ip]

DF = data.frame(x = x, y = y, trend = z, jumps = q)

plt1 = ggplot(data = DF) +
  geom_line(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = trend), col = 'red', size = 1.2,
    linetype = "dashed") +
  geom_line(aes(x = x, y = trend), col = 'red', size = 0.7) +
  geom_line(aes(x = x, y = trend + jumps), col = 'blue', size = 1) +
  JOPS_theme() + xlab("") + ylab("") +
  ggtitle('Data (simulated), trend and fit')

plt2 = ggplot(data = DF) +
  geom_line(aes(x = x, y = y - trend)) +
  geom_line(aes(x = x, y = jumps), col = 'blue', size = 1) +
  JOPS_theme() + xlab("") + ylab("") +
  ggtitle('Detrended with fitted jumps')


plts = grid.arrange(plt1, plt2, nrow = 2, ncol = 1)
plot(plts)

