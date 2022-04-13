# Harmonic penalty smoothing (Variable star data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)

# Get the data
data(Varstar)
x = Varstar$V1
y = Varstar$V2
sel = y < 16 & (x < 3000 | x > 3300)
x = x[sel]
y = y[sel]
y = y - mean(y)

# B-spline basis
nseg = 100
xl = 2400
xr = 4200

B = bbase(x, xl, xr, nseg = nseg)
n = ncol(B)

# Construct penalty matrix
D = diff(diag(n), diff = 2)
dk = (xr - xl)/nseg
per = 110
phi = 2 * pi * dk/per
dp = c(1, -2 * cos(phi), 1)
Dp = matrix(0, n - 2, n)
for (j in 1:(n - 2)) Dp[j, j:(j + 2)] = dp
lambda = 0.1
Ph = lambda * t(Dp) %*% Dp
P = lambda * t(D) %*% D

a = solve(t(B) %*% B + P, t(B) %*% y)
ah = solve(t(B) %*% B + Ph, t(B) %*% y)
mu = B %*% ah

xg = seq(xl, xr, length = 500)
Bg = bbase(xg, xl, xr, nseg = nseg)
z = Bg %*% a
zh = Bg %*% ah

# Dataframes for ggplot
Fdat = data.frame(x = x, y = y)
Fz = data.frame(x = xg, z, zh)
lim = c(-1, 1) * 1.5

# Make and save plots
clr = 'tomato'
plt1 = ggplot(data = Fdat, aes(x = x, y= y)) +
  geom_point(color = clr, size = 0.9) +
  xlab("Day") + ylab("Centered magnitude") +
  ggtitle('Variable star, standard penalty, d = 2') +
  ylim(lim) +
  geom_line(aes(x = x, y = z), data = Fz, color = 'blue') +
  JOPS_theme()

plt2 = ggplot(data = Fdat, aes(x = x, y= y)) +
  geom_point(color = clr, size = 1) +
  xlab("Day") + ylab("Centered magnitude") +
  ggtitle('Variable star, harmonic penalty') +
  ylim(lim) +
  geom_line(aes(x = x, y = zh), data = Fz, color = 'blue') +
  JOPS_theme()

grid.arrange(plt1, plt2, nrow = 2, ncol = 1)
