# P-spline fit using mixed model and fast Harville algorithm
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(MASS)
library(JOPS)

# Get the data
data(mcycle)
x = mcycle$times
y = mcycle$accel
m = length(y)
mmin = min(x)
mmax = max(x)

# Set P-spline parameters
nseg = 20
pord = 2
bdeg = 3

# Compute basis matrix and inner products
B = bbase(x, bdeg = bdeg, nseg = nseg)
n = ncol(B)
D = diff(diag(n), diff = 2)
P = t(D) %*% D
BtB = t(B) %*% B
Bty = t(B) %*% y
lambda = 1
for (it in 1:10) {
  G = BtB + lambda * P
  a = solve(G, Bty)
  mu = B %*% a
  r = y - mu
  H = solve(G, BtB)
  ed = sum(diag(H))
  sig2 = sum(r ^ 2) / (m - ed)
  tau2 = sum((D %*% a) ^ 2) / ed
  lanew = sig2 / tau2
  dla = (lanew - lambda) / lambda
  lambda = lanew
  cat(it, ed, dla, "\n")
}

xg = seq(min(x), max(x), length = 200)
Bg = bbase(xg, bdeg = 3, nseg = nseg)
yg = Bg %*% a

titl = bquote("HFS algorithm:" ~ lambda == .(round(lambda, 2)))

# Make the plot
F1 = data.frame(x,y)
F2 = data.frame(xg1 = xg, yg1 = yg)
plt1 = ggplot(F1, aes(x = x, y = y)) +
  geom_point(data = F1, size = 1.5, color = "darkgray") +
  geom_line(aes(x = xg1, y = yg1), data = F2, size = 1, colour = I("blue"), lty = 1) +
  xlab("Time (ms)") + ylab("Acceleration (g)") +
  ggtitle(titl) +
  JOPS_theme()

# Plots and save graphs
graphics.off()
grid.arrange(plt1, ncol = 1, nrow = 1)

