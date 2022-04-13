# Quantile curves for body weight of 1000 boys (Dutch boys data)
# Paul Eilers and Brian Marx 2019
# A graph in "Practical Smoothing. The Joys of P-splines"

library(colorspace)
library(AGD)
library(ggplot2)
library(JOPS)

# Get the data
data(boys7482)
Dat = subset(boys7482, age > 5, select = c(age, wgt))
Dat = na.omit(Dat)
m0 = nrow(Dat)
m = 1000
set.seed(2013)
sel = sample(1:m0, m)
x = (Dat$age[sel])
y = Dat$wgt[sel]

# P-spline paremeters
xl = min(x)
xr = max(x)
nsegx = 20
bdeg = 3
pp = seq(0.05, 0.95, by = 0.1)
#pp = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95)
np = length(pp)

# Compute bases
B = bbase(x, xl, xr, nsegx, bdeg)
nbx = ncol(B)
xg = seq(xl, xr, length = 100)
Bg = bbase(xg, xl, xr, nsegx, bdeg)

# Compute penalty matrices
D = diff(diag(nbx), diff = 2)
lambdax = 1
P = lambdax * t(D) %*% D

Zg = NULL
for (p in pp) {
  b = 0.001
  dzmax = 1e-3 * max(y)
  z = 0
  for (it in 1:500) {
    r = y - z
    w = ifelse(r > 0, p, 1- p)  / sqrt(b ^ 2 + r ^ 2)
    W = diag(c(w))
    Q = t(B) %*% W %*% B
    a = solve(Q + P, t(B) %*% (w * y))
    znew = B %*% a
    dz = sum(abs(z - znew))
#    cat(it, dz, '\n')
    if (dz < dzmax) break
    z = znew
  }
  Zg = cbind(Zg, Bg %*% a)
  cat(p, '\n')
}

# Prepare data for plots
DF1 = data.frame(x = x, y = y)
np = length(pp)
ng = length(xg)
DF2 = data.frame(x = rep(xg, np), y = as.vector(Zg), tau = as.factor(rep(pp, each = ng)))

plt1 =  ggplot(DF1,  aes(x = x, y = y)) +
  geom_point(color = I('grey75'), size = 0.9) +
  xlab('Age (yr)') + ylab('Weight (kg)') +
  ggtitle('Quantile curves for weights of 1000 Dutch boys') +
  geom_line(data = DF2, aes(x = x, y = y, group = tau, color = tau), size = 1) +
  JOPS_theme()  +
  labs(color =  bquote(tau) ) +
  scale_color_manual(values=rainbow_hcl(np, start = 10, end = 350))

plot(plt1)


