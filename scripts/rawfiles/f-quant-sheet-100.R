# Quantile sheet for body weight of 100 boys (Dutch boys data)
# Paul Eilers and Brian Marx 2019
# A graph in "Practical Smoothing. The Joys of P-splines"

library(rgl)
library(AGD)
library(colorspace)
library(JOPS)
library(ggplot2)

# Get the datadata(boys7482)
Dat = subset(boys7482, age > 5, select = c(age, wgt))
Dat = na.omit(Dat)
m0 = nrow(Dat)
m = 100
set.seed(2013)
sel = sample(1:m0, m)
x = (Dat$age[sel])
y = Dat$wgt[sel]

# P-spline paremeters
xl = min(x)
xr = max(x)
nsegx = 10
nsegp = 10
bdeg = 3
pp = seq(0.05, 0.95, by = 0.1)
np = length(pp)

# Compute bases
Bx = bbase(x, xl, xr, nsegx, bdeg)
Bp = bbase(pp, 0, 1, nsegp, bdeg)
nbx = ncol(Bx)
nbp = ncol(Bp)
Tx = rowtens(Bx)
Tp = rowtens(Bp)

# Compute penalty matrices
Dx = diff(diag(nbx), diff = 2)
Dp = diff(diag(nbp), diff = 2)
lambdax = 0.3
lambdap = 1
Px = lambdax * t(Dx) %*% Dx
Pp = lambdap * t(Dp) %*% Dp
P = kronecker(Pp, diag(nbx)) + kronecker(diag(nbp), Px)
kappa = 0
P = P + kappa * diag(nrow(P))

# Initialize
Y = outer(y, rep(1, np))
Z = 0 * Y + mean(Y)
OP = outer(rep(1, m), pp)

# Iterate
b = 0.001
dzmax = 1e-3 * max(y)
for (it in 1:500) {
  R = Y - Z
  W = ifelse(R > 0, OP, 1- OP) / sqrt(b ^ 2 + R ^ 2)
  Q = t(Tx) %*% W %*% Tp
  dim(Q) = c(nbx, nbx, nbp, nbp)
  Q = aperm(Q, c(1, 3, 2, 4))
  dim(Q) = c(nbx * nbp, nbx * nbp)
  r = t(Bx) %*% (Y * W) %*% Bp
  dim(r) = c(nbx * nbp, 1)
  A = solve(Q + P, r)
  dim(A) = c(nbx, nbp)
  Znew = Bx %*% A %*% t(Bp)
  dz = sum(abs(Z - Znew))
  #cat(it, dz, '\n')
  if (dz < dzmax) break
  Z = Znew
}

xg = seq(xl, xr, length = 100)
Bg = bbase(xg, xl, xr, nsegx, bdeg)
Zg = Bg %*% A %*% t(Bp)

np = length(pp)
cols = diverge_hcl(np, h = c(240, -30), c= 100, l = 60, power = 1)

# Dataframes for ggplot
DF1 = data.frame(x = x, y = y)
np = length(pp)
ng = length(xg)
DF2 = data.frame(x = rep(xg, np), y = as.vector(Zg), tau = as.factor(rep(pp, each = ng)))

# Make and save plots
plt1 =  ggplot(DF1,  aes(x = x, y = y)) +
  geom_point(color = I('grey70'), size = 1.5) +
  xlab('Age (yr)') + ylab('Weight (kg)') +
  ggtitle('Quantile sheet for weights of 100 Dutch boys') +
  geom_line(data = DF2, aes(x = x, y = y, group = tau, color = tau), size = 1) +
  JOPS_theme()  +
  scale_color_manual(values=rainbow_hcl(np, start = 10, end = 350))

plot(plt1)
