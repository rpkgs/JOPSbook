# GAMLSS smoothing for percentile curves (Motorcycle data)
# A graph in the book "Practical Smoothing. The Joys of P-splines"
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(colorspace)
library(gamlss)
library(JOPS)
library(gridExtra)

# Get the data
data(mcycle)
x = mcycle$times
y = mcycle$accel
dat = data.frame(x =x, y = y)

# Fit GAMLSS
gx = gamlss(y ~ pb(x), sigma.formula = ~ pb(x), data = dat)
g1 = gamlss(y ~ pb(x), sigma.formula = ~ 1, data = dat)
g1 = gamlss(y ~ pb(x), data = dat)

# Compute centile curves on grid
xg = seq(2.5,  57.5, by = 0.1)
pp = c(0.05, 0.5, 0.95)
np = length(pp)
cnt = 100 * pp
Zx = centiles.pred(gx, xname = "x", xvalues = xg, type = "centiles", cent = cnt)
Z1 = centiles.pred(g1, xname = "x", xvalues = xg, type = "centiles", cent = cnt)

# Data frames for ggplot
nz = ncol(Zx) - 1
Zmat = as.matrix(Zx[, -1])
fcnt = as.factor(rep(cnt, each = length(xg) ))
Zmu = data.frame(x = rep(Zx[, 1], nz), z = as.vector(Zmat), Centile = fcnt)
Data = data.frame(x, y)

# Build the graph
cols = c('red','blue', 'red')
plt1 = ggplot(Data, aes(x = x, y = y)) +
  geom_point(size = 1.2, color = I("black")) +
  geom_line(data = Zmu, aes(x = x, y = z, linetype = Centile != 50,
            color = Centile), size =0.7)+
  ggtitle("Varying standard deviation") +
  xlab("Time(ms)") + ylab("Acceleration (g)") +
  JOPS_theme() +
  theme(legend.position = "none")  +
  scale_color_manual(values = cols)

# Data frames for ggplot
nz = ncol(Z1) - 1
Zmat = as.matrix(Z1[, -1])
fcnt = as.factor(rep(cnt, each = length(xg) ))
Zmu = data.frame(x = rep(Z1[, 1], nz), z = as.vector(Zmat), Centile = fcnt)
Data = data.frame(x, y)

plt2 = ggplot(Data, aes(x = x, y = y)) +
  geom_point(size = 1.2, color = I("black")) +
  geom_line(data = Zmu, aes(x = x, y = z, linetype = Centile != 50,
            color = Centile), size =0.7)+
  ggtitle("Constant standard deviation") +
  xlab("Time(ms)") + ylab("Acceleration (g)") +
  JOPS_theme() +
  theme(legend.position = "none")  +
  scale_color_manual(values = cols)

# Plot and save
grid.arrange(plt1, plt2, nrow = 2, ncol = 1)
