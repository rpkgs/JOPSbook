# Smoothing of height and its derivative against age (boys7482 data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(AGD)
library(JOPS)

# Get the data
data(boys7482)
agemax = 20
Data = subset(boys7482, !is.na(age) & !is.na(hgt) & age < agemax & age > 0)
x = Data$age
y = Data$hgt
m0 = length(y)
set.seed(2019)
sel = sample(1:m0, 1000)
x = x[sel]
y = y[sel]

# Set P-spline parameters
nseg = 50
xmin = 0
xmax = 20
B = bbase(x, xmin, xmax, nseg = nseg, bdeg = 3)
n = ncol(B)
D = diff(diag(n), diff = 2)
BtB = t(B) %*% B
Bty = t(B) %*% y
m = length(y)

# Smooth
lambda = 100
a = solve(BtB + lambda * t(D) %*% D , Bty)

# Compute trend on grid and its derivative
xg = seq(xmin, xmax, by = 0.1)
Bg = bbase(xg, xmin, xmax, nseg = nseg)
h = (xmax - xmin) / nseg
B1 = bbase(xg, xl = xmin, xr = xmax, nseg = nseg, bdeg = 2)
B1 = B1 %*% diff(diag(n)) / h
z = Bg %*% a
g = B1 %*% a

# Create data frame for ggplot
XY = data.frame(x = x, y = y)
Z = data.frame(x = xg, z = z, g = g )

plt1 =
  ggplot(XY, aes(x = x, y = y)) +
  geom_point(color = I("darkgrey"), size = 1) +
  xlab("Age") + ylab("Height (cm)") +
  xlim(0, 20) +
  geom_line(data = Z, aes(x = x, y = z), size = 1,
            color = I('blue')) +
  ggtitle('Heights of Dutch boys') +
  JOPS_theme()

plt2 =
  ggplot(Z, aes(x = x, y = g)) +
  xlab("Age") + ylab("Growth speed (cm/y)") +
  ylim(0, 25) +  xlim(0, 20) +
  geom_line(data = Z, aes(x = x, y = g), size = 1,
            color = I('blue')) +
  ggtitle('Growth speed of Dutch boys') +
  JOPS_theme()

grid.arrange(plt1, plt2, nrow = 2, ncol = 1)

