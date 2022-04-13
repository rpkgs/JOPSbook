# Smoothing of BMI against age with zero slope constraint (Dutch boys data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(AGD)
library(JOPS)

# Get the data
data(boys7482)
agemax = 18
Data = subset(boys7482, !is.na(age) & !is.na(hgt) & age < agemax & age > 10)
x = Data$age
y = Data$hgt
m0 = length(y)
set.seed(2019)
sel = sample(1:m0, 500)
x = x[sel]
y = y[sel]

# Set P-spline parameters
nseg = 50
xmin = 10
xmax = 22
B = bbase(x, xmin, xmax, nseg = nseg, bdeg = 3)
n = ncol(B)
D = diff(diag(n), diff = 2)
BtB = t(B) %*% B
Bty = t(B) %*% y
m = length(y)

# Specify the subdomain where the slope should be zero
xc = seq(21, xmax, by = 0.2)
B1 = bbase(xc, xl = xmin, xr = xmax, nseg = nseg, bdeg = 2)
B1 = B1 %*% diff(diag(n))

# Smooth
lambda = 1000
kappa = 1e8
a = solve(BtB + lambda * t(D) %*% D + kappa * t(B1) %*% B1, Bty)
a0 = solve(BtB + lambda * t(D) %*% D , Bty)

# Compute trend on grid
xg = seq(xmin, xmax, by = 0.1)
Bg = bbase(xg, xmin, xmax, nseg = nseg)
z = Bg %*% a
z0 = Bg %*% a0

 # Create data frames for ggplot
XY = data.frame(x = x, y = y)
Z = data.frame(x = xg, z = z, z0 = z0 )

plt1 =
  ggplot(XY, aes(x = x, y = y)) +
  geom_point(color = I("darkgrey"), size = 1) +
  xlab("Age") + ylab("Height (cm)") +
  geom_line(data = Z, aes(x = x, y = z), size = 1,
            color = I('blue')) +
  geom_line(data = Z, aes(x = x, y = z0), size = 1,
            color = I('red'), lty = 2) +
  xlim(xmin, xmax) +  ylim(120, 205) +
  ggtitle('Heights of Dutch boys') +
  geom_vline(xintercept = min(xc), lty = 3, size = 1) +
  geom_vline(xintercept = max(xc), lty = 3, size = 1) +
  JOPS_theme()

# Make and save pdf
plot(plt1)
