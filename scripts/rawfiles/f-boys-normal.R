# GAMLSS smoothing of weights of 100 boys (Dutch boys data)
# A graph in the book "Practical Smoothing. The Joys of P-splines"
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(colorspace)
library(AGD)
library(gamlss)
library(JOPS)

# Get the data
Dat = subset(boys7482, age > 5, select = c(age, wgt))
Dat = na.omit(Dat)
m0 = nrow(Dat)
m = 1000
set.seed(2013)
sel = sample(1:m0, m)
x = (Dat$age[sel])
y = Dat$wgt[sel]
Data = data.frame(x = x, y = y)

# Fit GAMLSS
g = gamlss(y ~ pb(x), sigma.formula = ~ pb(x), data = Data)

# Compute centile curves on grid
xg = seq(5, 22, by = 0.1)
pp = seq(0.05, 0.95, by = 0.1)
np = length(pp)
cnt = 100 * pp
Z = centiles.pred(g, xname = "x", xvalues = xg, type = "centiles", cent = cnt)

# Data frames for ggplot
nz = ncol(Z) - 1
Zmat = as.matrix(Z[, -1])
fcnt = as.factor(rep(cnt, each = length(xg) ))
Zmu = data.frame(x = rep(Z[, 1], nz), z = as.vector(Zmat), Centile = fcnt)
Data = data.frame(x, y)

# Build the graph
cols = rainbow_hcl(np, start = 20, end = 340)
plt1 = ggplot(Data, aes(x = x, y = y)) +
  geom_point(size = 1, color = "darkgrey") +
  geom_line(data = Zmu, aes(x = x, y = z, group = Centile, color = Centile), size =1)+
  ggtitle("GAMLSS (Normal) fit for Dutch boys") +
  xlab("Age (yr)") + ylab("Weight (kg)") +
  JOPS_theme() +
  scale_color_manual(values=cols)

# Plot and save
print(plt1)
