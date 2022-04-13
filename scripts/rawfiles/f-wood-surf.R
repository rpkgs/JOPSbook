# Optimal smoothing with L- and V-curve tuning (Wood surface data)
# A graph in the book 'The Joy of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(spam)
library(JOPS)

y = Woodsurf$y
m = length(y)
x = 1:m
w = 0 * y + 1

# Prepare for smoothing
E = diag.spam(m)
D = diff(E, diff = 2)

# Note: lambda is chosen from L- and V- curves
lambda = 2800000
E = diag(m)
D = diff(E, diff = 2)
P = lambda * t(D) %*% D
z = solve(E + P, y)
# y = Y[,1]
m = length(y)
x = 1:m
w = 0 * y + 1

# Prepare dataframes for ggplot
F1 = data.frame(x = x, y = y, z= z)
pl = ggplot(F1) +
  geom_point(aes(x = x, y = y), colour = 'darkgrey')  +
  geom_line(aes(x = x, y = z),  color = 'blue', size = 1, lty = 1) +
  geom_line(aes(x = x, y = y), colour = 'darkgrey') +
  xlab('Position') + ylab('Height (unknown units)') +
  ggtitle('Wood surface') +
  JOPS_theme()

# Plot and save pdf
print(pl)

