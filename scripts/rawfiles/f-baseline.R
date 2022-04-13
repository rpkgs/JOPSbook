# Estimation of baseline by expectile smoothing (indiumoxide data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

# library(diffractometry) [Currently no longer on CRAN]
library(ggplot2)
library(gridExtra)
library(JOPS)

# Get the data and make selection
data(indiumoxide)
x0 <- indiumoxide[,1]
y0 <- indiumoxide[,2]
sel = x0 > 47 & x0 < 65
x = x0[sel]
y = y0[sel]
m = length(x)

# Construct basis and penalty matrix
B = bbase(x, nseg = 40)
n = ncol(B)
D = diff(diag(n), diff = 2)
lambda = 1000
P = lambda * t(D) %*% D

# Estimate baseline
tau = 0.04
w = rep(0.5, m)
for (it in 1:10) {
  W = diag(c(w))
  a = solve(t(B) %*% W %*% B + P, t(B) %*% W %*% y)
  z = B %*% a
  wnew = ifelse(y > z, tau, 1 - tau)
  neq = sum(wnew == w)
  if (neq == m) break
  w = wnew
  cat(it, m, neq, '\n')
}

# Subtract baseline
r = y - z
Dat = data.frame(x = x, y = y, z = z, r = r)

plt1 = ggplot(Dat, aes(x = x, y = y) ) +
       geom_line(col = 'darkgrey', size = 0.2) +
       geom_line(aes(x = x, y = z), col = 'blue', size = 0.6) +
       JOPS_theme() +
       ggtitle('X-ray diffraction scan of indium oxide with baseline') +
       xlab('Diffraction angle') +
       ylab('Counts')

plt2 = ggplot(Dat, aes(x = x, y = r) ) +
       geom_line(col = 'darkgrey', size = 0.2) +
       JOPS_theme() +
       ggtitle('Baseline subtracted') +
       xlab('Diffraction angle') +
       ylab('Counts')

grid.arrange(plt1, plt2, ncol = 1, nrow = 2)



