# External prediction performance plot (Biscuit data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(JOPS)

# Get the data
library(fds)
data(nirc)

iindex = nirc$x
X = nirc$y
sel = 50:650  #1200 <= x & x<= 2400
X = X[sel, ]
iindex = iindex[sel]
dX = diff(X)
diindex = iindex[-1]
y = as.vector(labc[1, 1:40])  #fat
oout = 23  # Denoted as outlier by Brown et al. (2001)
dX = t(dX[, -oout])
y = y[-oout]

# Optimal fit
fit1 = psSignal(y, dX, diindex, nseg = 25, lambda = 3.2e-07, pord = 3)

# Dataframes and ggplot
F1=data.frame(x=fit1$press_mu, y=y)
plt1 = ggplot(F1, aes(x = x, y = y)) +
  geom_point(size = 1.2) +
  geom_abline(intercept = 0, slope=1, size = 0.3) +
  xlab("Predicted %fat") +
  ylab("Observed %fat") +
  theme(
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16)) +
  ggtitle(" ") +
  JOPS_theme()

# Make and save pdf
print(plt1)
