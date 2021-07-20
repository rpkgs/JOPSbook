# Binomial smoothing, optimal on AIC (Kyphosis data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)
library(rpart)

# Extract data
y = kyphosis$Kyphosis
x = kyphosis$Age
y = 1 * (kyphosis$Kyphosis == "present")  # make y 0/1

# Explore lambdas
llamb = seq(from = -2, to = 3, by = 0.25)
lambin = 10^llamb
r = length(llamb)
aics2 = aics3 = rep(0, r)
for (i in 1:r) {
    pden = psBinomial(x, y, nseg = 10, bdeg = 3, pord = 2, lambda = lambin[i],
        show = F)
    aics2[i] = pden$aic
    pden = psBinomial(x, y, nseg = 10, bdeg = 3, pord = 3, lambda = lambin[i],
        show = F)
    aics3[i] = pden$aic
}

# Plot AIC traces
F1 = data.frame(x = llamb, y2 = aics2, y3 = aics3)

plt1 = ggplot(F1, aes(x = x, y = aics2)) +
  geom_point(size = 1.5, color="blue") +
  xlab(expression(log10(lambda))) + ylab('AIC') +
  ggtitle('Second order penalty') +
  JOPS_theme()

plt2 = ggplot(F1, aes(x = x, y = aics3)) +
  geom_point(size = 1.5, color="blue") +
  xlab(expression(log10(lambda))) + ylab('AIC') +
  ggtitle('Third order penalty') +
  JOPS_theme()

k2 = which.min(aics2)
lambda2 = lambin[k2]
pden2 = psBinomial(x, y, nseg = 10, bdeg = 3, pord = 2, lambda = lambda2, show = F)
k3 = which.min(aics3)
lambda3 = lambin[k3]
pden3 = psBinomial(x, y, nseg = 10, bdeg = 3, pord = 3, lambda = lambda3, show = F)

F2 = data.frame(x = x, y = y)

plt3 = ggplot(F2, aes(x = x, y = y)) +
  geom_point(size = 1.5, color="blue") +
  xlab('Age (months)') + ylab('P(Kyphosis)') +
  ggtitle('Second order penalty') +
  JOPS_theme()

F3 = data.frame(x = pden2$xgrid, y2 = pden2$pgrid, y3 = pden3$pgrid)
plt3 = ggplot(F2, aes(x = x, y = y)) +
  geom_point(size = 1.5, color="blue") +
  geom_line(data = F3, color = 'red', size= 1, aes(x = x, y = y2)) +
  xlab('Age (months)') + ylab('P(Kyphosis)') +
  ggtitle('Second order penalty') +
  JOPS_theme()

plt4 = ggplot(F2, aes(x = x, y = y)) +
  geom_point(size = 1.5, color="blue") +
  geom_line(data = F3, color = 'red', size= 1, aes(x = x, y = y3)) +
  xlab('Age (months)') + ylab('P(Kyphosis)') +
  ggtitle('Third order penalty') +
  JOPS_theme()

# Save graphs
grid.arrange(plt1, plt2, plt3, plt4, ncol = 2, nrow = 2)
