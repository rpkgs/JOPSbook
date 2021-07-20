# Polynomial fits with differing support (Motorcycle data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(MASS)
library(ggplot2)
library(JOPS)

# Get the data
data(mcycle)
x = times = mcycle$times
y = accel = mcycle$accel
Data = data.frame(x, y)
showse = F
ltyp = c(1, 6)
lsize = 0.9

make_grid = function(x, n = 100) xg = seq(min(x), max(x), length = n)

# Predictions at different subsets of x
mc1 = subset(mcycle, times > 5)
lm1 = lm(accel ~ poly(times, 9), data = mc1)
xg1 = make_grid(mc1$times)
fit1 = predict(lm1, data.frame(times = xg1))
F1 = data.frame(times = xg1, accel = fit1)
mc2 = subset(mcycle, times > 0)
lm2 = lm(accel ~ poly(times, 9), data = mc2)
xg2 = make_grid(mc2$times)
fit2 = predict(lm2, data.frame(times = xg2))
F2 = data.frame(times = xg2, accel = fit2)

# Plot and save graph
qp = ggplot(mcycle, aes(x = times, y = accel)) +
  geom_point() +
  xlab('Time (ms)') +
  ylab('Acceleration (g)') +
  ggtitle('Polynomial fits to motor cycle helmet data') +
  geom_hline(yintercept = 0) +
  geom_line(data = F1, color = I("red"), size = 1, lty = 6) +
  geom_line(data = F2, color = I("blue"), size = 1, lty = 1) +
  JOPS_theme()

# Make plot and save pdf
print(qp)


