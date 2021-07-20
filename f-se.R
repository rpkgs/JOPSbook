# P-spline fit with twice se bands, optimal on CV (Motorcyle data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(colorspace)
library(MASS)
library(JOPS)

# Get the data
data(mcycle)
x = mcycle$times
y = mcycle$accel
Data = data.frame(x, y)

# Fitting, the tuning parameter chosen based on min CV
fit = psNormal(x, y, nseg = 20, bdeg = 3, pord = 2, lambda = 0.8)

# Twice se bands
selow=fit$ygrid - 2*fit$se_eta
seup=fit$ygrid + 2*fit$se_eta

# Dataframe for ggplot
G = data.frame(x = fit$xgrid, y = fit$ygrid, selow = selow, seup=seup)

# Make and save plot
qp = qplot(x, y, color = I('grey50')) +
  xlab('Time (ms)') +
  ylab('Acceleration (g)') +
  ggtitle('P-spline fit to motorcycle helmet data') +
  geom_hline(yintercept = 0) +
  geom_line(data = G, aes(x, y), color = col_mean, linetype = lty_mean, size = 1) +
  geom_line(data = G, aes(x, selow), color = col_se, size = 0.8, linetype = lty_se) +
  geom_line(data = G, aes(x, seup), color = col_se, size = 0.8, linetype = lty_se)  +
  JOPS_theme()

# Plot and save
print(qp)
