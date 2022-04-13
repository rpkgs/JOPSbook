# Whittaker trend and sinus, signal separation (ECG data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

# Read Matlab data from Physiobank
library(spam)
library(ggplot2)
library(gridExtra)
library(JOPS)

# Get the data
data(ECG)
X = ECG

m0 = nrow(X)
time = (1:m0) * 0.004
sel = which(0.3 < time & time < 1.3)
y = X[sel, 3]
time = time[sel]

{
  period <- (1/60) / 0.004
  r = whit_pen2(y, period, lamb_anom = 1, lamb_season = 1e3)
  par(mfrow = c(2, 1), mar = c(3, 3, 1, 1))
  plot(r$z1, type = "l")
  plot(r$z2, type = "l")
}

# Simple smooth
lambda3 = 1
z3 = solve(E + lambda3 * P1, y)
lambda4 = 10
z4 = solve(E + lambda4 * P1, y)

#Plot data and fits
lw = 0.3
lwd = 0.6
Data = data.frame(x = time, y = y)
Dfit = data.frame(x = time, y = z4)
titl =  bquote("Whittaker smoother" ~ lambda== 10)
plt1 = ggplot(Data, aes(x = x, y = y)) +
  geom_line(aes(x = x, y = y), size = lw, color = grey(0.60)) +
  geom_line(data = Dfit, aes(x = x, y = y), size = lwd, color = 'blue') +
  xlab('Time (s)') + ylab('ECG') +
  ggtitle(titl) +
  JOPS_theme()

Data = data.frame(x = time, y = y)
Dfit = data.frame(x = time, y = z3)
titl =  bquote("Whittaker smoother" ~ lambda== 1)
plt2 = ggplot(Data, aes(x = x, y = y)) +
  geom_line(aes(x = x, y = y), size = lw, color = grey(0.60)) +
  geom_line(data = Dfit, aes(x = x, y = y), size = lwd, color = 'blue') +
  xlab('Time (s)') + ylab('ECG') +
  ggtitle(titl) +
  JOPS_theme()

Data = data.frame(x = time, y = y)
Dfit = data.frame(x = time, y = z1)
titl =  bquote("Signal estimate" ~ lambda== .(lambda1))
plt3 = ggplot(Data, aes(x = x, y = y)) +
  geom_line(aes(x = x, y = y), size = lw, color = grey(0.60)) +
  geom_line(data = Dfit, aes(x = x, y = y), size = lwd, color = 'blue') +
  xlab('Time (s)') + ylab('ECG') +
  ggtitle(titl) +
  JOPS_theme()

Data = data.frame(x = time, y = z2)
Dfit = data.frame(x = time, y = z4)
titl =  bquote("Interference" ~ lambda== .(lambda2))
plt4 = ggplot(Data, aes(x = x, y = y)) +
  geom_line(aes(x = x, y = y), size = lw, color = "blue") +
  xlab('Time (s)') + ylab('ECG') +
  ggtitle(titl) +
  ylim(500 * c(-1, 1)) +
  JOPS_theme()

# Create plots and pdf
grid.arrange(plt1, plt3, plt2, plt4, ncol = 2, nrow = 2)



