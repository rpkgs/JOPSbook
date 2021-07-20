# Original and differenced NIR spectra (Biscuit data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(colorspace)
library(JOPS)

# Get the data
library(fds)
data(nirc)
x = nirc$x
Y = nirc$y
sel = 1200 <= x & x <= 2400
Y = Y[sel, ]
x = x[sel]
Yd = diff(Y)
dx = x[-1]
nr = ncol(Y)
m = length(x)

# Make data frames for ggplot
id =  as.factor(rep(1:nr, each = m))
F1 = data.frame(x = rep(x, nr), y = as.vector(Y), id = id)
idd = as.factor(rep(1:nr, each = m - 1))
F2 = data.frame(x = rep((dx), nr), y = as.vector(Yd), id = idd)
plt1 = ggplot()  +
  geom_line(aes(x = x, y = y,  col = id), data = F1, size = .1) +
  ggtitle("Observed spectra") +
  xlab("Wavelength (nm)") + ylab("-log(reflectance)") +
  JOPS_theme() +
  theme(legend.position = "none") +
  scale_color_manual(values = rainbow_hcl(40))

plt2 = ggplot()  +
  geom_line(aes(x = x, y = y, col = id), data = F2, size = .1) +
  ggtitle("Differences of spectra") + xlab("Wavelength (nm)") + ylab("") +
  JOPS_theme() +
  theme(legend.position = "none") +
  scale_color_manual(values = rainbow_hcl(40))


# Make and save the graphs
grid.arrange(plt1, plt2, ncol = 1, nrow = 2)
