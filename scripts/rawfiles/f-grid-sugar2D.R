# Image of fluorescence measurements (Sugar data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(JOPS)
library(ggplot2)
library(reshape2)

# Define indices and dimensions
ex = c(340, 325, 305, 290, 255, 240, 230)
em = seq(from = 275, to = 560, by = 0.5)
p1 = length(ex)
p2 = length(em)
z = Sugar$X[1, ]

DF = data.frame(x = rep(em, p1), y = rep(ex, each = p2), z = z)

plt = ggplot(DF,  aes(x, y, fill = z)) +
  geom_tile(show.legend = T) +
  scale_fill_gradientn(colours = terrain.colors(100)) +
  #scale_fill_gradient(high = 'darkgreen', low = 'white') +
  xlab('Emission wavelength (nm)') +
  ylab('Excitation wavelength (nm)') +
  labs(fill = "Fluorescence") +
  JOPS_theme()

# Make and save pdf
plot(plt)
