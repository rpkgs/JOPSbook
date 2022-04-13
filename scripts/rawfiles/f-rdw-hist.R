# Histogram of red blood cell distribution width (linear and log scales) (RDW data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)

# Make the histogram plots
data(rdw)
mycol = "wheat3"
DF = data.frame(RDW = rdw)
cfill = 'wheat3'
plt1 = ggplot(DF, aes(x = RDW)) + ylab("Count") +
  geom_histogram(binwidth = 0.1, fill = cfill, color = 'white') +
  ggtitle('Histogram of observations') +
  JOPS_theme()

plt2 = ggplot(DF, aes(x = log10(RDW))) + ylab("Count") +
  geom_histogram(binwidth = 0.005, fill = cfill, color = 'white') +
  ggtitle('Histogram of their logarithms') +
  xlim(1, 1.25) +
  JOPS_theme()

# Make and save plot
grid.arrange(plt1, plt2, nrow = 2, ncol = 1)
