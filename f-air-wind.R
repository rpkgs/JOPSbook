# New York air quality data polynomial fits
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(JOPS)

# Get the data
data(airquality)
Dat = data.frame(x = airquality$Wind, y = airquality$Ozone)
Dat = na.omit(Dat)

# Generate the graph
pl = qplot(x, y, data = Dat)  +
  geom_smooth(method = 'lm', formula = 'y ~ poly(x, 2)',
              colour = 'red', se = F) +
  geom_smooth(method = 'lm', formula = 'y ~ x',
              colour = 'blue', lty = 2, se = F) +
  xlab('Wind speed (mph)') + ylab('Ozone concentration (ppb)') +
  ggtitle('New York air quality') +
  JOPS_theme()

# Plot graph and save pdf
print(pl)

