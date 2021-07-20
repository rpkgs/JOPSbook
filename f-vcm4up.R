# Illustration of varying slope (Hard disk price data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)

Dat1 = Disks
Month = Dat1$Month+12*(Dat1$Year==2000)
Size = Dat1$Size
Price = Dat1$PriceDG/2.2

# Generate the plots
plts = list()
k = 0
slopes = c(16.8, 20.5, 14, 7.6)
months = c("February", "May", "September", "December")

for (j in c(2, 5, 9, 12)) {

    # Get the data
    k = k + 1
    Dat = na.omit(data.frame(Size, Price))
    names(Dat) = c("x", "y")
    x = Dat$x[Month == j]
    y = Dat$y[Month == j]

    # Create data frames for ggplot
    Data = data.frame(x, y)
    titl = paste(months[k], "Slope=", slopes[k])

    # Build the graph
    plt1 = ggplot(Data, aes(x = x, y = y)) +
      ylim(c(100, 500)) + xlim(c(5, 35)) +
      geom_point(color = "grey40") +
      geom_smooth(method = "lm", formula = "y~x", color = "red", se = FALSE) +
      xlab("Hard drive size (GB)") +
      ylab("Price (Euro)") +
      ggtitle(titl) +
      JOPS_theme() +
      theme(plot.title = element_text(size = rel(0.9)))

    # Add to the list of plots
    plts[[k]] = plt1
}

# Plot and save
grid.arrange(grobs = plts, nrow = 2, ncol = 2)
