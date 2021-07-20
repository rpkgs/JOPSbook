# Example of two-dimensonal image regressors (Sugar data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(fields)
library(JOPS)

# Get the data
x = Sugar$X
y = as.vector(Sugar$y[, 3])

# Center the image
x0 = x - apply(x, 1, mean)

# Define indices and dimensions
M1_index = rev(c(340, 325, 305, 290, 255, 240, 230))
M2_index = seq(from = 275, to = 560, by = 0.5)
p1 = length(M1_index)
p2 = length(M2_index)

# Get the image associated with obs 1
X1 = matrix(x0[1, ], p1, p2, byrow = T)

par(mfrow = c(1, 1))
image.plot(M2_index, M1_index, t(X1), col = terrain.colors(100), ylab = "Excitation wavelength (nm)",
    xlab = "Emission wavelength (nm)", cex=1.5)


