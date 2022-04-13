# Non-adaptive smoothing with SOP (indiumoxide data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(SOPExamples) # (see https://bitbucket.org/mxrodriguez/sop)
library(JOPS)
library(ggplot2)
library(gridExtra)

# Load the indiumoxide data
data(indiumoxide)
sel <- 1:2000
x <- indiumoxide[sel, 1]
y <- indiumoxide[sel, 2]
n <- length(x)

# SOP method
bdeg = 3
pord = 2
fitpar = list(nseg = 200, bdeg = bdeg, pord = pord)
adpar = list(nseg = 80, bdeg = bdeg)
fit <- adaptive1D.SOP(x = x, y = y, smooth.cov = fitpar, adaptive = F,
         smooth.sp = adpar, trace = F, family = poisson(link = "log"))
predx <- data.frame(x = x)
predmu <- predict(fit, newdata = predx)

Pf = data.frame(x = x, y = y, mu = fit$mu)

plt1 = ggplot(Pf, aes(x = x, y = y) ) +
       geom_line(col = 'darkgrey', size = 0.2) +
       geom_line(aes(x = x, y = mu), col = 'blue', size = 0.6) +
       JOPS_theme() +
       ggtitle('X-ray diffraction of indium oxide with non-adaptive smooth') +
       xlab('Diffraction angle') +
       ylab('Counts')

plt2 = ggplot(Pf, aes(x = x,  y = y-mu) ) +
       geom_line(col = 'darkgrey', size = 0.2) +
       JOPS_theme() +
       ggtitle('Resiuals') +
       xlab('Diffraction angle') +
       ylab('Residuals')

grid.arrange(plt1, plt2, nrow = 2, ncol = 1)



