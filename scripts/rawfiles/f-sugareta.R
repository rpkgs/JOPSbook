# External prediction performance plot (Sugar data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(fields)
library(JOPS)
library(ggplot2)
library(gridExtra)

# Get the data
x0 = Sugar$X
x0 = x0 - apply(x0, 1, mean)  #center Signal
y = as.vector(Sugar$y[, 3])  #Response is Ash

# Inputs for two-dimensional signal regression
nseg = c(7, 37)
pord = c(3, 3)
min_ = c(230, 275)
max_ = c(340, 560)
M1_index = rev(c(340, 325, 305, 290, 255, 240, 230))
M2_index = seq(from = 275, to = 560, by = 0.5)
p1 = length(M1_index)
p2 = length(M2_index)
int <- T  # intercept in model

# Fit optimal model based on LOOCV, found via svcm
opt_lam=c(8858.6679, 428.1332)
Pars_opt=rbind(c(min_[1],max_[1],nseg[1],3,opt_lam[1],pord[1]),
               c(min_[2],max_[2],nseg[2],3,opt_lam[2],pord[2]))
fit=fit_opt=ps2DSignal(y, x0, p1, p2, "unfolded", M1_index, M2_index, Pars_opt,
                    int=int, ridge_adj=1e-6)

# Set up dataframes for ggplot
F1=data.frame(x=fit_opt$press_mu, y=y)
plt1 = ggplot(F1, aes(x = x, y = y)) +
  geom_point(size = 1.2) +
  geom_abline(intercept = 0,slope=1,size = 0.3) +
  xlab("Predicted Ash") + ylab("Observed Ash") +
  ggtitle(" ") +
  JOPS_theme()

# Make and save pdf
grid.arrange(plt1, ncol = 1, nrow = 1)
