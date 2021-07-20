# Play with mixed model for P-splines (Motorcycle data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(JOPS)
library(nlme)
library(MASS)
library(ggplot2)

makeZ <- function(x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3) {
# construct smoothing random effect
  B <- bbase(x, xl, xr, nseg)
  D <- diff(diag(ncol(B)), diff = 2)
  Q <- solve(D %*% t(D), D)
  Z <- B %*% t(Q)
  Z
}

# Get the data
data(mcycle)
x = mcycle$times
y = mcycle$accel
grp = 0 * x + 1   # A 'grouping variable', needed for lme()

# Set P-spline parameters
xmin = 0
xmax = max(x)
nseg = 20
bdeg = 3
Z = makeZ(x, xmin, xmax, nseg, bdeg)

# Mixed model fitting with lme
lm1 <- lme(y ~ x, random = list(grp = pdIdent(~ Z - 1)))
cfix = lm1$coefficients$fixed
cran = c(lm1$coefficients$random$grp)

# Compute fit on grid
xg = seq(xmin, xmax, length = 100)
Zg = makeZ(xg, xmin, xmax, nseg, bdeg)
fitg = cfix[1] + cfix[2] * xg + Zg %*% cran

# Create data frames for ggplot
Data = data.frame(x = x, y = y)
Df1 = data.frame(x = xg, y = fitg, id = as.factor(1))

# Build the graph
plt1 = ggplot(Data, aes(x = x, y = y)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_line(data = Df1, size = 1, colour = I("blue")) +
  xlab("Time (ms)") + ylab("Acceleration (g)") +
  ggtitle("Motorcycle data with lme fit") +
  ylim(c(-150, 100)) +
  JOPS_theme()

print(plt1)
