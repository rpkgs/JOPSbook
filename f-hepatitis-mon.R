# Monotonic smoothing by asymmetric penalty (Hepatitis Bulgaria data)
# Data source: Keiding, JRSSA, 1991
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(JOPS)

# Get the data
data(Hepatitis)
x = Hepatitis$Age
y = Hepatitis$Infected
v = Hepatitis$Sampled
m = length(y)

# Preparations
D = diff(diag(m), diff = 2)
D1 = diff(diag(m))
lambda = 100
P = lambda * t(D) %*% D
kappa = 1e5

# Fit without constraint
p = (y + 1) / (v + 2)
eta = log(p / (1 - p))
u = rep(0, m - 1)
for (it in 1:10) {
   p = 1 / (1 + exp(-eta))
   mu = v * p
   w = mu * (1 - p)
   W = diag(w)
   enew = solve(W + P, y - mu + w * eta)
   dz = max(abs(enew - eta))
   eta = enew
   u = diff(eta) < 0
   if (dz < 1e-4) break
}
p1 = p

# Fit with constraint
p = (y + 1) / (v + 2)
eta = log(p / (1 - p))
u = rep(0, m - 1)
for (it in 1:10) {
   p = 1 / (1 + exp(-eta))
   mu = v * p
   w = mu * (1 - p)
   W = diag(w)
   Q = t(D1) %*% diag(u) %*% D1
   enew = solve(W + P + kappa * Q, y - mu + w * eta)
   dz = max(abs(enew - eta))
   eta = enew
   u = diff(eta) < 0
   if (dz < 1e-4) break
}
p2 = p

# Plot results
DF1 = data.frame(x = x, ratio = y / v, p1 = p1, p2 = p2)
plt1 = ggplot(data = DF1) +
  geom_line(aes(x = x, y = ratio, color = 'black')) +
  geom_line(aes(x = x, y = p1, color = 'red'), size = 1) +
  geom_line(aes(x = x, y = p2, color = 'blue'), size = 1) +
  ggtitle('Hepatitis B prevalence') +
  xlab('Age') + ylab('Fraction')  +
  scale_color_identity(name = "",
                          breaks = c("black", "red", "blue"),
                          labels = c("Raw fraction", "Smooth", "Monotone smooth"),
                          guide = "legend") +
  theme(
    plot.title = element_text(size=16),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16)) +
  JOPS_theme() +
  theme(legend.position="bottom")

# Make and save pdf
plot(plt1)

