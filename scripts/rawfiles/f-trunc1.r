# Construction of linear B-splines from truncated power basis
# A graph in the book "Practical Smoothing. The Joys of P-splines"
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)

# Compute the trunctaed linear functions
m = 200;
u = seq(0.0, 1, length = m)
knt = seq(0.2, 0.6, by  = 0.2)
n = length(knt);
U = outer(u, rep(1, n))
K = outer(rep(1, m ), knt);
P = (U - K) * (U > K);

# Partial sums
f0 = P[, 1]
f1 = P[, 1] - 2 * P[, 2];
f2 = P[, 1] - 2 * P[, 2] + 1 * P[, 3];

# Fill data frames for ggplot
ind = as.factor(rep(1:3, each = m))
ff1 = c(P[, 3], P[, 2], P[, 1])
V1 = data.frame(x = rep(u, 3), y = ff1, ind = ind)
ff2 = c(P[, 3], -2 *  P[, 2], P[, 1])
V2 = data.frame(x = rep(u, 3), y = ff2, ind = ind)
cols = c('blue', 'green', 'red')
V3 = data.frame(x = u, y = f2)

# Make the plots
plt1 = ggplot(V1, aes(x = x, y = y, group = ind, color = ind))+
  ggtitle('Three truncated linear functions') +
  geom_line(size = 1)  +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_vline(xintercept = knt, col = 'orange', linetype = 2, size = 0.8) +
	xlim(0, 1) + ylim(0, 0.6)  +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
	JOPS_theme() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "none")

plt2 = ggplot() + #V2, aes(x = x, y = y, group = ind2, color = ind2) ) +
  ggtitle('Scaled truncated linear functions, and their sum') +
  geom_line(data = V3, aes(x = x, y = y), color = 'grey', size = 3) +
	geom_hline(yintercept = 0) +
  geom_vline(xintercept = knt, col = 'orange', linetype = 2, size = 0.8) +
#  geom_point(aes(x = knots, y = 0 * knots), size = 3) +
	geom_line(data = V2, aes(x = x, y = y, group = ind, color = ind), size= 0.9) +
	ylim(-0.6, 0.6) + xlim(0, 1) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
	JOPS_theme() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "none")

# Show the plots and save pdf
grid.arrange(plt1, plt2, nrow = 2, ncol = 1)
