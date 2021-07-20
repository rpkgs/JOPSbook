# Construction of quadratic B-splines from truncated power basis
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)

m = 200
u = seq(0, 1, length = m)

# Compute the trunctaed quadratic functions
knt = seq(0.2, 0.8, by = 0.2)
n = length(knt)
U = outer(u, rep(1, n))
K = outer(rep(1, m), knt)
p = 2
P = (U - K) ^ p * (U > K)

# Partial sums
f4 = P[, 1]
f3 = P[, 1] - 3 * P[, 2]
f2 = P[, 1] - 3 * P[, 2] + 3 * P[, 3]
f1 = P[, 1] - 3 * P[, 2] + 3 * P[, 3] - P[, 4]
ind = as.factor(rep(1:4, each = m))
R = data.frame(x = u, y = c(f1, f2, f3, f4), ind = ind)
Q = data.frame(x = rep(u, 4), y = as.vector(P[, 4:1]), ind = ind)
S = data.frame(x = u, y = f1)
# Make the plots
plt1 = ggplot(Q, aes(x = x, y = y,  color = ind)) +
  ggtitle('Four truncated quadratic functions') +
  geom_line(size = 1)  +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_vline(xintercept = knt, col = 'orange', linetype = 2, size = 0.8) +
	xlim(0, 1) + ylim(0, 0.6)  +
	JOPS_theme() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "none")

plt2 = ggplot(R, aes(x = x, y = y, color = ind)) +
  ggtitle('Steps in the construction of one quadratic B-spline') +
  geom_line(data = S, aes(x = x, y = y), color = 'grey', size = 3) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_line(size = 1) +
  geom_vline(xintercept = knt, col = 'orange', linetype = 2, size = 0.8) +
  ylim(-0.1, 0.1)+   xlim(0, 1) +
  JOPS_theme() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "none")


# Show the plots and save pdf
grid.arrange(plt1, plt2, nrow = 2, ncol = 1)


