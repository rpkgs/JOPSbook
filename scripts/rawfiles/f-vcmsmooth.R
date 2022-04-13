# Plot of smoothly varying slope (Hard disk price data)
# A graph in the book "Practical Smoothing. The Joys of P-splines"
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(gamlss.data)
library(JOPS)

# Get the data
Dat1 = Disks
Month = Dat1$Month+12*(Dat1$Year==2000)
Size = Dat1$Size
Price = Dat1$PriceDG/2.2
y = Price
m = length(y)
x = Month

# Fit a few SLR by month
slope1 = rep(0, 6)
for (ii in 2:7) {
    slope1[ii - 1] = lsfit(Size[Month == ii], y[Month == ii])$coef[2]
}
kk = 0
slope2 = rep(0, 5)
for (ii in 9:13) {
    kk = kk + 1
    slope2[kk] = lsfit(Size[Month == ii], y[Month == ii])$coef[2]
}
slopeframe = data.frame(slopes = c(slope1, slope2), iindex = c(2:7, 9:13))

# Construct basis for VCM
B = bbase(x, nseg = 10)
n = ncol(B)
BtB = t(B) %*% B
D = diff(diag(n), diff = 2)
lambda = 100
P = t(D) %*% D
mon = F
xg = seq(2, 13, length = 100)
Bgrid = bbase(xg, nseg = 10)

# Scaled basis and additive structure
C1 = diag(Size) %*% B
C = cbind(B, C1)
lambda = c(100, 100)  # Values determined with grid search
Q = kronecker(diag(lambda), P)
K = diag(length(diag(Q)))
Q = Q + 1e-6 * K

# Solution to VCM
z = y
S = t(C) %*% C
a = solve(S + Q, t(C) %*% z)
mu = C %*% a

# Indices for coefficients
r1 = 1:n
r2 = r1 + n

G = solve(S + Q, S)
Br = solve(S + Q)
Brsub = Br[r2, r2]
ed = sum(diag(G))
trend = C %*% a
trendgrid = Bgrid %*% a[r2]
sigma = sqrt(sum((y - mu)^2) / (m - ed))
vara = Bgrid %*% Brsub %*% t(Bgrid)
sea = sigma * sqrt(diag(vara))
se2up = trendgrid + 2 * sea
se2lo = trendgrid - 2 * sea

# Set up plot
Fdat = data.frame(x, y, mu)
F2dat = data.frame(xg, trendgrid, se2up, se2lo)

plt1 = ggplot(slopeframe, aes(x = iindex, y = slopes)) +
    geom_point(data = slopeframe, size = 2.2) +
    scale_x_continuous(breaks = c(2:13),
      labels = c("Feb", "Mar", "Apr", "May", "Jun", "Jul",
                 "Aug", "Sep", "Oct", "Nov", "Dec", "Jan")) +
    geom_line(aes(x = xg, y = trendgrid), data = F2dat, size = 1,
      colour = col_mean, lty = 1) +
    geom_line(aes(x = xg, y = se2up), data = F2dat, size = 1,
      colour = col_se, lty = 2) +
    geom_line(aes(x = xg,  y = se2lo), data = F2dat, size = 1,
      colour = col_se, lty = 2) + xlab("Month (1999-2000)") +
    geom_hline(yintercept = 0, color = "darkgrey") +
    ylab("Slope: Euro/GB") +
    ggtitle('Trend of hard disk prices') +
    ylim(c(0, 30)) +
    JOPS_theme()

# Plot and save file
grid.arrange(plt1, ncol = 1, nrow = 1)
