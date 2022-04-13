# Smoothing of the circular frequency distribution of wind
# directions A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(JOPS)
library(circular)
library(ggplot2)

# Get the data and xpress as degrees (0 is North)
data(wind)
wind = wind * 360/(2 * pi)
wind = wind - 360 * (wind > 180)

# Make the histogram
brks = seq(-180, 180, by = 5)
hst = hist(wind, breaks = brks, plot = F)
x = hst$mids
y = hst$counts
m = length(y)

# Prepare basis and penalty
B0 = cbase(x, -180, 180, nseg = 20)
B = cbind(B0, rep(1, m))
n = ncol(B0)
D = cdiff(n)
P = matrix(0, n + 1, n + 1)
P[1:n, 1:n] = t(D) %*% D

# Set a range of log-lambdas
llas = seq(-3, 6, by = 0.1)
nla = length(llas)
aics = NULL
Mus = NULL

# Fit for all lambdas
eta = log(y + 1)
a = 0
for (lla in llas) {
    lambda = 10^lla
    # Fit smooth distribution
    for (it in 1:20) {
        mu = exp(eta)
        W = diag(c(mu))
        r = y - mu + mu * eta
        G = t(B) %*% W %*% B
        anew = solve(G + lambda * P, t(B) %*% r)
        da = max(abs(anew - a))
#        cat(it, da, "\n")
        if (da < 1e-04)
            break
        a = anew
        eta = B %*% a
    }
    # Compute AIC
    H = solve(G + lambda * P, G)
    ed = sum(diag(H))
    ok = y > 0
    aic = 2 * sum(y[ok] * log(y[ok]/mu[ok])) + 2 * ed
    aics = c(aics, aic)
    Mus = cbind(Mus, mu)
}

# Save last mu
mum = mu

# Select lowest AIC
k = which.min(aics)
cat("Log-lambda = ", llas[k], "AIC = ", aics[k], "\n")

Data = data.frame(x = x, y = y, muk = Mus[, k], mum = mum)
Data2= data.frame(x=wind)


plt2=ggplot(Data2, aes(x = x))+
  geom_histogram(fill = 'wheat3', binwidth = 4)+
  geom_hline(yintercept = 0) +
  xlab("Angle (degrees)") + ylab("Frequency") +
  ggtitle("Wind directions") +
  geom_line(data = Data, aes(x = x, y = muk), col = 'blue', size = 1) +
  geom_line(data = Data, aes(x = x, y = mum), col = 'red', size = 0.7) +
  theme(
    plot.title = element_text(size=16),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16)) +
  JOPS_theme()

# Make and save graph
print(plt2)



