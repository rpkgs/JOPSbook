# PS-ANOVA for precipitation anomaly data. Plot with and without clipping (USA precip data)
# A graph from the book "Practical Smoothing. The Joys of P-splines"
# Paul Eilers and Brian Marx, 2019

library(fields)
library(data.table)
library(SpATS)
library(ggplot2)
library(JOPS)
library(gridExtra)
library(maps)
library(sp)
library(reshape2)

# Get USA outline
usa = map("usa", region = "main", plot = F)

# Get precipitation data
data(USprecip)
dats = data.frame(USprecip)

# Fit PS-ANOVA model
t0 = proc.time()
ctrl =  list(tolerance = 1e-05, monitoring = F, maxit = 100, update.psi = F)
fit = SpATS.nogeno(response = "anomaly",
						spatial = ~PSANOVA(lat, lon,  nseg = c(30,30), nest.div = c(2,2)),
            data = dats, control = ctrl)
summary(fit)
t1 = proc.time() - t0
cat(t1, '\n')

# Extract the surfaces
Tr = obtain.spatialtrend(fit, grid = c(200, 200))

# Turn matrix into a "long" data frame
Mu = Tr$fit
rownames(Mu) = Tr$row.p
colnames(Mu) = Tr$col.p
dens = melt(Mu)
names(dens) = c('x', 'y', 'z')

# Plot fit with contours
sl = T
ccol = 'blue'
capsize = 15
plt2 = ggplot(dens,  aes(x, y, fill = z)) +
       geom_raster(show.legend = sl) +
       scale_fill_gradientn(colours = terrain.colors(100))+
       xlab('Longitude') + ylab('Latitude') +
       ggtitle('Precipitation anomaly, smoothed') +
       geom_contour(aes(z = z), color = ccol, show.legend = T) +
       JOPS_theme() +
       coord_fixed() +
       theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             plot.title = element_text(size = capsize),
             axis.title.x = element_text(size = capsize),
             axis.title.y = element_text(size = capsize))


# Find points within US boundary for clipping
v = point.in.polygon(dens$x, dens$y, usa$x, usa$y)
dens = subset(dens, v == 1)

plt1 = ggplot(dens,  aes(x, y, fill = z)) +
       geom_raster(show.legend = sl) +
       scale_fill_gradientn(colours = terrain.colors(100))+
       xlab('Longitude') + ylab('Latitude') +
       ggtitle('Precipitation anomaly, smoothed and clipped') +
       geom_contour(aes(z = z), color = ccol, show.legend = T) +
       JOPS_theme() +
       coord_fixed() +
       theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             plot.title = element_text(size = capsize),
             axis.title.x = element_text(size = capsize),
             axis.title.y = element_text(size = capsize))

grid.arrange(plt1, plt2, nrow = 2, ncol = 1)
