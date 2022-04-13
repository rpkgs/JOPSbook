# PS-ANOVA (Ethanol data)
# A graph from the book "Practical Smoothing. The Joys of P-splines"
# Paul Eilers and Brian Marx, 2019

library(data.table)
library(SpATS)
library(ggplot2)
library(JOPS)
library(gridExtra)
library(SemiPar)
library(rlist)
library(reshape2)

# Get the data
data(ethanol)

# Fit the PS-ANOVA model
ctrl =  list(tolerance = 1e-05, monitoring = F, maxit = 100, update.psi = F)
fit <- SpATS.nogeno(response = "NOx",
						spatial = ~PSANOVA(E, C,  nseg = c(20, 20), nest.div = c(2, 2)),
            data = ethanol,	control = ctrl)

# Report effective dimensions
print(summary(fit))

# Extract the surfaces
Tr = obtain.spatialtrend(fit)
G = Tr$pfit

plot_surf = function(Mu, Tr, name = '') {
  # Helper function to plot and annotate a surface
  rownames(Mu) = Tr$row.p
  colnames(Mu) = Tr$col.p
  dens <- melt(Mu)
  names(dens) = c('x', 'y', 'z')
  sl = T
  ccol = 'blue'
  plt = ggplot(dens,  aes(x, y, fill = z)) +
      geom_raster(show.legend = sl) +
      scale_fill_gradientn(colours = terrain.colors(100))+
      xlab('Compression') + ylab('Equivalence') +
      ggtitle(name) +
      geom_contour(aes(z = z), color = ccol, show.legend = T) +
      JOPS_theme() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  return(plt)
}

# Make panels for the component surfaces
plts = list()
Mu = G$bc
plt = plot_surf(Mu, Tr, 'Bilinear')
plts = list.append(plts, plt)
Mu = G$uhv + G$fv
plt = plot_surf(Mu, Tr, 'Additive plus vaying-coefficient')
plts = list.append(plts, plt)
Mu = G$vhu + G$fu
plt = plot_surf(Mu, Tr, 'Additive plus vaying-coefficient')
plts = list.append(plts, plt)
Mu = G$fuv
plt = plot_surf(Mu, Tr, 'Interaction')
plts = list.append(plts, plt)
Mu = Tr$fit
plt = plot_surf(Mu, Tr, 'Fitted surface')
plts = list.append(plts, plt)

# Panel for effective dimensions
ed = fit$eff.dim
enm = names(ed)
De = data.frame(name = factor(enm, ordered = T, levels = enm), ed = ed)

plt = ggplot(data = De, aes(name, ed)) +
      geom_col(fill ='wheat3') +
      coord_flip() +
      ggtitle('Effective dimensions') +
      xlab('') +  ylab('') +
      JOPS_theme()
plts = list.append(plts, plt)

# Print the panels
pg = grid.arrange(grobs = plts, nrow = 3, ncol = 2)
plot(pg)



