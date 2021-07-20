# Sample log-periodogram spectra (Phoneme data)
# A graph in the book 'Practical Smoothing. The Joys of P-splines'
# Paul Eilers and Brian Marx, 2019

library(ggplot2)
library(gridExtra)
library(JOPS)
library(ElemStatLearn)

# Get the data
data(phoneme)
phon = phoneme[, 257]
aa = which(phon == "aa")
ao = which(phon == "ao")
n = 150
X = t(as.matrix(phoneme[, 1:n]))
freq = (1:n) * 1000 / 32

# Selecting example spectra with 'aa' and 'ao'
sel = c(aa[1], ao[1])
X = X[, sel]

# Prepare dataframe for ggplot
id = as.factor(rep(c("aa", "ao"), each = n))
F1 = data.frame(x = freq, y = as.vector(X), id = id)

plt1 = ggplot(data = F1)  +
  geom_line(aes(x = x, y = y, col = id), size = .4) +
  ggtitle("Examples of phoneme spectra") + ylab("") + xlab("Frequency (Hz)") +
  JOPS_theme() +
  labs(color = 'Vowel') +
  scale_color_manual(values = c('blue', 'red'))

# Make and save pdf
grid.arrange(plt1, ncol = 1, nrow = 1)
