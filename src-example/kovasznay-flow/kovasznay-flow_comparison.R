#============================#
# Kovasznay Flow
# Alwin Wang
#----------------------------#

#--- Set Up ----
# Use rstudioapi to get saved location of this file; use str to print structures
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
# Libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
# Custom ggplot2 setup
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

dns <- read.table("kovas1.flddump", skip = 10)
mesh <- read.table("kovas1.msh", skip = 1)

numeric <- cbind(mesh, dns)
colnames(numeric) <- c("x", "y", "u", "v", "p", "k", "l", "m", "n", "o")

numericplot <- numeric[!duplicated(numeric[,1:2]),]

ggplot(numericplot, aes(x, y, z = v)) +
  geom_point(aes(colour = o)) +
  # stat_contour(aes(x, y, z = v)) +
  scale_colour_gradientn(colours = spectralpalette(20))

Re = 40
analytic <- numeric[,1:2] 
analytic %<>%
  mutate(lambda = Re/2 - sqrt(Re^2/4 + 4*pi^2)) %>%
  mutate(u = 1 - exp(lambda*x)*cos(2*pi*y),
         v = lambda/(2*pi)*exp(lambda*x)*sin(2*pi*y),
         p = 1/2*(1-exp(2*lambda*x))) %>%
  mutate(k = (lambda^3/(2*pi)-2*pi*lambda)*exp(lambda*x)*sin(2*pi*y) * -1,
         l = (lambda^2-4*pi^2)*exp(lambda*x)*cos(2*pi*y) * -1,
         m = -lambda/2*exp(2*lambda*x),
         n = 0,
         o = (lambda^2/(2*pi)-2*pi)*exp(lambda*x)*sin(2*pi*y) * -1) %>% # Why do I need to x-1 ??
  select(-lambda)

ggplot(analytic, aes(x , y)) +
  geom_point(aes(colour = o)) +
  scale_colour_gradientn(colours = spectralpalette(20))

summary(abs(numeric - analytic))
