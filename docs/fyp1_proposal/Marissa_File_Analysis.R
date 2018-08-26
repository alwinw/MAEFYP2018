# Airfoils
library(stats)
library(ggplot2)
library(dplyr)
library(tidyr)
# library(pspline)
theme_set(theme_bw())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # requres devtools

heav <- function(t) ifelse(t>0,1,0)

t = seq(0, 1, length.out = 1000)

MOD_ALPHA_X = -((5*pi/2*sin(5*pi*t))*heav(t)+(-5*pi/2*sin(5*pi*t))*heav(t-0.2)+
                  (-5*pi/2*sin(5*pi*t))*heav(t-0.8)+(5*pi/2*sin(5*pi*t))*heav(t-1))

inflow_u = (1/2-1/2*cos(5*pi*t))*heav(t)+(1/2-1/2*cos(5*pi*t+pi))*heav(t-0.2)+(-1/2-1/2*cos(5*pi*t+pi))*heav(t-0.8)+(-1/2-1/2*cos(5*pi*t))*heav(t-1)

plot <- data.frame(t = t, MOD_ALPHA_X = MOD_ALPHA_X, inflow_u = inflow_u)

ggplot(plot) +
  geom_path(aes(x = t, y = MOD_ALPHA_X))

ggplot(plot) +
  geom_path(aes(x = t, y = inflow_u))


mesh = read.table("sp5_nodes.txt")
colnames(mesh) <- c("node", "x", "y", 'z')
bcs = read.table("sp5_bcs.txt")
colnames(bcs) <- c("n", "node", "surf","bc", "group", "bc2")

meshplot <- full_join(mesh, bcs, by="node")

ggplot(mesh, aes(x, y)) + 
  geom_point(alpha=0.1) + 
  # geom_path(alpha=0.1) +
  geom_text(aes(label=node),hjust=0.5, vjust=0.5) +
  coord_fixed()

ggplot(meshplot, aes(x, y)) + 
  geom_point(alpha=0.1) + 
  # geom_path(alpha=0.1) +
  geom_text(aes(label=group),hjust=0.5, vjust=0.5) +
  coord_fixed()

