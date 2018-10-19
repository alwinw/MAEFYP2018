#============================#
# Kovasznay Flow
# Alwin Wang
#----------------------------#

#--- Set Up ----
# Use rstudioapi to get saved location of this file; use str to print structures
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Requres devtools, rstudioapi
# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
# Custom ggplot2 setup
theme_set(theme_bw())                                           # Set black and white theme
spectralpalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

resultslist <- unlist(strsplit(
  list.files("results", pattern = "flddump"), "*.flddump"))

convergence <- data.frame()
Re = 40
for (result in resultslist) {
  # Numerical result
  dns  <- read.table(paste0("results/", result, ".flddump"), skip = 10)
  mesh <- read.table(paste0("results/", result, ".mshi"),    skip = 1)
  numeric <- cbind(mesh[,1:2], dns)
  colnames(numeric) <- c("x", "y", "u", "v", "p", "k", "l", "m", "n", "o")
  # Analytical result
  analytic <- numeric[,1:2] 
  analytic %<>%
    mutate(lambda = Re/2 - sqrt(Re^2/4 + 4*pi^2)) %>%
    mutate(u = 1 - exp(lambda*x)*cos(2*pi*y),
           v = lambda/(2*pi)*exp(lambda*x)*sin(2*pi*y),
           p = 1/2*(1-exp(2*lambda*x))) %>%
    mutate(k = (lambda^3/(2*pi)-2*pi*lambda)*exp(lambda*x)*sin(2*pi*y),
           l = (lambda^2-4*pi^2)*exp(lambda*x)*cos(2*pi*y),
           m = -lambda/2*exp(2*lambda*x),
           n = 0,
           o = (lambda^2/(2*pi)-2*pi)*exp(lambda*x)*sin(2*pi*y)) %>% # Why do I need to x-1 ??
    select(-lambda)
  # Comparison
  summary <- as.data.frame(do.call(cbind, lapply(abs(analytic - numeric), summary)))
  summary$num  <- rownames(summary)
  summary$name <- result 
  summary$N_P  <- as.numeric(gsub('[a-z]+', '', result))
  rownames(summary) <- NULL
  print(summary)
  # Convergence
  convergence <- rbind(convergence, summary)
}

# Plot results
norminf <- filter(convergence, num == "Max.") %>%
  select(-x, -y, -num, -name) %>%
  gather(var, norm_inf, -N_P) %>%
  mutate(pressure = var %in% c("p", "m", "n")) %>% 
  mutate(pressure = ifelse(pressure, "Pressure", "Velocity")) %>% 
  mutate(var = ifelse(var == "k", "dodx", var),
         var = ifelse(var == "l", "dody", var),
         var = ifelse(var == "m", "dpdx", var),
         var = ifelse(var == "n", "dpdy", var) )

norminf$var = factor(norminf$var, levels = c("p", "dpdx", "dpdy", "u", "v", "o", "dodx", "dody"))

ggplot(norminf, aes(N_P, norm_inf, group = var, colour = var)) +
  geom_line(aes(linetype = var)) +
  geom_point(aes(shape = var)) +
  scale_y_continuous(trans='log10') +
  # scale_colour_discrete(name = "Var", )
  facet_wrap(~pressure, scales = "free_y") +
  ggtitle("Kovasznay Flow Convergence")

ggsave(paste0("kovas.png"),
       scale = 2.5, width = 10, height = 6, units = "cm", dpi = 300)
