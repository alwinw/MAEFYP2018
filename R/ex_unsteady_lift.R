#============================#
# Plot Boundary Conditions
# Alwin Wang
#----------------------------#

# source("plot_BCs.R")
source("../../R/plot_BCs.R")

# t = seq(0, D_T*N_STEP, length.out = N_STEP + 1)
t = seq(0, D_T*N_STEP, length.out = 10)
flow <- data.frame(
  t = t,
  Vinf = BC_inflow_u(t)
)
theta = seq(0, pi, length.out = 10)
alpha = 4*pi/180



# Needs to be a discretised problem ugh. Back to Gamma and x


flx <- read.table(file = "../analysis/whole_10000/results/NACA0012r.flx")
# flx <- read.table(file = "results/NACA0012r.flx")
colnames(flx) <- c("step", "t", 
                   paste0(c("Fpre", "Fvis", "Ftot"), "x"),
                   paste0(c("Fpre", "Fvis", "Ftot"), "y"),
                   paste0(c("Fpre", "Fvis", "Ftot"), "z"))
flx <- flx %>% 
  mutate(
    Vinf   = BC_inflow_u(t),
    Gamma  = pi*alpha*Vinf,
    Gammax = Gamma/4) %>% 
  mutate(
    qs = Vinf*Gamma,
    am = (lead(Gammax) - lag(Gammax))/(lead(t) - lag(t)))

ggplot(flx, aes(t)) +
  annotate("rect", xmin=0.05, xmax=0.25, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  annotate("rect", xmin=0.85, xmax=1.05, ymin=-Inf, ymax=Inf,
           fill = "grey90", alpha = 0.5) +
  geom_line(aes(y = Ftoty, colour = "DNS Lift"), size = 0.8) +
  geom_line(aes(y = qs, colour = "Quasi-Steady"), linetype = "dashed") +
  geom_line(aes(y = am, colour = "Apparent Mass")) +
  # geom_hline(aes(colour = "XFOIL", yintercept = 0.1315)) +
  # geom_hline(aes(colour = "Thin Airfoil", yintercept = 2*pi*alpha))
  ylab("Lift Force") +
  xlab("Time (s)") +
  scale_x_continuous(breaks = c(0, seq(0.05, 1.5, 0.2), 1.5),
                     minor_breaks = c(0, seq(0.05, 1.5, 0.1), 1.5),
                     labels = function(x) ifelse(x==0 | x==1.5, "", sprintf("%.2f", x)),
                     limits = c(0, 1.5))
ggsave("Unsteady.png", scale = 1,
       width = 15, height = 7, units = "cm", dpi = 300)

ggsave("Unsteady.png", scale = 1.5,
       width = 15, height = 6, units = "cm", dpi = 600)

# WTF is viscous lift negative

ggplot(flx, aes(t)) +
  geom_line(aes(y = Fvisy, colour = "DNS Viscous Lift")) +
  ylab("Lift Force") +
  xlab("Time (s)") +
  scale_x_continuous(breaks = c(0, seq(0.05, 2.05, 0.2)),
                     minor_breaks = c(0, seq(0.05, 2.05, 0.1)),
                     labels = function(x) ifelse(x==0, "", sprintf("%.2f", x)),
                     limits = c(0, 2))
ggsave("Viscous Lift.png", scale = 1.5,
       width = 15, height = 7, units = "cm", dpi = 600)
             