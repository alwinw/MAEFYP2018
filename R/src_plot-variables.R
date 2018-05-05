#============================
# Plot Variables of Interest
# Alwin Wang
#----------------------------

#--- Sample plot of long_meshdata ----
long_meshplot <- ggplot(long_meshdata, aes(x, y, colour = enum)) +
  # Element edges
  geom_path(aes(group = enum)) +
  # Element labels
  geom_text(aes(x = elabx, y = elaby, label = enum, size = area), alpha = 0.5) +
  # Node labels
  geom_text(aes(x = x, y = y, label = nnum, size = area*0.2), 
            data = long_meshdata %>% ungroup() %>% group_by(nnum) %>% top_n(-1, area)) +
  scale_size(guide="none") +
  # Labels
  xlab("x") + ylab("y")
# Plots
long_meshplot + # Full mesh
  coord_fixed(
    xlim = c(min(long_meshdata$x), max(long_meshdata$x)),
    ylim = c(min(long_meshdata$y), max(long_meshdata$y)))
long_meshplot + # Airfoil only
  coord_fixed(
    xlim = c(-0.55, 0.65),
    ylim = c(-0.3, 0.3))
long_meshplot + # LE
  coord_fixed(
    xlim = c(-0.4, -0.3),
    ylim = c(-0.04, 0.08))
long_meshplot + # TE
  coord_fixed(
    xlim = c(0.5, 0.65),
    ylim = c(-0.1, 0.02))

#--- Determine mesh layout ----
ggplot() +
  # Element edges
  geom_path(aes(x, y, colour = enum, group = enum), data = long_meshdata) +
  # Element labels
  geom_text(aes(x = elabx, y = elaby, label = enum, size = area), alpha = 0.5, data = long_meshdata) +
  # Node labels
  geom_text(aes(x = x, y = y, label = nnum, size = area*0.2), 
            data = long_meshdata %>% ungroup() %>% group_by(nnum) %>% top_n(-1, area)) +
  scale_size(guide="none") + 
  geom_text(aes(x = x, y = y, label = j), data = mesh) +
  # Labels
  xlab("x") + ylab("y")
