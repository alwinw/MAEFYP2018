#============================
# Plot Variables of Interest
# Alwin Wang
#----------------------------

# Leading Edge
ggplot(long$threaddata, aes(x, y, colour = local)) +
  geom_polygon(aes(x, y), fill = NA, colour = "black", alpha = 0.5,
               data = airfoildata$bndry) +
  geom_polygon(aes(x, y, group = enum, colour = local), fill = NA,
               data = long$threaddata %>% filter(!is.na(nnum)) %>% arrange(enum, ncorner)) +
  geom_point(alpha = 0.2) +
  geom_text(aes(elabx, elaby, label = enum, size = area), alpha = 0.5) +
  geom_label(aes(x, y, label = nnum), size = 3, alpha = 1) +
  coord_fixed(
    xlim = c(-0.45, -0.25),
    ylim = c(-0.055, 0.085),
    expand = FALSE) +
  scale_colour_gradientn(colours = spectralpalette(20)) +
  scale_size(guide = "none", range = c(1*3, 6*10)) +
  ggtitle("Mesh: Leading Edge")

# Trailing Edge
ggplot(long$threaddata, aes(x, y, colour = local)) +
  geom_polygon(aes(x, y), fill = NA, colour = "black", alpha = 0.5,
               data = airfoildata$bndry) +
  geom_polygon(aes(x, y, group = enum, colour = local), fill = NA,
               data = long$threaddata %>% filter(!is.na(nnum)) %>% arrange(enum, ncorner)) +
  geom_point(alpha = 0.2) +
  geom_text(aes(elabx, elaby, label = enum, size = area), alpha = 0.5) +
  geom_label(aes(x, y, label = nnum), size = 3, alpha = 1) +
  coord_fixed(
    xlim = c(0.5, 0.675),
    ylim = c(-0.12, 0.024),
    expand = FALSE) +
  scale_colour_gradientn(colours = spectralpalette(20)) +
  scale_size(guide = "none", range = c(1*3, 6*10)) +
  ggtitle("Mesh: Trailing Edge")






















myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
theme_set(theme_bw())

#--- Sample plot of long_seshdata ----
long_meshplot <- ggplot(long_seshdata, aes(x, y, colour = enum)) +
  # Element edges
  geom_polygon(aes(group = enum), fill = NA) +
  # Element labels
  geom_text(aes(x = elabx, y = elaby, label = enum, size = area), alpha = 0.5) +
  # Node labels
  geom_text(aes(x = x, y = y, label = nnum, size = area*0.2), 
            data = long_seshdata %>% ungroup() %>% group_by(nnum) %>% top_n(-1, area)) +
  scale_size(guide="none") +
  # Labels
  xlab("x") + ylab("y")
# Plots
long_meshplot + # Full mesh
  coord_fixed(
    xlim = c(min(long_seshdata$x), max(long_seshdata$x)),
    ylim = c(min(long_seshdata$y), max(long_seshdata$y)))
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
  geom_path(aes(x, y, colour = enum, group = enum), data = long_seshdata) +
  # Element labels
  geom_text(aes(x = elabx, y = elaby, label = enum, size = area), alpha = 0.5, data = long_seshdata) +
  # Node labels
  geom_text(aes(x = x, y = y, label = nnum, size = area*0.2), 
            data = long_seshdata %>% ungroup() %>% group_by(nnum) %>% top_n(-1, area)) +
  scale_size(guide="none") + 
  geom_text(aes(x = x, y = y, label = jnum), data = mesh) +
  # Labels
  xlab("x") + ylab("y")

ggplot(mesh, aes(x, y, group = enum, colour = jnum)) +
  geom_path() +
  geom_point() + # LE
  coord_fixed(
    xlim = c(-0.4, -0.3),
    ylim = c(-0.04, 0.08))

#--- Boundary ----
ggplot(long_bndry, aes(x = snum, y = c(0, diff(s)), group = is.na(wsnum), colour = s)) +
  geom_point() +
  geom_path() +
  scale_colour_gradientn(colours = myPalette(100)) +
  ylim(0.00001, NA)

#--- Dump Files ----
dumpfileplot <- filter(dumpfile$flowfield)
dumpfileplot <- cbind(dumpfileplot, mesh)
# Plot
ggplot(dumpfileplot, aes(x, y)) +
  geom_point(aes(colour = p)) +
  coord_fixed() +
  scale_colour_gradientn(colours = myPalette(100))

ggplot(dumpfileplot, aes(x, y)) +
  geom_point(aes(colour = t, alpha = abs(t))) +
  coord_fixed(
    xlim = c(-0.55, 0.65),
    ylim = c(-0.3, 0.3)) +
  scale_colour_gradientn(colours = myPalette(100))

