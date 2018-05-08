#============================
# Pressure Field
# Alwin Wang
#----------------------------

#--- Determine smaller vector to plot over ----
# temp <- long_bndry %>% filter(!is.na(wnum)); temp = temp$s
# median(diff(temp)) = 0.005745228
# But high aspect ratio, cannot be used

pres_h <- long_seshdata %>% 
  # Remove non session node points, after original elements only
  group_by(enum) %>%
  # Count the number of nodes per element on the surface
  add_count(wall) %>%
  filter(n == 2, wall) %>%
  # Find the average 'height' of the rectangle with base on surface
  mutate(base = EucDist(x, y),
         aveh = area/base,
         ar = base/aveh) %>%
  # Filter results to one per element
  filter(!is.na(base))
median(pres_h$aveh)

# Inital nodes and elements
pres_mesh <- filter(long_seshdata, wall)
pres_mesh <- list(
  nodes = unique(pres_mesh$nnum),
  elements = unique(pres_mesh$enum))
# Plot
pres_mesh$mesh <- filter(long_meshdata, enum %in% pres_mesh$elements)
ggplot() + 
  geom_point(aes(x, y), pres_mesh$mesh, alpha = 0.2) + 
  geom_polygon(aes(x, y, group = enum), 
               filter(pres_mesh$mesh, node) %>% arrange(ncorner), fill = NA, colour = "black") + 
  coord_fixed()

# Next ndoes and elements
pres_mesh$nodes =unique(
  long_seshdata[long_seshdata$enum %in% pres_mesh$elements,]$nnum)
pres_mesh$elements =unique(
  long_seshdata[long_seshdata$nnum %in% pres_mesh$nodes,]$enum)
# Plot
pres_mesh$mesh <- filter(long_meshdata, enum %in% pres_mesh$elements)
ggplot() + 
  geom_point(aes(x, y), pres_mesh$mesh, alpha = 0.2) + 
  geom_polygon(aes(x, y, group = enum), 
               filter(pres_mesh$mesh, node) %>% arrange(ncorner), fill = NA, colour = "black") + 
  coord_fixed()
