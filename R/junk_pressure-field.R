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
  geom_path(aes(xdash, ydash, group = snum), data = long_walloffset, colour = "red") +
  coord_fixed()

# Next nodes and elements
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
  geom_path(aes(xdash, ydash, group = snum), data = long_walloffset, colour = "red") +
  coord_fixed()

# ACTUALLY I DO NOT THINK THAT I SHOULD INCLUDE THE WALL POINTS IN THE INTERPOLATION
# THESE MAY CAUSE COLINEARITY PROBLEMS AND I ALREADY KNOW THE VORTICITY AT THAT POINT
# NO NEED TO INTERPOLATE!!

long_walloffset <- rbind(
  long_walloffset %>% select(-xdash, -ydash),
  long_walloffset %>% select(-x, -y) %>% rename(x = xdash, y = ydash))

long_dump <- cbind(mesh, dumpfile$flowfield)
# Remove duplicated rows
long_dump <- long_dump[!duplicated(select(long_dump, -enum, -jnum)),]

ptm <- proc.time()
# Use whole field, ~2.66 seconds
long_interp <- as.data.frame(
  interpp(x = long_dump$x, y = long_dump$y, z = long_dump$p,
  xo = long_walloffset$x, yo = long_walloffset$y,
  linear = FALSE,
  duplicate = "strip"))
proc.time() - ptm


long_pres <- filter(long_dump, enum %in% pres_mesh$elements)

ptm <- proc.time()
# Use reduced field, ~0.112 seconds
long_interp_pres <- as.data.frame(
  interpp(x = long_pres$x, y = long_pres$y, z = long_pres$p,
          xo = long_walloffset$x, yo = long_walloffset$y,
          linear = FALSE,
          duplicate = "strip"))
proc.time() - ptm

# SOME ERRORS ARE 20+% !!
hist((long_interp$z - long_interp_pres$z)/long_interp$z * 100, breaks = 50)

ggplot(long_interp, aes(x = x, y = y, colour = -z)) +
  geom_point() +
  coord_fixed()

test_linear_interp <- as.data.frame(
  interpp(x = long_dump$x, y = long_dump$y, z = long_dump$p,
  xo = long_walloffset$x, yo = long_walloffset$y,
  linear = TRUE,
  duplicate = "strip"))

# Interpolate back onto original (x, y)
check_interp <- as.data.frame(
  interpp(x = long_interp$x, y = long_interp$y, long_interp$z,
          xo = long_pres$x, yo = long_pres$y,
          linear = FALSE,
          duplicate = "strip"))
check_interp = (check_interp$z - long_pres$p)/long_pres$p
check_interp = check_interp[!is.na(check_interp)]
hist(check_interp, breaks = 50)


# Interpolate back onto original (x, y)
check_interp_pres <- as.data.frame(
  interpp(x = long_interp_pres$x, y = long_interp_pres$y, long_interp_pres$z,
          xo = long_pres$x, yo = long_pres$y,
          linear = FALSE,
          duplicate = "strip"))
check_interp_pres = (check_interp_pres$z - long_pres$p)/long_pres$p
check_interp_pres = check_interp_pres[!is.na(check_interp_pres)]
hist(check_interp_pres, breaks = 50)
mean(check_interp_pres)
