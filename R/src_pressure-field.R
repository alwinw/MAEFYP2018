#============================
# Pressure Field
# Alwin Wang
#----------------------------

#--- Determine smaller vector to plot over ----
# temp <- long_bndry %>% filter(!is.na(wnum)); temp = temp$s
# median(diff(temp)) = 0.005745228
# But high aspect ratio, cannot be used

avecellh <- long_meshdata %>% 
  # Remove non session node points, after original elements only
  filter(!is.na(nnum)) %>%
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

  
