#============================
# Pre-process File Data
# Alwin Wang
#----------------------------

#---- Elements ----
# Determine node points for the elements
# Check that all elements are quadrangles
if (as.logical(sum(elements$shapetag != "<Q>"))) warning("Not all elements are Quadrangles") 
# Manipulate data
long_meshdata <- elements %>%
  # Remove shape tag since checked shape already
  select(-shapetag) %>%
  # Gather nodes 1-5 into a single column
  gather(ncorner, nnum, -enum) %>%
  # Combine with nodes to get (x,y) for each element corner
  left_join(., nodes, "nnum") %>%
  # Order data
  group_by(enum) %>%
  # Add columns for labels
  mutate(
    elabx = mean(x),
    elaby = mean(y)) %>%
  # Close the loop (without effecting mean)
  rbind(., mutate(filter(., ncorner=="n1"), ncorner = "n5")) %>%
  arrange(enum, ncorner)
  # Calculate area later

#--- Mesh ----
# Combine the mesh (N order poly) with original elements
mesh$mnum = 1:nrow(mesh)
temp <- left_join(mesh, long_meshdata, by = c("x", "y"))
# NOTE: THERE WILL BE MANY DUPLICATES e.g. five elements meet at one node
nrow(filter(temp, is.na(ncorner))) + nrow(long_meshdata) - nrow(mesh)

#--- Airfoil Surface ----
# Clean up wall mesh
wallmsh$wnum = 1:nrow(wallmsh)
# Determine the chord line
chordlm = lm(y~x, chord)    # class lm
# Determine upper and lower surfaces
wallmsh$up = wallmsh$y >= predict(chordlm, wallmsh)
# Order wall file
wallmsh <- arrange(wallmsh, up, x * (up * 2 -1))
wallmsh$wsnum = 1:nrow(wallmsh)
# Make the same number of columns for bndry
bndry$wnum <- NA
bndry$wsnum <- NA
bndry$up = bndry$y >= predict(chordlm, bndry)
# Combine wall and bndry files
long_bndry <- rbind(wallmsh, bndry) %>%
  arrange(-up, -x * (up * 2 -1))
long_bndry$snum = 1:nrow(long_bndry)
ggplot(long_bndry, aes(x, y, color = up)) + geom_path() + geom_point(aes(shape=up)) +
  coord_cartesian(xlim = c(-0.4, -0.3), ylim = c(0.0, 0.08))

#--- Spline Length s ----


# Eventually compare spline length to XFOIL output
