#============================
# Pre-process File Data
# Alwin Wang
#----------------------------

#---- Elements ----
# Determine node points for the elements
# Check that all elements are quadrangles
if (as.logical(sum(elements$shapetag != "<Q>"))) warning("Not all elements are Quadrangles") 
# Manipulate data
meshlongdata <- elements %>%
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
temp <- left_join(mesh, meshlongdata, by = c("x", "y"))
# NOTE: THERE WILL BE MANY DUPLICATES e.g. five elements meet at one node
nrow(filter(temp, is.na(ncorner))) + nrow(sessionlongdata) - nrow(mesh)

#--- Airfoil Surface ----
