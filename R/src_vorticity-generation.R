#============================#
# Vorticity Generation
# Alwin Wang
#----------------------------#

#--- Convert session into a session long mesh ----
LongSession <-function(session) {
  # Elements
  # Check that all elements are quadrangles
  if (as.logical(sum(session$elements$shapetag != "<Q>"))) warning("Not all elements are Quadrangles") 
  # Manipulate data
  long_seshdata <- session$elements %>%
    # Remove shape tag since checked shape already
    select(-shapetag) %>%
    # Gather nodes 1-5 into a single column
    gather(ncorner, nnum, -enum) %>%
    # Combine with nodes to get (x,y) for each element corner
    left_join(., session$nodes, "nnum") %>%
    # Order data
    group_by(enum) %>%
    # Add columns for labels
    mutate(
      elabx = mean(x),
      elaby = mean(y)) %>%
    # Close the loop (without effecting mean)
    rbind(., mutate(filter(., ncorner=="n1"), ncorner = "n5")) %>%
    arrange(enum, ncorner) %>% 
    # Calculate area
    mutate(area = x*lead(y) - lead(x)*y) %>%
    mutate(area = 1/2*abs(sum(ifelse(is.na(area), 0, area)))) %>%
    # Remove extra ncorner
    filter(ncorner != "n5") %>%
    ungroup()
  # Return the output
  return(long_seshdata)
}

#--- Convert mesh into a long mesh ----
LongMesh <- function(mesh, airfoildata) {
  # Determine which mesh data belong to which
  mesh$mnum = 1:nrow(mesh)
  long_meshdata = mesh
  # Return the output
  return(long_meshdata)
}

#--- Determine local mesh ----
MeshLocal <- function(long_seshdata, long_meshdata) {
  # Determine mesh offset
  long_local <- long_seshdata %>% 
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
  # Determine the minimum average height of the 
  localave <- data.frame(min = min(long_local$aveh), 
                         median = median(long_local$aveh),
                         mean = mean(long_local$aveh))
  # 
  # Inital nodes and elements
  long_localmesh <- filter(long_seshdata, wall)
  long_localmesh <- list(
    nodes = unique(long_localmesh$nnum),
    elements = unique(long_localmesh$enum))
  # Plot
  # long_localmesh$mesh <- filter(long_meshdata, enum %in% long_localmesh$elements)
  # Next nodes and elements
  long_localmesh$nodes =unique(
    long_seshdata[long_seshdata$enum %in% long_localmesh$elements,]$nnum)
  long_localmesh$elements =unique(
    long_seshdata[long_seshdata$nnum %in% long_localmesh$nodes,]$enum)
  # Plot
  long_localmesh$mesh <- filter(long_meshdata, enum %in% long_localmesh$elements)
  ggplot() + 
    geom_point(aes(x, y), long_localmesh$mesh, alpha = 0.2) + 
    geom_polygon(aes(x, y, group = enum), 
                 filter(long_localmesh$mesh, node) %>% arrange(ncorner), fill = NA, colour = "black") + 
    # geom_path(aes(xdash, ydash, group = snum), data = long_walloffset, colour = "red") +
    coord_fixed()
  return(list(localave = localave, mesh = long_localmesh$mesh))
}
