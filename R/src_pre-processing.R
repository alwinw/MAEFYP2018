#============================
# Pre-process File Data
# Alwin Wang
#----------------------------

# In wall mesh, the 5th node of the 1st element = 1st node of the 2nd element
# To join over (x, y), these duplicate coordinates need to be removed!
unixy_wallmsh <- long_wall[!duplicated(select(long_wall, x, y)),]


#--- Elements ----
# Determine node points for the elements
# Check that all elements are quadrangles
if (as.logical(sum(elements$shapetag != "<Q>"))) warning("Not all elements are Quadrangles") 
# Manipulate data
long_seshdata <- elements %>%
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
  arrange(enum, ncorner) %>% 
  # Calculate area
  mutate(area = x*lead(y) - lead(x)*y) %>%
  mutate(area = 1/2*abs(sum(ifelse(is.na(area), 0, area)))) %>%
  # Join with wall mesh
  left_join(., unixy_wallmsh, by = c("x", "y")) %>%
  mutate(wall = !is.na(wnum)) %>%
  filter(ncorner != "n5")

#--- Mesh ----
# Combine the mesh (N order poly) with original elements
# Determine which mesh data belong to which
mesh$mnum = 1:nrow(mesh)
# Determine which mesh nodes are wall mesh nodes
# Join with unique (x, y) for temp_wallmsh
long_meshdata <- left_join(mesh, unixy_wallmsh, by = c("x", "y"))
long_meshdata$wall = !is.na(long_meshdata$wnum)
# Check the number of rows has not changed
if (nrow(long_meshdata) != nrow(mesh)) {
  warning("Number of nodes has changed!")}
# Check that all wall mesh nodes were found
if (nrow(unixy_wallmsh) != 
    nrow(long_meshdata %>% filter(wall) %>% select(wnum) %>% unique())) {
  # Note: need unique because there are duplicate (x, y) points in the mesh file
  #       since it is a print out for each element and N x N
  warning("Not all wall mesh nodes found")}

# Determine which mesh nodes are session node numbers
long_meshdata <- left_join(long_meshdata, long_seshdata, by = c("enum", "x", "y"))
long_meshdata$node = !is.na(long_meshdata$nnum)
if (nrow(filter(long_seshdata)) !=
    nrow(long_meshdata %>% select(enum, ncorner) %>% filter(!is.na(ncorner)) %>% unique())) {
  warning("Not all mesh nodes found")}
# Double check mesh node order correct (ONLY IF N = 5 FOR SPACE SPACING)
# long_meshdata %<>% cbind(., data.frame(
#     check = rep(c("n1", rep(NA,3), "n2", rep(NA, 15), "n4", rep(NA, 3), "n3"), nrow(long_meshdata)/25))) %>%
#   mutate(check = as.character(check)) %>%
#   mutate(error = is.na(check) & is.na(ncorner)) %>%
#   mutate(error = ifelse(is.na(check), FALSE, ncorner != check))
# if (sum(long_meshdata$error) != 0) {
#   warning("Node spacing not as expected")}

#--- History File ----
# long_his <- left_join(his, long_bndry, by = "bnum")

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
ggplot(long_his, aes(s, p, group = t, colour = t)) + geom_line() +
  scale_colour_gradientn(colours = myPalette(100))
ggplot(long_his, aes(s, t, group = t, colour = p)) + geom_line() +
  scale_colour_gradientn(colours = myPalette(100))

#--- Step Variables ----
# acceleration value

# Pressue
