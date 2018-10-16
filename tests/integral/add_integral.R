#============================#
# Integration Comparision
# Alwin Wang
#----------------------------#

#--- * Integral Output                                            ----
# Load semtex integration output
dump$inte <- LoadIntegral(data_dump$folder, data_dump$dumpfile,
                          vars = c("u", "v", "p", letters[11:15]))
dump$inte$elem <- dump$inte$elem %>% 
  select(enum, u, v, p, k, l, m, n, o) %>% 
  rename(dodx = k, dody = l, dpdx = m, dpdy = n)
# Create comparison using mass matrix
comp <- list()
comp$comp <- dump$dump %>% 
  select(u, v, p, dodx, dody, dpdx, dpdy, o)
comp$comp      = comp$comp*dump$dump$mass
comp$comp$enum = dump$dump$enum
# Summarise
comp$elem <- comp$comp %>% 
  group_by(enum) %>% 
  summarise_all(.funs = "sum")
comp$total <- comp$elem %>% 
  select(-enum) %>% 
  summarise_all(.funs = "sum")
# Compare output
max(abs(dump$inte$elem - comp$elem ))
dump$inte$total$intvar - comp$total