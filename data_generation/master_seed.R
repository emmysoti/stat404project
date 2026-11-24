
# master will be a new number every time this is run. master_seed is simply one of the earlier values it took

# This package connects to a true random number generator at random.org
renv::install("random")
library(random)

# master is 979 
master <- randomNumbers(n = 1, min = 1, max = 1000, col = 1)[1,1]
master_seed <- 979

