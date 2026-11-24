
##############################################################################################
## Please only run this script once! It generates the master seed for the rest of the seeds!
##############################################################################################

# This package connects to a true random number generator at random.org
renv::install("random")
library(random)

# master seed is 297
master_seed <- randomNumbers(n = 1, min = 1, max = 1000, col = 1)[1,1]

