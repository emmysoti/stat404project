# renv::install("tidyverse")
library(tidyverse)
library(rio)

# Set master seed from data_generation/master_seed.R
source("src/data_generation/master_seed.R")
set.seed(master_seed)

# Generate random seeds for the forty datasets:
seeds <- sample(1:1000, size = 40, replace = FALSE)

# Separate the seeds into four vectors for later implementation
seeds1 <- seeds[1:10] # linear, low
seeds2 <- seeds[11:20] # linear, high
seeds3 <- seeds[21:30] # quadratic, low
seeds4 <- seeds[31:40] # quadratic, high

# Set n, mean, and SD of generated data
n <- 1000
mean_x <- 0
sd_x <- 3

# Set mean and SD of noise w/ level
low_mean <- 0
low_sd <- 0.5
high_mean <- 0
high_sd <- 5

# Define data generating function
generate.data <-
  function(seeds, relationship_type, noise_type) {
    if (relationship_type != "linear" && relationship_type != "quadratic") {
      stop("Invalid data relationship type")
    }
    if (noise_type != "low" && noise_type != "high") {
      stop("Invalid noise relationship type")
    }

    data.list <- list()

    for (i in 1:10) {
      # Set seeds from randomly generated seeds
      set.seed(seeds[i])

      # Generate data
      x <- rnorm(n, mean = mean_x, sd = sd_x)
      noise <- rnorm(
        n,
        if_else(noise_type == "low", low_mean, high_mean),
        if_else(noise_type == "low", low_sd, high_sd)
      )
      if (relationship_type == "linear") {
        y <- 2 * x + 3 + noise
      } else {
        y <- 2 * x^2 + 3 * x + 4 + noise
      }

      # Assign x and y values to dataframe, save in appropriate list
      df <- data.frame(x, y)
      data.list[[i]] <- df
    }

    return(data.list)
  }

##### Generate 10 datasets: linear & low noise
linear_low <- generate.data(seeds1, "linear", "low")


##### Generate 10 datasets: linear & high noise
linear_high <- generate.data(seeds2, "linear", "high")


##### Generate 10 datasets: non-linear & low noise
quadratic_low <- generate.data(seeds3, "quadratic", "low")


##### Generate 10 datasets: non-linear & high noise
quadratic_high <- generate.data(seeds4, "quadratic", "high")
