# _targets.R

#setup
library(targets)
source("R/functions.R")

# target-specific options
tar_option_set(packages = c(
  "dplyr"
))

# list of target objects
list(
  tar_target(data, data.frame(x = sample.int(100), y = sample.int(100))),
  tar_target(summary, summ(data)) # Call your custom functions as needed.
)
