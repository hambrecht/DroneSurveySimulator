# Before running this script to make sure that the environment is set up correctly.

# Check if renv is installed
if (!requireNamespace("renv", quietly = TRUE)) {
  message("The 'renv' package is not installed. Installing it now...")
  install.packages("renv")
} else {
  message("The 'renv' package is already installed.")
}
library(renv)

# Restore the environment and install required packages
renv::restore(prompt = interactive())