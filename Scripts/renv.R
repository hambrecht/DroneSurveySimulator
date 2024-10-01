# Load renv package
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
library(renv)

# Initialize renv in the project
renv::init(bare=TRUE)

renv::install(c("here", "sf", "tidyr", "units", "purrr", "dplyr", "dsims", "Distance", "dsm", "knitr", "jsonlite", "rlang"))


# Snapshot the current state of the library
renv::snapshot()

# Print a message indicating the setup is complete
cat("renv setup complete. The project dependencies have been installed and snapshotted.\n")