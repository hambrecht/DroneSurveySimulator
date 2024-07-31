# DroneSurveySimulator
This project analyses aerial survey data of moose populations using distance sampling methods to estimate population density. The estimated density is then used in a simulated aerial survey. The project is structured to separate different stages of the workflow, ensuring clarity and modularity.

## Project Structure

 README.md # Project overview and instructions
├── Ruserdata/ # Directory for raw and processed data
├── Scripts/ # R scripts for different stages of the analysis
│ ├── data_preparation.R # Scripts for data cleaning and preparation
│ ├── density_analysis.R # Scripts for statistical analysis and modelling
│ ├── simulation.R # Script for running simulations
│ ├── visualisation.R # Script for visualisations and reporting
└── Output/ # intermediate files, plots, and reports

## Description of Files and Directories

- **Ruserdata/**: Contains raw data (`raw/`) and processed data (`processed/`).
- **scripts/**: Organised by functionality:
  - **data_processing/**: Includes scripts for cleaning and preparing the raw data (`data_cleaning.R`).
  - **analysis/**: Contains scripts for deriving population density using distance sampling (`density_analysis.R`).
  - **simulation/**: Includes scripts to simulate aerial surveys using the derived population density (`simulation.R`).
  - **visualisation/**: Scripts for generating plots, tables, and reports (`visualisation.R`).
- **Output/**: Stores the outputs from analyses, simulations, and visualisations.

## Getting Started

1. **Install Required Packages**: Ensure you have the required R packages installed by running:
   ```R
   install.packages(c("Distance", "dplyr", "dsims", "dsm", "here", "knitr", "pbapply", "purrr", "sf", , "tidyr", "units")) # List all required packages
   ```
2. **Running the Scripts**: Start with `data_cleaning.R` to prepare the data.
- Proceed with `density_preparation.R` for density estimation.
- Run `simulation.R` to perform the simulated survey.
- Use `visualisation.R` to generate visual outputs.
