library(here)
library(knitr)
library(dplyr)
library(dsims)
library(scales)

# Extract key metrics from each simulation summary
extract_metrics <- function(sim) {
    summary_data <- summary(sim, description.summary = FALSE)

    list(
        mean_estimate_density = round(summary_data@individuals$D$mean.Estimate* 1e9, 3), # mean estimate of abundance across simulation
        mean_N = round(summary_data@individuals$summary$mean.n, 0),
        true_density = round(summary_data@individuals$D$Truth* 1e9, 3),
        relative_mean_estimate_density = round(summary_data@individuals$D$mean.Estimate / summary_data@individuals$D$Truth, 3), # mean estimate of abundance across simulation
        percent_bias = round(summary_data@individuals$D$percent.bias, 3), # the percentage of bias in the estimates
        rrmse = round(summary_data@individuals$D$RMSE / summary_data@individuals$D$mean.Estimate, 3), # root mean squared error/no. successful reps
        ci_coverage_prob = round(summary_data@individuals$D$CI.coverage.prob, 3), # proportion of times the 95% confidence interval contained the true value.
        mean_rse = round(summary_data@individuals$D$mean.se / summary_data@individuals$D$mean.Estimate, 3), # the mean standard error of the estimates of abundance
        sd_of_means = round(summary_data@individuals$D$sd.of.means / summary_data@individuals$D$mean.Estimate, 3), # the standard deviation of the estimates
        mean_ER = summary_data@individuals$summary$mean.ER, # mean standard error of the encounter rates cross simulation
        mean_se_ER = summary_data@individuals$summary$mean.se.ER, # standard deviation of the encounter rates across simulation
        N = summary_data@individuals$N$Truth, # number of individuals in the simulation
        cover = NA, # area covered by each design
        sampleblocks = length(summary_data@design.summary$design.type),
        effortKM = round((summary_data@individuals$summary$mean.Effort) / 1000, 2)
    )
}

# List all files containing 'cover' in folder Output/Simulation
files <- list.files(path = here("Output", "Simulation"), pattern = "density_sim-A.*\\.RData$", full.names = TRUE)


# Extract the required numbers from each file name and load the file to extract the summary value
stats <- lapply(files, function(file) {
    # Regular expression to extract the desired components
    matches <- regmatches(file, regexec("QC_([^_]+)_design_([^\\-]+)-density_sim-A(\\d+)", file))

    # Convert to data frame for clarity
    simConfigs <- do.call(rbind, lapply(matches, function(x) data.frame(Overlap = x[2], NSubplots = x[3], Abundance = x[4])))
    # if(is.na(simConfigs[1,1])){simConfigs[1,1] <- "gimbal"}

    # Load the file
    load(file)

    if (!is.na(QC_Sys_sim_density)) {
        summary_data <- extract_metrics(QC_Sys_sim_density)
        Mean_estimated_Density <- summary_data$mean_estimate_density
        mean_N <- summary_data$mean_N
        True_Density <- summary_data$true_density 
        Mean_relative_Estimate <- summary_data$relative_mean_estimate
        Percent_Bias <- summary_data$percent_bias
        RRMSE <- summary_data$rrmse
        CI_Coverage_Prob <- summary_data$ci_coverage_prob
        Mean_SE <- summary_data$mean_rse
        CV <- summary_data$sd_of_means
        Mean_ER <- summary_data$mean_ER
        Mean_se_ER <- summary_data$mean_se_ER
        N <- summary_data$N
        SampleBlocks <- summary_data$sampleblocks
        Effort <- summary_data$effortKM
    } else {
        Mean_estimated_Density <- NA
        mean_N <- NA
        True_Density <- NA
        Mean_relative_Estimate <- NA
        Percent_Bias <- NA
        RRMSE <- NA
        CI_Coverage_Prob <- NA
        Mean_SE <- NA
        CV <- NA
        Mean_ER <- NA
        Mean_se_ER <- NA
        N <- NA
        SampleBlocks <- NA
        Effort <- NA
    }
    data.frame(subplots = SampleBlocks, Mean_estimated_Density = Mean_estimated_Density, mean_N = mean_N, trueDensity = True_Density, Mean_relative_Estimate = Mean_relative_Estimate, Percent_Bias = Percent_Bias, RRMSE = RRMSE, CI_Coverage_Prob = CI_Coverage_Prob, Mean_SE = Mean_SE, CV = CV, Mean_ER = Mean_ER, Mean_se_ER = Mean_se_ER, N = N, NSubplots = simConfigs[1, 2], Abundance = simConfigs[1, 3], EffortKM = Effort, Overlap = simConfigs[1, 1])
})

# Combine the list of data frames into a single data frame
sum_stats_df <- do.call(rbind, stats)

# Convert all columns in extracted_df to numeric
sum_stats_df[, 1:8] <- lapply(sum_stats_df[, 1:8], as.numeric)


# Scale the area column so that 0 remains 0 and other values are scaled between 0 and 1
# max_area <- max(sum_stats_df$area, na.rm = TRUE)
# sum_stats_df$area_p <- round(sum_stats_df$area / max_area, 2)

# Print the resulting data frame to verify the changes
kable(sum_stats_df)


# Write the dataframe to a CSV file
output_path <- here("Output", "Simulation", "sumSimQC.csv")

# Check if sum_stats_df is not empty before writing to CSV
if (nrow(sum_stats_df) > 0) {
    write.csv(sum_stats_df, file = output_path, row.names = FALSE)
    print(paste0("Data written to ", output_path))
} else {
    message("sum_stats_df is empty. No CSV file was written.")
}