library(here)
library(knitr)
library(dplyr)


# Extract key metrics from each simulation summary
extract_metrics <- function(sim) {
  summary_data <- summary(sim, description.summary = FALSE)

  list(
    mean_estimate_density = round(summary_data@individuals$D$mean.Estimate*1e6, 3), # mean estimate of abundance across simulation
    true_density =  round(summary_data@individuals$D$Truth*1e6, 3),
    relative_mean_estimate_density = round(summary_data@individuals$D$mean.Estimate / summary_data@individuals$D$Truth, 3), # mean estimate of abundance across simulation
    percent_bias = round(summary_data@individuals$D$percent.bias, 3), # the percentage of bias in the estimates
    rrmse = round(summary_data@individuals$D$RMSE / summary_data@individuals$D$mean.Estimate, 3), # root mean squared error/no. successful reps
    ci_coverage_prob = round(summary_data@individuals$D$CI.coverage.prob, 3), # proportion of times the 95% confidence interval contained the true value.
    mean_rse = round(summary_data@individuals$D$mean.se / summary_data@individuals$D$mean.Estimate, 3), # the mean standard error of the estimates of abundance
    sd_of_means = round(summary_data@individuals$D$sd.of.means / summary_data@individuals$D$mean.Estimate, 3), # the standard deviation of the estimates
    mean_ER = summary_data@individuals$summary$mean.ER, # mean standard error of the encounter rates cross simulation
    mean_se_ER = summary_data@individuals$summary$mean.se.ER  # standard deviation of the encounter rates across simulation
  )
}

# List all files containing 'cover' in folder Output/Simulation
files <- list.files(path = here("Output", "Simulation"), pattern = "density_sim.*\\.RData$", full.names = TRUE)
inputFilePaths <- list.files(path = here("Output", "Simulation"), pattern = "^simulation.*\\.RData$", full.names = TRUE)
# Filter inputFilePaths to include only entries with '503' and '517'
filteredPaths <- grep("503|517", inputFilePaths, value = TRUE)

# Add filtered paths to files list
files <- c(files, filteredPaths)

# wmu_path_list <- wmu_path_list[-4]

# Extract the last three digits before the file extension for each file
wmu_ids <- sapply(files, function(file) {
    match <- regmatches(file, regexpr("(?<=.{15})\\d{3}", file, perl = TRUE))
    return(match)
})
# Extract one to three digits between two points (e.g., *.06.*, *.123.*) for each file
density_ids <- sapply(files, function(file) {
    match <- regmatches(file, regexpr("\\.\\d{1,3}\\.", file))
    return(sub("\\.", "", sub("\\.$", "", match)))  # Remove leading and trailing dots
})


# Initialize an empty list to store the metrics
metrics_list <- list()

# Loop through each file, load it, extract the metrics, and store them in the list
for (i in seq_along(files)) {
    before_load <- ls()
    load(files[i])
    w_suffix <- wmu_ids[[i]]
    if (length(density_ids[[i]]) == 0) {
      d_suffix <- "1"
    } else {
      d_suffix <- paste0("0.", density_ids[[i]])
    }
    sim_name <- paste0("WMU",w_suffix,"_D0.",d_suffix)
    H_SG_metric <- extract_metrics(H_SG_sim_density)
    FW_Sys_G_metric <- extract_metrics(FW_Sys_G_sim_density)
    QC_Sys_metric <- extract_metrics(QC_Sys_sim_density)

    # Combine metrics into a single dataframe
    metrics <- data.frame(
        WMU = rep(w_suffix,3),
        Density = rep(d_suffix,3),
        Simulation = c("H-SG", "FW-Sys_G", "QC-Sys"),
        Mean_estimated_Density = c(H_SG_metric$mean_estimate_density, FW_Sys_G_metric$mean_estimate_density, QC_Sys_metric$mean_estimate_density),
        True_Density = c(H_SG_metric$true_density,  FW_Sys_G_metric$true_density,  QC_Sys_metric$true_density),
        Mean_relative_Estimate = c(H_SG_metric$relative_mean_estimate, FW_Sys_G_metric$relative_mean_estimate, QC_Sys_metric$relative_mean_estimate),
        Percent_Bias = c(H_SG_metric$percent_bias, FW_Sys_G_metric$percent_bias,  QC_Sys_metric$percent_bias),
        RRMSE = c(H_SG_metric$rrmse, FW_Sys_G_metric$rrmse, QC_Sys_metric$rrmse),
        CI_Coverage_Prob = c(H_SG_metric$ci_coverage_pro, FW_Sys_G_metric$ci_coverage_prob,  QC_Sys_metric$ci_coverage_prob),
        Mean_SE = c(H_SG_metric$mean_rse,  FW_Sys_G_metric$mean_rse, QC_Sys_metric$mean_rse),
        CV = c(H_SG_metric$sd_of_means,  FW_Sys_G_metric$sd_of_means,  QC_Sys_metric$sd_of_means),
        Mean_ER = c(H_SG_metric$mean_ER, FW_Sys_G_metric$mean_ER,  QC_Sys_metric$mean_ER),
        Mean_se_ER = c(H_SG_metric$mean_se_ER, FW_Sys_G_metric$mean_se_ER, QC_Sys_metric$mean_se_ER)
        )
    metrics_list[[sim_name]] <- metrics
}

# Combine the metrics into a single dataframe
comparison_df <- do.call(rbind, lapply(metrics_list, as.data.frame))
kable(comparison_df)

output_path <- here("Output", "Simulation", "density_sim_comparison.csv")
write.csv(comparison_df, file = output_path, row.names = TRUE)
