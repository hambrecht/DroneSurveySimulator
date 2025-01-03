library(here)
library(knitr)
library(dplyr)


# Extract key metrics from each simulation summary
extract_metrics <- function(sim) {
  summary_data <- summary(sim, description.summary = FALSE)

  list(
    relative_mean_estimate = round(summary_data@individuals$D$mean.Estimate / summary_data@individuals$D$Truth, 3), # mean estimate of abundance across simulation
    percent_bias = round(summary_data@individuals$D$percent.bias, 3), # the percentage of bias in the estimates
    rrmse = round(summary_data@individuals$D$RMSE / summary_data@individuals$D$mean.Estimate, 3), # root mean squared error/no. successful reps
    ci_coverage_prob = round(summary_data@individuals$D$CI.coverage.prob, 3), # proportion of times the 95% confidence interval contained the true value.
    mean_rse = round(summary_data@individuals$D$mean.se / summary_data@individuals$D$mean.Estimate, 3), # the mean standard error of the estimates of abundance
    sd_of_means = round(summary_data@individuals$D$sd.of.means / summary_data@individuals$D$mean.Estimate, 3), # the standard deviation of the estimates
    mean_ER = round(summary_data@individuals$summary$mean.ER, 3), # mean standard error of the encounter rates cross simulation
    mean_se_ER = summary_data@individuals$summary$mean.se.ER  # standard deviation of the encounter rates across simulation
  )
}

inputFilePaths <- list.files(path = here("Output", "Simulation"), pattern = "^simulation.*\\.RData$", full.names = TRUE)
load(inputFilePaths[2])

H_SG_metric <- extract_metrics(H_SG_sim)
Rnd_metric <- extract_metrics(Rnd_sim)
Sys_metric <- extract_metrics(Sys_sim)
ZZ_metric <- extract_metrics(ZZ_sim)
ZZC_metric <- extract_metrics(ZZC_sim)
FW_Sys_2C_metric <- extract_metrics(FW_Sys_2C_sim)
FW_ZZ_2C_metric <- extract_metrics(FW_ZZ_2C_sim)
FW_Sys_G_metric <- extract_metrics(FW_Sys_G_sim)
FW_ZZ_G_metric <- extract_metrics(FW_ZZ_G_sim)
QC_Sys_metric <- extract_metrics(QC_Sys_sim)


# Combine metrics into a single dataframe
comparison_df <- data.frame(
  Simulation = c("H-SG", "Rnd", "Sys", "ZZ", "ZZC", "FW-Sys_2C", "FW-ZZ_2C", "FW-Sys_G", "FW-ZZ_G", "QC-Sys"),
  # Relative_Avalibility = c(H_SG_metric$percent_available, Rnd_metric$percent_available, Sys_metric$percent_available, ZZ_metric$percent_available, ZZC_metric$percent_available, FW_Sys_2C_metric$percent_available, FW_ZZ_2C_metric$percent_available, FW_Sys_G_metric$percent_available, FW_ZZ_G_metric$percent_available, QC_Sys_metric$percent_available),
  Mean_relative_Estimate = c(H_SG_metric$relative_mean_estimate, Rnd_metric$relative_mean_estimate, Sys_metric$relative_mean_estimate, ZZ_metric$relative_mean_estimate, ZZC_metric$relative_mean_estimate, FW_Sys_2C_metric$relative_mean_estimate, FW_ZZ_2C_metric$relative_mean_estimate, FW_Sys_G_metric$relative_mean_estimate, FW_ZZ_G_metric$relative_mean_estimate, QC_Sys_metric$relative_mean_estimate),
  Percent_Bias = c(H_SG_metric$percent_bias, Rnd_metric$percent_bias, Sys_metric$percent_bias, ZZ_metric$percent_bias, ZZC_metric$percent_bias, FW_Sys_2C_metric$percent_bias, FW_ZZ_2C_metric$percent_bias, FW_Sys_G_metric$percent_bias, FW_ZZ_G_metric$percent_bias, QC_Sys_metric$percent_bias),
  RRMSE = c(H_SG_metric$rrmse, Rnd_metric$rrmse, Sys_metric$rrmse, ZZ_metric$rrmse, ZZC_metric$rrmse, FW_Sys_2C_metric$rrmse, FW_ZZ_2C_metric$rrmse, FW_Sys_G_metric$rrmse, FW_ZZ_G_metric$rrmse, QC_Sys_metric$rrmse),
  CI_Coverage_Prob = c(H_SG_metric$ci_coverage_pro, Rnd_metric$ci_coverage_pro, Sys_metric$ci_coverage_prob, ZZ_metric$ci_coverage_prob, ZZC_metric$ci_coverage_prob, FW_Sys_2C_metric$ci_coverage_prob, FW_ZZ_2C_metric$ci_coverage_prob, FW_Sys_G_metric$ci_coverage_prob, FW_ZZ_G_metric$ci_coverage_prob, QC_Sys_metric$ci_coverage_prob),
  Mean_SE = c(H_SG_metric$mean_rse, Rnd_metric$mean_rse, Sys_metric$mean_rse, ZZ_metric$mean_rse, ZZC_metric$mean_rse, FW_Sys_2C_metric$mean_rse, FW_ZZ_2C_metric$mean_rse, FW_Sys_G_metric$mean_rse, FW_ZZ_G_metric$mean_rse, QC_Sys_metric$mean_rse),
  CV = c(H_SG_metric$sd_of_means, Rnd_metric$sd_of_means, Sys_metric$sd_of_means, ZZ_metric$sd_of_means, ZZC_metric$sd_of_means, FW_Sys_2C_metric$sd_of_means, FW_ZZ_2C_metric$sd_of_means, FW_Sys_G_metric$sd_of_means, FW_ZZ_G_metric$sd_of_means, QC_Sys_metric$sd_of_means),
  Mean_n = c(H_SG_metric$mean_rN, Rnd_metric$mean_rN, Sys_metric$mean_rN, ZZ_metric$mean_rN, ZZC_metric$mean_rN, FW_Sys_2C_metric$mean_rN, FW_ZZ_2C_metric$mean_rN, FW_Sys_G_metric$mean_rN, FW_ZZ_G_metric$mean_rN, QC_Sys_metric$mean_rN),
  Mean_ER = c(H_SG_metric$mean_ER, Rnd_metric$mean_ER, Sys_metric$mean_ER, ZZ_metric$mean_ER, ZZC_metric$mean_ER, FW_Sys_2C_metric$mean_ER, FW_ZZ_2C_metric$mean_ER, FW_Sys_G_metric$mean_ER, FW_ZZ_G_metric$mean_ER, QC_Sys_metric$mean_ER),
  Mean_se_ER = c(H_SG_metric$mean_se_ER, Rnd_metric$mean_se_ER, Sys_metric$mean_se_ER, ZZ_metric$mean_se_ER, ZZC_metric$mean_se_ER, FW_Sys_2C_metric$mean_se_ER, FW_ZZ_2C_metric$mean_se_ER, FW_Sys_G_metric$mean_se_ER, FW_ZZ_G_metric$mean_se_ER, QC_Sys_metric$mean_se_ER)
)

# Print the comparison dataframe
# print(comparison_df)
kable(comparison_df)


# Initialize an empty list to store the metrics
metrics_list <- list()

# Loop through each file, load it, extract the metrics, and store them in the list
for (file in inputFilePaths) {
  load(file)
  sim_name <- sub(".*simulation-(.*)\\.RData", "\\1", file)
  H_SG_metric <- extract_metrics(H_SG_sim)
  Rnd_metric <- extract_metrics(Rnd_sim)
  Sys_metric <- extract_metrics(Sys_sim)
  ZZ_metric <- extract_metrics(ZZ_sim)
  ZZC_metric <- extract_metrics(ZZC_sim)
  FW_Sys_2C_metric <- extract_metrics(FW_Sys_2C_sim)
  FW_ZZ_2C_metric <- extract_metrics(FW_ZZ_2C_sim)
  FW_Sys_G_metric <- extract_metrics(FW_Sys_G_sim)
  FW_ZZ_G_metric <- extract_metrics(FW_ZZ_G_sim)
  QC_Sys_metric <- extract_metrics(QC_Sys_sim)

  # Combine metrics into a single dataframe
  metrics <- data.frame(
      Simulation = c("H-SG", "Rnd", "Sys", "ZZ", "ZZC", "FW-Sys_2C", "FW-ZZ_2C", "FW-Sys_G", "FW-ZZ_G", "QC-Sys"),
      # Relative_Avalibility = c(H_SG_metric$percent_available, Rnd_metric$percent_available, Sys_metric$percent_available, ZZ_metric$percent_available, ZZC_metric$percent_available, FW_Sys_2C_metric$percent_available, FW_ZZ_2C_metric$percent_available, FW_Sys_G_metric$percent_available, FW_ZZ_G_metric$percent_available, QC_Sys_metric$percent_available),
      Mean_relative_Estimate = c(H_SG_metric$relative_mean_estimate, Rnd_metric$relative_mean_estimate, Sys_metric$relative_mean_estimate, ZZ_metric$relative_mean_estimate, ZZC_metric$relative_mean_estimate, FW_Sys_2C_metric$relative_mean_estimate, FW_ZZ_2C_metric$relative_mean_estimate, FW_Sys_G_metric$relative_mean_estimate, FW_ZZ_G_metric$relative_mean_estimate, QC_Sys_metric$relative_mean_estimate),
      Percent_Bias = c(H_SG_metric$percent_bias, Rnd_metric$percent_bias, Sys_metric$percent_bias, ZZ_metric$percent_bias, ZZC_metric$percent_bias, FW_Sys_2C_metric$percent_bias, FW_ZZ_2C_metric$percent_bias, FW_Sys_G_metric$percent_bias, FW_ZZ_G_metric$percent_bias, QC_Sys_metric$percent_bias),
      RRMSE = c(H_SG_metric$rrmse, Rnd_metric$rrmse, Sys_metric$rrmse, ZZ_metric$rrmse, ZZC_metric$rrmse, FW_Sys_2C_metric$rrmse, FW_ZZ_2C_metric$rrmse, FW_Sys_G_metric$rrmse, FW_ZZ_G_metric$rrmse, QC_Sys_metric$rrmse),
      CI_Coverage_Prob = c(H_SG_metric$ci_coverage_pro, Rnd_metric$ci_coverage_pro, Sys_metric$ci_coverage_prob, ZZ_metric$ci_coverage_prob, ZZC_metric$ci_coverage_prob, FW_Sys_2C_metric$ci_coverage_prob, FW_ZZ_2C_metric$ci_coverage_prob, FW_Sys_G_metric$ci_coverage_prob, FW_ZZ_G_metric$ci_coverage_prob, QC_Sys_metric$ci_coverage_prob),
      Mean_SE = c(H_SG_metric$mean_rse, Rnd_metric$mean_rse, Sys_metric$mean_rse, ZZ_metric$mean_rse, ZZC_metric$mean_rse, FW_Sys_2C_metric$mean_rse, FW_ZZ_2C_metric$mean_rse, FW_Sys_G_metric$mean_rse, FW_ZZ_G_metric$mean_rse, QC_Sys_metric$mean_rse),
      CV = c(H_SG_metric$sd_of_means, Rnd_metric$sd_of_means, Sys_metric$sd_of_means, ZZ_metric$sd_of_means, ZZC_metric$sd_of_means, FW_Sys_2C_metric$sd_of_means, FW_ZZ_2C_metric$sd_of_means, FW_Sys_G_metric$sd_of_means, FW_ZZ_G_metric$sd_of_means, QC_Sys_metric$sd_of_means),
      Mean_n = c(H_SG_metric$mean_rN, Rnd_metric$mean_rN, Sys_metric$mean_rN, ZZ_metric$mean_rN, ZZC_metric$mean_rN, FW_Sys_2C_metric$mean_rN, FW_ZZ_2C_metric$mean_rN, FW_Sys_G_metric$mean_rN, FW_ZZ_G_metric$mean_rN, QC_Sys_metric$mean_rN),
      Mean_ER = c(H_SG_metric$mean_ER, Rnd_metric$mean_ER, Sys_metric$mean_ER, ZZ_metric$mean_ER, ZZC_metric$mean_ER, FW_Sys_2C_metric$mean_ER, FW_ZZ_2C_metric$mean_ER, FW_Sys_G_metric$mean_ER, FW_ZZ_G_metric$mean_ER, QC_Sys_metric$mean_ER),
      Mean_se_ER = c(H_SG_metric$mean_se_ER, Rnd_metric$mean_se_ER, Sys_metric$mean_se_ER, ZZ_metric$mean_se_ER, ZZC_metric$mean_se_ER, FW_Sys_2C_metric$mean_se_ER, FW_ZZ_2C_metric$mean_se_ER, FW_Sys_G_metric$mean_se_ER, FW_ZZ_G_metric$mean_se_ER, QC_Sys_metric$mean_se_ER)
    )
  metrics_list[[sim_name]] <- metrics
}

# Combine the metrics into a single dataframe
comparison_df <- do.call(rbind, lapply(metrics_list, as.data.frame))
kable(comparison_df)
output_path <- here("Output", "Simulation", "sim_comparison.csv")
write.csv(comparison_df, file = output_path, row.names = TRUE)

# Assuming comparison_df is your dataframe
averaged_df <- comparison_df %>%
  group_by(Simulation) %>%
  summarize(
    # Relative_Avalibility = round(mean(Relative_Avalibility, na.rm = TRUE), 3),
    Mean_relative_Estimate = round(mean(Mean_relative_Estimate, na.rm = TRUE), 3),
    Percent_Bias = round(mean(Percent_Bias, na.rm = TRUE), 3),
    RRMSE = round(mean(RRMSE, na.rm = TRUE), 3),
    CI_Coverage_Prob = round(mean(CI_Coverage_Prob, na.rm = TRUE), 3),
    relative_Mean_SE = round(mean(Mean_SE, na.rm = TRUE), 3),
    CV = round(mean(CV, na.rm = TRUE), 3),
    Mean_n = round(mean(Mean_n, na.rm = TRUE), 3),
    Mean_ER = round(mean(Mean_ER, na.rm = TRUE), 3),
    Mean_se_ER = round(mean(Mean_se_ER, na.rm = TRUE), 6)
  )
# Manually set the order of the groups
averaged_df$Simulation <- factor(averaged_df$Simulation, levels = c(
  "H-SG",
  "Sys",
  "Rnd",
  "ZZ",
  "ZZC",
  "FW-Sys_2C",
  "FW-ZZ_2C",
  "FW-Sys_G",
  "FW-ZZ_G",
  "QC-Sys"
))
# Sort the dataframe by the Simulation column
averaged_df <- averaged_df %>% arrange(Simulation)

# Transpose the dataframe
transposed_df <- as.data.frame(t(averaged_df))

# Set the column names to the first row and remove the first row
colnames(transposed_df) <- transposed_df[1,]
transposed_df <- transposed_df[-1,]

# Print the transposed dataframe
kable(transposed_df)

# Save the transposed dataframe as a CSV file
output_path <- here("Output", "Simulation", "sim_average.csv")
write.csv(transposed_df, file = output_path, row.names = TRUE)
