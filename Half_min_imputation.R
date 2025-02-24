# Load required library
library(tidyverse)
library(ggplot2)

# Part 1 ------
# Perform Half-min imputation

#load data 
final_data_2023 <- read.csv("/Users/marcinebessire/Desktop/project/Final_Data_2023.csv", check.names = FALSE)
final_data_2024 <- read.csv("/Users/marcinebessire/Desktop/project/Final_Data_2024.csv", check.names = FALSE)

#function for Half-minimum imputation
half_min_imputation <- function(df,year){
  #create copy
  df_imputed <- df 

  for (col in names(df_imputed)){ 
    if (is.numeric(df_imputed[[col]])){ #check if column is numeric 
      min_val <- min(df_imputed[[col]], na.rm = TRUE) #find min value 
      df_imputed[[col]][is.na(df_imputed[[col]])] <- 0.5 * min_val #calculate half min and replace NA
    }
  }
  #save new imputed dataset with same structure as original
  output_file <- paste0("/Users/marcinebessire/Desktop/project/Half_Min_", year, ".csv")
  write_csv(df_imputed,output_file)
  
  return(df_imputed)
}

#apply function to data 
imputed_data_23 <- half_min_imputation(final_data_2023, "2023")
imputed_data_24 <- half_min_imputation(final_data_2024, "2024")


# Part 2 ----
# Make new df with the same column names (same metabolites for comparison)
#columns to exclude
exclude_cols <- c('Name', 'ID', 'Year', 'MonthDay', 'Trial')

#identify same columns
cols_df1 <- setdiff(colnames(imputed_data_23), exclude_cols)
cols_df2 <- setdiff(colnames(imputed_data_24), exclude_cols)

#find common columns
common_cols <- intersect(cols_df1, cols_df2)

#visualize results
cat("Number of common columns:", length(common_cols), "\n") #84
cat("Common columns:", paste(common_cols, collapse=", ")) 

#make new dataframes for each year with those 84 columns + metadata
final_col <- c(exclude_cols, common_cols)

#new dataframes with same metabolites
half_min_23 <- imputed_data_23 %>% select(all_of(final_col))
half_min_24 <- imputed_data_24 %>% select(all_of(final_col))

# Part 3
# Wilcoxon rank-sum test (for independent data)

#select only metabolite columns (starting from column 6)
half_min_23_metabolites <- half_min_23[, 6:ncol(half_min_23)]
half_min_24_metabolites <- half_min_24[, 6:ncol(half_min_24)]

#get metabolite column names
metabolite_cols <- colnames(half_min_23_metabolites)

#dataframe for the results
results_Wilcoxon <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)

#loop through each metabolite and run Wilcoxon rank-sum test
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  test_result <- wilcox.test(half_min_23_metabolites[[metabolite]], half_min_24_metabolites[[metabolite]], paired = FALSE, exact = FALSE)
  
  #save results to dataframe
  results_Wilcoxon$p_value[i] <- test_result$p.value
  results_Wilcoxon$statistic[i] <- test_result$statistic
  
}

#adjust p-value for multiple testing (Benjamini-Hochberg
results_Wilcoxon$adj_p_value <- p.adjust(results_Wilcoxon$p_value, method = "BH")

print(results_Wilcoxon)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_Wilcoxon <- results_Wilcoxon %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_Wilcoxon)

# Part 4 
# Run unpaired t-test for each metabolite 
#dataframe for the results
results_ttest <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)

#loop through each metabolite and run t-test for each metabolites
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  #use Welchs test here (not same variance, to test H0 where two groups have equal mean)
  test_result2 <- t.test(half_min_23_metabolites[[metabolite]], half_min_24_metabolites[[metabolite]], paired = FALSE, var.equal = FALSE)
  
  #save results to dataframe
  results_ttest$p_value[i] <- test_result2$p.value
  results_ttest$statistic[i] <- test_result2$statistic
}

#adjust p-value for multiple testing (Benjamini-Hochberg)
results_ttest$adj_p_value <- p.adjust(results_ttest$p_value, method = "BH")

print(results_ttest)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_ttest <- results_ttest %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_ttest)

# Part 3 ----
# Check Half-Min imputation 
# by comparing the  mean median and standard deviation for and after imputation

# #identify columns with missing values for 2023 and 2024
# cols_with_na_23 <- colnames(final_data_2023)[colSums(is.na(final_data_2023)) > 0]
# cols_with_na_24 <- colnames(final_data_2024)[colSums(is.na(final_data_2024)) > 0]
# 
# #function to compute summary statistics
# compute_stats <- function(original_df, imputed_df, imputed_columns){
#   #original data
#   stats_original <- original_df %>%
#     select(all_of(imputed_columns)) %>%
#     summarise(across(everything(), list(
#       Median = ~median(.x, na.rm = TRUE),
#       Mean = ~mean(.x, na.rm = TRUE),
#       SD = ~sd(.x, na.rm = TRUE)
#     ))) %>%
#     mutate(Dataset = "Original") %>%
#     relocate(Dataset)
#   
#   #imputed data
#   stats_imputed <- imputed_df %>%
#     select(all_of(imputed_columns)) %>%
#     summarise(across(everything(), list(
#       Median = ~median(.x, na.rm = TRUE),
#       Mean = ~mean(.x, na.rm = TRUE),
#       SD = ~sd(.x, na.rm = TRUE)
#     ))) %>%
#     mutate(Dataset = "Imputed") %>%
#     relocate(Dataset)
#   
#   #calculate % change ((imputed-original)/orignal) * 100
#   stats_change <- stats_original
#   stats_change[,-1] <- ((stats_imputed[,-1] - stats_original[,-1])/stats_original[,-1]) * 100
#   stats_change$Dataset <- "% Change"
#   
#   final_stats <- bind_rows(stats_original, stats_imputed, stats_change)
#   
#   return(final_stats)
#   
# }
# 
# #compute statistics
# summary_comparison23 <- compute_stats(final_data_2023, imputed_data_23, cols_with_na_23)
# summary_comparison24 <- compute_stats(final_data_2024, imputed_data_24, cols_with_na_24)

# #save csv 
# output_file23 <- "/Users/marcinebessire/Desktop/project/Summary_Statistics_HalfMin_23.csv"
# write_csv(summary_comparison23, output_file23)
# output_file24 <- "/Users/marcinebessire/Desktop/project/Summary_Statistics_HalfMin_24.csv"
# write_csv(summary_comparison24, output_file24)
# 
# # Part 4 -----
# # Visualization of original and imputed data 
# 
# #for 2023
# #reshape data for plotting
# plot_data <- final_data_2023 %>%
#   select(all_of(cols_with_na_23)) %>%
#   pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
#   mutate(Dataset = "Original") %>%
#   bind_rows(
#     imputed_data_23 %>%
#       select(all_of(cols_with_na_23)) %>%
#       pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
#       mutate(Dataset = "Imputed")
#   )
# 
# # Boxplot Comparison
# ggplot(plot_data, aes(x = Dataset, y = Value, fill = Dataset)) +
#   geom_boxplot() +
#   facet_wrap(~Variable, scales = "free") +
#   theme_minimal() +
#   labs(title = "Boxplot: Original vs Imputed", x = "Dataset", y = "Value")
# 
# #for 2024
# #reshape data for plotting
# plot_data <- final_data_2024 %>%
#   select(all_of(cols_with_na_24)) %>%
#   pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
#   mutate(Dataset = "Original") %>%
#   bind_rows(
#     imputed_data_24 %>%
#       select(all_of(cols_with_na_24)) %>%
#       pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
#       mutate(Dataset = "Imputed")
#   )
# 
# # Boxplot Comparison
# ggplot(plot_data, aes(x = Dataset, y = Value, fill = Dataset)) +
#   geom_boxplot() +
#   facet_wrap(~Variable, scales = "free") +
#   theme_minimal() +
#   labs(title = "Boxplot: Original vs Imputed", x = "Dataset", y = "Value")
# 
# 
