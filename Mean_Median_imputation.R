#load library
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

# Part 1 ------
# Mean imputation

#load data 
data23 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites23.csv", check.names = FALSE)
data24 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites24.csv", check.names = FALSE)

#numeric data 
numeric23 <- data23[, 6:ncol(data23)]
numeric24 <- data24[, 6:ncol(data24)]

#function for Mean imputation
mean_imputation <- function(df,year){
  #create copy
  df_imputed <- df 
  
  for (col in names(df_imputed)){ 
    if (is.numeric(df_imputed[[col]])){ #check if column is numeric 
      mean_val <- mean(df_imputed[[col]], na.rm = TRUE) #mean value of column
      df_imputed[[col]][is.na(df_imputed[[col]])] <- mean_val #calculate half min and replace NA
    }
  }
  
  #save new imputed dataset with same structure as original
  #output_file <- paste0("/Users/marcinebessire/Desktop/project/Half_Min_", year, ".csv")
  #write_csv(df_imputed,output_file
  return(df_imputed)
}

#apply function to data 
mean_imputed_data23 <- mean_imputation(numeric23, "2023")
mean_imputed_data24 <- mean_imputation(numeric24, "2024")

#function for Median imputation
median_imputation <- function(df,year){
  #create copy
  df_imputed <- df 
  
  for (col in names(df_imputed)){ 
    if (is.numeric(df_imputed[[col]])){ #check if column is numeric 
      median_val <- median(df_imputed[[col]], na.rm = TRUE) #mean value of column
      df_imputed[[col]][is.na(df_imputed[[col]])] <- median_val #calculate half min and replace NA
    }
  }
  
  #save new imputed dataset with same structure as original
  #output_file <- paste0("/Users/marcinebessire/Desktop/project/Half_Min_", year, ".csv")
  #write_csv(df_imputed,output_file
  return(df_imputed)
}

#apply function to data 
median_imputed_data23 <- median_imputation(numeric23, "2023")
median_imputed_data24 <- median_imputation(numeric24, "2024")

# Part 2 --------
# Wilcoxon rank-sum test (for independent data)

wilcoxon_test <- function(data1, data2) {
  #get metabolite column names
  metabolite_cols <- colnames(data1)
  
  #dataframe to store Wilcoxon test results
  results <- data.frame(
    Metabolite = metabolite_cols,
    p_value = numeric(length(metabolite_cols)),
    statistic = numeric(length(metabolite_cols))
  )
  
  #loop through each metabolite and run Wilcoxon rank-sum test
  for (i in seq_along(metabolite_cols)) {
    metabolite <- metabolite_cols[i]
    
    test_result <- wilcox.test(data1[[metabolite]], data2[[metabolite]], paired = FALSE, exact = FALSE)
    
    # Save results
    results$p_value[i] <- test_result$p.value
    results$statistic[i] <- test_result$statistic
  }
  
  return(results)
}

#call wilcoxon function
results_Wilcoxon_mean <- wilcoxon_test(mean_imputed_data23, mean_imputed_data24)
results_Wilcoxon_median <- wilcoxon_test(median_imputed_data23, median_imputed_data24)


#now adjust p-value and check for significance
adjust_p_values <- function(results, alpha = 0.05) {
  
  #adjust p-values using the Benjamini-Hochberg (BH) method
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  
  #filter significant results based on adjusted p-value threshold
  significant_results <- results %>%
    filter(adj_p_value < alpha) %>%
    arrange(adj_p_value)
  
  return(significant_results)
}

#Call significant function
significant_Wilcoxon_mean <- adjust_p_values(results_Wilcoxon_mean) #82/84
significant_Wilcoxon_median <- adjust_p_values(results_Wilcoxon_median) #82/84


# Part 3 -------
# Run unpaired t-test for each metabolite 

#dataframe for the results
t_test <- function(data1, data2) {
  #get metabolite column names
  metabolite_cols <- colnames(data1)
  
  #dataframe to store Wilcoxon test results
  results <- data.frame(
    Metabolite = metabolite_cols,
    p_value = numeric(length(metabolite_cols)),
    statistic = numeric(length(metabolite_cols))
  )
  
  #loop through each metabolite and run Wilcoxon rank-sum test
  for (i in seq_along(metabolite_cols)) {
    metabolite <- metabolite_cols[i]
    
    test_result <- t.test(data1[[metabolite]], data2[[metabolite]], paired = FALSE, var.equal = FALSE)
    
    # Save results
    results$p_value[i] <- test_result$p.value
    results$statistic[i] <- test_result$statistic
  }
  
  return(results)
}

#call t-test function
results_t_test_mean <- t_test(mean_imputed_data23, mean_imputed_data24)
results_t_test_median <- t_test(median_imputed_data23, median_imputed_data24)

#check significance 
significant_t_test_mean <- adjust_p_values(results_t_test_mean) #80/84
significant_t_test_median <- adjust_p_values(results_t_test_median) #79/84

# Part 3 -------
# Check normality of data 

#check normality using shapiro.test (becuase t-test assumes normality)
shapiro_results_mean23 <- apply(mean_imputed_data23, 2, function(x) shapiro.test(x)$p.value)
shapiro_results_mean24 <- apply(mean_imputed_data24, 2, function(x) shapiro.test(x)$p.value)
shapiro_results_median23 <- apply(median_imputed_data23, 2, function(x) shapiro.test(x)$p.value)
shapiro_results_median24 <- apply(median_imputed_data24, 2, function(x) shapiro.test(x)$p.value)

#convert to a dataframe for easy viewing
shapiro_df_mean23 <- data.frame(Metabolite = names(shapiro_results_mean23), p_value = shapiro_results_mean23) #total 84 metabolites
shapiro_df_mean24 <- data.frame(Metabolite = names(shapiro_results_mean24), p_value = shapiro_results_mean24)
shapiro_df_median23 <- data.frame(Metabolite = names(shapiro_results_median23), p_value = shapiro_results_median23) #total 84 metabolites
shapiro_df_median24 <- data.frame(Metabolite = names(shapiro_results_median24), p_value = shapiro_results_median24)

#if p-value < 0.05 then not normal distribution
non_normal_count_mean23 <- sum(shapiro_df_mean23$p_value < 0.05)
non_normal_count_mean23 #60 metabolites are non-normal distributed 
non_normal_count_mean24 <- sum(shapiro_df_mean24$p_value < 0.05)
non_normal_count_mean24 #60 metabolites are non-normal distributed 28 with CV

non_normal_count_median23 <- sum(shapiro_df_mean23$p_value < 0.05)
non_normal_count_median23 #60 metabolites are non-normal distributed 
non_normal_count_median24 <- sum(shapiro_df_mean24$p_value < 0.05)
non_normal_count_median24 #60 metabolites are non-normal distributed 28 with CV

# Part 4 -----
# Distirbution plot before and after Imputation and Imputed Values Only 

plot_imputation_distribution <- function(original_data, imputed_data, year, output_file) {
  #convert data to long format
  imputed_long <- imputed_data %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")
  
  original_long <- original_data %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")
  
  #identify imputed values
  imputed_only <- original_long %>%
    mutate(Imputed = is.na(Original_Data)) %>%
    filter(Imputed) %>%
    select(Metabolite) %>%
    inner_join(imputed_long, by = "Metabolite") %>%
    mutate(Dataset = "Imputed_Values")
  
  #merge original and imputed datasets
  comparison <- original_long %>%
    left_join(imputed_long, by = "Metabolite") %>%
    pivot_longer(cols = c("Original_Data", "Imputed_Data"), 
                 names_to = "Dataset", values_to = "Value")
  
  #add Imputed_Only as a separate dataset
  imputed_only <- imputed_only %>%
    mutate(Value = Imputed_Data, Dataset = "Imputed_Values") %>%
    select(Metabolite, Dataset, Value)
  
  #combine both datasets
  comparison <- bind_rows(comparison, imputed_only)
  
  # Open a PDF device to save the plot
  pdf(output_file, width = 8, height = 6)
  
  #generate the plot
  p <- ggplot(comparison, aes(x = Value, fill = Dataset)) +
    geom_density(alpha = 0.5) +  # Transparency for overlapping
    labs(title = paste("Distribution of Original Data, Imputed Data, and Imputed Values (", year, ")", sep = ""),
         x = "Metabolite Value",
         y = "Density") +
    theme_minimal() + 
    scale_fill_manual(values = c("Original_Data" = "lightblue", 
                                 "Imputed_Data" = "red", 
                                 "Imputed_Values" = "green")) +
    xlim(-10, 50)
  
  print(p)
  dev.off()
}

density_plot_mean_23 <- plot_imputation_distribution(numeric23, mean_imputed_data23, "2023", "/Users/marcinebessire/Desktop/project/Mean_Distribution_Comparison23.pdf")
density_plot_mean_24 <- plot_imputation_distribution(numeric24, mean_imputed_data24, "2024", "/Users/marcinebessire/Desktop/project/Mean_Distribution_Comparison24.pdf")
density_plot_median_23 <- plot_imputation_distribution(numeric23, median_imputed_data23, "2023", "/Users/marcinebessire/Desktop/project/Median_Distribution_Comparison23.pdf")
density_plot_median_24 <- plot_imputation_distribution(numeric24, median_imputed_data24, "2024", "/Users/marcinebessire/Desktop/project/Median_Distribution_Comparison24.pdf")


# Part 5 ------
# calculate normalized difference of each imputation (before and after) and plot

calculate_normalized_difference <- function(original_data, imputed_data, year, output_file) {
  #mean before imputation
  mean_before <- original_data %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_Before")
  
  #mean after imputation
  mean_after <- imputed_data %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_After")
  
  #merge before and after mean values
  mean_comparison <- left_join(mean_before, mean_after, by = "Metabolite")
  
  #normalized difference
  mean_comparison <- mean_comparison %>%
    mutate(Normalized_Difference = (Mean_After - Mean_Before) / Mean_Before)
  
  # save the plot
  pdf(output_file, width = 10, height = 6)
  
  #plot normalized difference
  p1 <- ggplot(mean_comparison, aes(x = Metabolite, y = Normalized_Difference, fill = Normalized_Difference)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = paste("Normalized Difference in Mean Before and After with Mean Imputation (", year, ")", sep = ""),
         x = "Metabolite",
         y = "Normalized Difference") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
  
  print(p1)
  
  #plot density of the normalized difference
  p2 <- ggplot(mean_comparison, aes(x = Normalized_Difference)) +
    geom_density(fill = "blue", alpha = 0.4, color = "black") +  # Density plot
    theme_minimal() +
    labs(title = paste("Density Plot of Normalized Difference with Mean Imputation (", year, ")", sep = ""),
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.4, 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")  # Reference line at 0
  
  print(p2)
  
  dev.off()
}

normalized_difference_mean23 <- calculate_normalized_difference(numeric23, mean_imputed_data23, "2023", "/Users/marcinebessire/Desktop/project/Mean_normalized_difference23.pdf")
normalized_difference_mean24 <- calculate_normalized_difference(numeric24, mean_imputed_data24, "2024", "/Users/marcinebessire/Desktop/project/Mean_normalized_difference24.pdf")

normalized_difference_median23 <- calculate_normalized_difference(numeric23, median_imputed_data23, "2023", "/Users/marcinebessire/Desktop/project/Median_normalized_difference23.pdf")
normalized_difference_median24 <- calculate_normalized_difference(numeric24, median_imputed_data24, "2024", "/Users/marcinebessire/Desktop/project/Median_normalized_difference24.pdf")


