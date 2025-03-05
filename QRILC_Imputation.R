#install packages required for QRILC
#install.packages("imputeLCMD")

#load necessary libraries
library(imputeLCMD)
library(tidyverse)
library(corrplot)

# Part 1 -------
# QRILC Imputation

#whole data 
data23 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites23.csv", check.names = FALSE)
data24 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites24.csv", check.names = FALSE)

#remove metadata from whole data to get numeric data 
numeric23 <- data23[, 6:ncol(data23)]
numeric24 <- data24[, 6:ncol(data24)]

#count missing values 
mv23 <- sum(is.na(numeric23)) #99 MV (26 with CV)
mv24 <- sum(is.na(numeric24)) #526 MV (333 with CV)

#total values 
tot23 <- nrow(numeric23) * ncol(numeric24) #6972
tot24 <- nrow(numeric24) * ncol(numeric24) #6048

#percentage of missing values 
percentage_mv23 <- (mv23 / tot23) * 100 #1.419% (0.522)
percentage_mv24 <- (mv24 / tot24) * 100 #8.697% (7.708)

#needs transformation to matrix 
numeric23 <- as.matrix(numeric23)
numeric24 <- as.matrix(numeric24)

#Log transformation (helps for very small values, prevents negative imputed values)
log_data23 <- log(numeric23 + 1e-6) #add small constant to avoid log(0)
log_data24 <- log(numeric24 + 1e-6) 

#QRILC imputation 
QRILC_impuation <- function(data, ...) {
  imputed_data <- impute.QRILC(data)[[1]] #returns list, extract imputed matrix
  return(imputed_data)
}

#call function for imputation 
imputed_QRILC23 <- QRILC_impuation(log_data23)
imputed_QRILC24 <- QRILC_impuation(log_data24)

#convert back to data frame and exponentiate back 
exp_QRILC23 <- exp(imputed_QRILC23) - 1e-6
exp_QRILC24 <- exp(imputed_QRILC24) - 1e-6
QRILC23 <- as.data.frame(exp_QRILC23)
QRILC24 <- as.data.frame(exp_QRILC24)

# #function for QRILC Imputation with Lower and Upper Bound
# QRILC_imputation_with_bounds <- function(data, lower_percentile = 0, upper_percentile = 0.99) {
#   #ensure data is a matrix
#   data_matrix <- as.matrix(data)
#   
#   #log transformation (helps prevent negative imputed values)
#   log_data <- log(data_matrix + 1e-6)  #small constant to avoid log(0)
#   
#   # Perform QRILC Imputation
#   imputed_log_data <- impute.QRILC(log_data)[[1]]  #extract matrix
#   
#   #convert back to original scale
#   imputed_data <- exp(imputed_log_data) - 1e-6
#   
#   #calculate lower and upper bounds based on observed (non-missing) values
#   lower_bound <- quantile(data, lower_percentile, na.rm = TRUE)  #min observed value
#   upper_bound <- quantile(data, upper_percentile, na.rm = TRUE)  #99th percentile
#   
#   #apply bounds to prevent extreme imputation
#   imputed_data[imputed_data < lower_bound] <- lower_bound
#   imputed_data[imputed_data > upper_bound] <- upper_bound
#   
#   #return as data frame
#   return(as.data.frame(imputed_data))
# }
# 
# # Apply the Function for QRILC Imputation with Bounds
# QRILC23_2 <- QRILC_imputation_with_bounds(numeric23)
# QRILC24_2 <- QRILC_imputation_with_bounds(numeric24)

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
results_Wilcoxon_QRILC <- wilcoxon_test(QRILC23, QRILC24)
#results_Wilcoxon_QRILC2 <- wilcoxon_test(QRILC23_2, QRILC24_2) 

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

significant_Wilcoxon_QRILC <- adjust_p_values(results_Wilcoxon_QRILC) #80 and 59 with CV
#significant_Wilcoxon_QRILC2 <- adjust_p_values(results_Wilcoxon_QRILC2) #81

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
results_t_test_QRILC <- t_test(QRILC23, QRILC24)
#results_t_test_QRILC2 <- t_test(QRILC23_2, QRILC24_2)

#check significance 
significant_t_test_QRILC <- adjust_p_values(results_t_test_QRILC) #80 and 58 with CV
#significant_t_test_QRILC2 <- adjust_p_values(results_t_test_QRILC2) #81

#grouped bar plot 
#count total number of metabolites
total_metabolites <- length(colnames(QRILC23)) #84

#nr of significant metabolites for each test
significant_wilcoxon <- nrow(significant_Wilcoxon_QRILC) #80
significant_ttest <- nrow(significant_t_test_QRILC) #80

#create summary dataframe
test_results <- data.frame(
  Test = rep(c("Wilcoxon", "T-test"), each = 2),
  Category = rep(c("Significant", "Non-Significant"), 2),
  Count = c(significant_wilcoxon, total_metabolites - significant_wilcoxon,
            significant_ttest, total_metabolites - significant_ttest)
)

pdf("/Users/marcinebessire/Desktop/project/QRILC_Significance.pdf", width = 10, height = 6)

#plot grouped bar chart
ggplot(test_results, aes(x = Test, y = Count, fill = Category)) +
  geom_bar(position="dodge", stat="identity") +  #grouped bars
  theme_minimal() +
  labs(title = "Number of Significant vs Non-Significant Metabolites using QRILC Imputation",
       x = "Statistical Test",
       y = "Count of Metabolites") +
  scale_fill_manual(values = c("Significant" = "blue", "Non-Significant" ="lightblue")) +
  theme(legend.position = "top") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5)  #add labels

dev.off()

# Part 3.2 -------
# Check normality of data 

#check normality using shapiro.test (becuase t-test assumes normality)
shapiro_results23 <- apply(QRILC23, 2, function(x) shapiro.test(x)$p.value)
shapiro_results24 <- apply(QRILC24, 2, function(x) shapiro.test(x)$p.value)


#convert to a dataframe for easy viewing
shapiro_df23 <- data.frame(Metabolite = names(shapiro_results23), p_value = shapiro_results23) #total 84 metabolites
shapiro_df24 <- data.frame(Metabolite = names(shapiro_results24), p_value = shapiro_results24)

#if p-value < 0.05 then not normal distribution
non_normal_count23 <- sum(shapiro_df23$p_value < 0.05)
non_normal_count23 #60 metabolites are non-normal distributed 43 with CV
non_normal_count24 <- sum(shapiro_df24$p_value < 0.05)
non_normal_count24 #46 metabolites are non-normal distributed 28 with CV


# Part 4 -----
# Distirbution plot before and after Imputation and Imputed Values Only 

#remove metadata from whole data to get numeric data 
numeric_23 <- data23[, 6:ncol(data23)]
numeric_24 <- data24[, 6:ncol(data24)]

#2023
#convert data to long format for visualization
imputed_23_long <- QRILC23 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")

original_23_long <- numeric_23 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")

#identify imputed values => missing values were replaced by half of the minimum observed value.
imputed_only_23 <- original_23_long %>%
  mutate(Imputed = is.na(Original_Data)) %>%
  filter(Imputed) %>%
  inner_join(imputed_23_long %>%
               distinct(Metabolite, .keep_all = TRUE), by = "Metabolite") %>%
  mutate(Dataset = "Imputed_Values")


#merge original and imputed datasets
comparison_23 <- original_23_long %>%
  left_join(imputed_23_long %>%
              distinct(Metabolite, .keep_all = TRUE), by = "Metabolite") %>%
  pivot_longer(cols = c("Original_Data", "Imputed_Data"), 
               names_to = "Dataset", values_to = "Value")

#add Imputed_Only as a separate dataset
imputed_only_23 <- imputed_only_23 %>%
  mutate(Value = Imputed_Data, Dataset = "Imputed_Values") %>%
  select(Metabolite, Dataset, Value)

#combine both datasets
comparison_23 <- bind_rows(comparison_23, imputed_only_23)

#2024
#convert data to long format for visualization
imputed_24_long <- QRILC24 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")

original_24_long <- numeric_24 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")

#identify imputed values => missing values were replaced by half of the minimum observed value.
imputed_only_24 <- original_24_long %>%
  mutate(Imputed = is.na(Original_Data)) %>%
  filter(Imputed) %>%
  inner_join(imputed_24_long %>%
               distinct(Metabolite, .keep_all = TRUE), by = "Metabolite") %>%
  mutate(Dataset = "Imputed_Values")

#merge original and imputed datasets
comparison_24 <- original_24_long %>%
  left_join(imputed_24_long %>%
              distinct(Metabolite, .keep_all = TRUE), by = "Metabolite") %>%
  pivot_longer(cols = c("Original_Data", "Imputed_Data"), 
               names_to = "Dataset", values_to = "Value")

#add Imputed_Only as a separate dataset
imputed_only_24 <- imputed_only_24 %>%
  mutate(Value = Imputed_Data, Dataset = "Imputed_Values") %>%
  select(Metabolite, Dataset, Value)

#combine both datasets
comparison_24 <- bind_rows(comparison_24, imputed_only_24)

#Now plot 
#open a PDF device to save multiple plots
pdf("/Users/marcinebessire/Desktop/project/QRILC_Distribution_Comparison.pdf", width = 8, height = 6)

#2023 plot 
ggplot(comparison_23, aes(x = Value, fill = Dataset)) +
  geom_density(alpha = 0.5) +  # Transparency for overlapping
  labs(title = "Distribution of Original/Imputed Data and Imputed Values with QRILC Imputation (2023)",
       x = "Metabolite Value",
       y = "Density") +
  theme_minimal() + 
  scale_fill_manual(values = c("Original_Data" = "lightblue", 
                               "Imputed_Data" = "red", 
                               "Imputed_Values" = "green")) +
  xlim(-10,50)

#2024 plot
#plot density distributions for original and imputed data separately
ggplot(comparison_24, aes(x = Value, fill = Dataset)) +
  geom_density(alpha = 0.5) +  # Transparency for overlapping
  labs(title = "Distribution of Original/Imputed Data, and Imputed Values with QRILC Imputation (2024)",
       x = "Metabolite Value",
       y = "Density") +
  theme_minimal() + 
  scale_fill_manual(values = c("Original_Data" = "lightblue", 
                               "Imputed_Data" = "red", 
                               "Imputed_Values" = "green")) +
  xlim(-10,50)

dev.off()


# Part 4.2 -----
# QQ plot comparing distribution

#2023
#convert data to long format for visualization
imputed_23_long <- QRILC23 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")

original_23_long <- numeric_23 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")

#identify imputed values => missing values were replaced by half of the minimum observed value.
imputed_only_23 <- original_23_long %>%
  mutate(Imputed = is.na(Original_Data)) %>%
  filter(Imputed) %>%
  inner_join(imputed_23_long %>%
               distinct(Metabolite, .keep_all = TRUE), by = "Metabolite") %>%
  mutate(Dataset = "Imputed_Values")


#merge original and imputed datasets
comparison_23 <- original_23_long %>%
  left_join(imputed_23_long %>%
              distinct(Metabolite, .keep_all = TRUE), by = "Metabolite") %>%
  pivot_longer(cols = c("Original_Data", "Imputed_Data"), 
               names_to = "Dataset", values_to = "Value")

#add Imputed_Only as a separate dataset
imputed_only_23 <- imputed_only_23 %>%
  mutate(Value = Imputed_Data, Dataset = "Imputed_Values") %>%
  select(Metabolite, Dataset, Value)

#combine both datasets
comparison_23 <- bind_rows(comparison_23, imputed_only_23)

#2024
#convert data to long format for visualization
imputed_24_long <- QRILC24 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")

original_24_long <- numeric_24 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")

#identify imputed values => missing values were replaced by half of the minimum observed value.
imputed_only_24 <- original_24_long %>%
  mutate(Imputed = is.na(Original_Data)) %>%
  filter(Imputed) %>%
  inner_join(imputed_24_long %>%
               distinct(Metabolite, .keep_all = TRUE), by = "Metabolite") %>%
  mutate(Dataset = "Imputed_Values")

#merge original and imputed datasets
comparison_24 <- original_24_long %>%
  left_join(imputed_24_long %>%
              distinct(Metabolite, .keep_all = TRUE), by = "Metabolite") %>%
  pivot_longer(cols = c("Original_Data", "Imputed_Data"), 
               names_to = "Dataset", values_to = "Value")

#add Imputed_Only as a separate dataset
imputed_only_24 <- imputed_only_24 %>%
  mutate(Value = Imputed_Data, Dataset = "Imputed_Values") %>%
  select(Metabolite, Dataset, Value)

#combine both datasets
comparison_24 <- bind_rows(comparison_24, imputed_only_24)


pdf("/Users/marcinebessire/Desktop/project/QRILC_QQ_Plots.pdf", width = 8, height = 6)

#QQ Plot Function
qq_plot <- function(data_x, data_y, x_label, y_label, title) {
  df <- data.frame(x = quantile(data_x, probs = seq(0, 1, 0.01), na.rm = TRUE),
                   y = quantile(data_y, probs = seq(0, 1, 0.01), na.rm = TRUE))
  
  ggplot(df, aes(x = x, y = y)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal()
}

#QQ Plot: Imputed Data 2023 vs. Original Data 2023
print(qq_plot(comparison_23$Value[comparison_23$Dataset == "Imputed_Data"], 
              comparison_23$Value[comparison_23$Dataset == "Original_Data"], 
              "Original Data (2023)", "Imputed Data (2023)", 
              "QQ Plot: Imputed Data vs. Original Data (2023)"))

#QQ Plot: Imputed Data 2024 vs. Original Data 2024
print(qq_plot(comparison_24$Value[comparison_24$Dataset == "Imputed_Data"], 
              comparison_24$Value[comparison_24$Dataset == "Original_Data"], 
              "Original Data (2024)", "Imputed Data (2024)", 
              "QQ Plot: Imputed Data vs. Original Data (2024)"))

#QQ Plot: Imputed Data 2024 vs. Imputed Data 2023
print(qq_plot(comparison_24$Value[comparison_24$Dataset == "Imputed_Data"], 
              comparison_23$Value[comparison_23$Dataset == "Imputed_Data"], 
              "Imputed Data (2023)", "Imputed Data (2024)", 
              "QQ Plot: Imputed Data 2024 vs. Imputed Data 2023"))

#QQ Plot: Original Data 2024 vs. Original Data 2023
print(qq_plot(comparison_24$Value[comparison_24$Dataset == "Original_Data"], 
              comparison_23$Value[comparison_23$Dataset == "Original_Data"], 
              "Original Data (2023)", "Original Data (2024)", 
              "QQ Plot: Original Data 2024 vs. Original Data 2023"))

dev.off()

#remove outliers using IQR methode (interquartile range to filter out extreme values)
remove_outliers <- function(data) {
  Q1 <- quantile(data, 0.25, na.rm = TRUE)
  Q3 <- quantile(data, 0.75, na.rm = TRUE)
  IQR_value <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_value
  upper_bound <- Q3 + 1.5 * IQR_value
  return(data[data >= lower_bound & data <= upper_bound])
}

#remove the outliers from the data
imputed_23_clean <- remove_outliers(comparison_23$Value[comparison_23$Dataset == "Imputed_Data"])
original_23_clean <- remove_outliers(comparison_23$Value[comparison_23$Dataset == "Original_Data"])

imputed_24_clean <- remove_outliers(comparison_24$Value[comparison_24$Dataset == "Imputed_Data"])
original_24_clean <- remove_outliers(comparison_24$Value[comparison_24$Dataset == "Original_Data"])

#pdf to save the QQ plots
pdf("/Users/marcinebessire/Desktop/project/QRILC_QQ_Plots_NoOutliers.pdf", width = 8, height = 6)

#QQ Plot: Imputed Data 2023 vs. Original Data 2023
print(qq_plot(original_23_clean, imputed_23_clean, 
              "Original Data (2023)", "Imputed Data (2023)", 
              "QQ Plot: Imputed Data vs. Original Data (2023)"))

# QQ Plot: Imputed Data 2024 vs. Original Data 2024
print(qq_plot(original_24_clean, imputed_24_clean, 
              "Original Data (2024)", "Imputed Data (2024)", 
              "QQ Plot: Imputed Data vs. Original Data (2024)"))

# QQ Plot: Imputed Data 2024 vs. Imputed Data 2023
print(qq_plot(imputed_23_clean, imputed_24_clean, 
              "Imputed Data (2023)", "Imputed Data (2024)", 
              "QQ Plot: Imputed Data 2024 vs. Imputed Data 2023"))

# QQ Plot: Original Data 2024 vs. Original Data 2023
print(qq_plot(original_23_clean, original_24_clean, 
              "Original Data (2023)", "Original Data (2024)", 
              "QQ Plot: Original Data 2024 vs. Original Data 2023"))

dev.off()

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
    labs(title = paste("Normalized Difference in Mean Before and After Imputation (", year, ")", sep = ""),
         x = "Metabolite",
         y = "Normalized Difference") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
  
  print(p1)
  
  #plot density of the normalized difference
  p2 <- ggplot(mean_comparison, aes(x = Normalized_Difference)) +
    geom_density(fill = "blue", alpha = 0.4, color = "black") +  # Density plot
    theme_minimal() +
    labs(title = paste("Density Plot of Normalized Difference with QRILC Imputation (", year, ")", sep = ""),
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.4, 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")  # Reference line at 0
  
  print(p2)
  
  dev.off()
}

normalized_difference23 <- calculate_normalized_difference(numeric_23, QRILC23, "2023", "/Users/marcinebessire/Desktop/project/QRILC_normalized_difference23.pdf")
normalized_difference24 <- calculate_normalized_difference(numeric_24, QRILC24, "2024", "/Users/marcinebessire/Desktop/project/QRILC_normalized_difference24.pdf")


