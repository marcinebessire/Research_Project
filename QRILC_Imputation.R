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

significant_Wilcoxon_QRILC <- adjust_p_values(results_Wilcoxon_QRILC) #80
#significant_Wilcoxon_QRILC2 <- adjust_p_values(results_Wilcoxon_QRILC2) #80

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
significant_t_test_QRILC <- adjust_p_values(results_t_test_QRILC) #80
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
non_normal_count23 #59 metabolites are non-normal distributed 
non_normal_count24 <- sum(shapiro_df24$p_value < 0.05)
non_normal_count24 #46 metabolites are non-normal distributed 


# Part 4 -----
# Distirbution plot before and after Imputation and Imputed Values Only 

#remove metadata from whole data to get numeric data 
numeric_QRILC23 <- data23[, 6:ncol(data23)]
numeric_QRILC24 <- data24[, 6:ncol(data24)]

plot_imputation_distribution <- function(original_data, imputed_data, year, output_file) {
  #convert data to long format
  imputed_long <- imputed_data %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Value")
  
  original_long <- original_data %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Value")
  
  #identify imputed values
  imputed_only <- original_long %>%
    mutate(Imputed = is.na(Original_Value)) %>%
    filter(Imputed) %>%
    select(Metabolite) %>%
    inner_join(imputed_long, by = "Metabolite") %>%
    mutate(Dataset = "Imputed_Only")
  
  #merge original and imputed datasets
  comparison <- original_long %>%
    left_join(imputed_long, by = "Metabolite") %>%
    pivot_longer(cols = c("Original_Value", "Imputed_Value"), 
                 names_to = "Dataset", values_to = "Value")
  
  #add Imputed_Only as a separate dataset
  imputed_only <- imputed_only %>%
    mutate(Value = Imputed_Value, Dataset = "Imputed_Only") %>%
    select(Metabolite, Dataset, Value)
  
  #combine both datasets
  comparison <- bind_rows(comparison, imputed_only)
  
  # Open a PDF device to save the plot
  pdf(output_file, width = 8, height = 6)
  
  #generate the plot
  p <- ggplot(comparison, aes(x = Value, fill = Dataset)) +
    geom_density(alpha = 0.5) +  # Transparency for overlapping
    labs(title = paste("Distribution of Original, Imputed, and Imputed_Only Values (", year, ")", sep = ""),
         x = "Metabolite Value",
         y = "Density") +
    theme_minimal() + 
    scale_fill_manual(values = c("Original_Value" = "lightblue", 
                                 "Imputed_Value" = "red", 
                                 "Imputed_Only" = "green")) +
    xlim(-10, 50)
  
  print(p)
  dev.off()
}

density_plot_QRILC_23 <- plot_imputation_distribution(numeric_QRILC23, QRILC23, "2023", "/Users/marcinebessire/Desktop/project/QRILC_Distribution_Comparison23.pdf")
density_plot_QRILC_24 <- plot_imputation_distribution(numeric_QRILC24, QRILC24, "2024", "/Users/marcinebessire/Desktop/project/QRILC_Distribution_Comparison24.pdf")

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
    labs(title = paste("Density Plot of Normalized Difference (", year, ")", sep = ""),
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.4, 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")  # Reference line at 0
  
  print(p2)
  
  dev.off()
}

normalized_difference23 <- calculate_normalized_difference(numeric_QRILC23, QRILC23, "2023", "/Users/marcinebessire/Desktop/project/QRILC_normalized_difference23.pdf")
normalized_difference24 <- calculate_normalized_difference(numeric_QRILC24, QRILC24, "2024", "/Users/marcinebessire/Desktop/project/QRILC_normalized_difference24.pdf")


# # Extra part with correlation ----- 
# # Check if Impuation method was good 
# 
# # Correlation analysis 
# cor_before_23 <- cor(numeric23, use = "pairwise.complete.obs")
# cor_after_23 <- cor(QRILC23)
# cor_before_24 <- cor(numeric24, use = "pairwise.complete.obs")
# cor_after_24 <- cor(QRILC24)
# 
# # Heat map
# par(mfrow = c(1,2)) #for side by side plots
# corrplot(cor_before_23, method = "color", tl.cex = 0.6, title = "Before Imputation (2023)")
# corrplot(cor_after_23, method = "color", tl.cex = 0.6, title = "After Imputation (2023)")
# dev.off()
# 
# # Scatter plot 
# #identify missing vlaues 
# missing_val <- is.na(numeric23) #total of 99 are missing
# missing_val2 <- is.na(numeric24)
# 
# #create logical matrix identifying if a pair of metabolites had missing values
# missing_pairs <- (missing_val %*% t(missing_val)) > 0 #True 
# missing_pairs2 <- (missing_val2 %*% t(missing_val2)) > 0 #True
# 
# #extract upper triangular part of matrix 2023
# imputed_indices <- missing_pairs[upper.tri(missing_pairs, diag = FALSE)]
# cor_values_before23 <- cor_before_23[upper.tri(cor_before_23, diag = FALSE)]
# cor_values_after23 <- cor_after_23[upper.tri(cor_after_23, diag = FALSE)]
# 
# plot(cor_values_before23, cor_values_after23,
#      xlab = "Before Imputation",
#      ylab = "After Imputation",
#      main = "Correlation Comparison",
#      col = ifelse(imputed_indices, "red", "blue"), #red for imputed pairs
#      pch = 19) 
# 
# abline(0,1,col = "black", lwd = 2) #ideally points should line on line
# 
# #extract upper triangular part of matrix 2024
# imputed_indices2 <- missing_pairs2[upper.tri(missing_pairs2, diag = FALSE)]
# cor_values_before24 <- cor_before_24[upper.tri(cor_before_24, diag = FALSE)]
# cor_values_after24 <- cor_after_24[upper.tri(cor_after_24, diag = FALSE)]
# 
# plot(cor_values_before24, cor_values_after24,
#      xlab = "Before Imputation",
#      ylab = "After Imputation",
#      main = "Correlation Comparison",
#      col = ifelse(imputed_indices2, "red", "blue"), #red for imputed pairs
#      pch = 19) 
# 
# abline(0,1,col = "black", lwd = 2) #ideally points should line on line

