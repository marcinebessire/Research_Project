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
non_normal_count23 #58 metabolites are non-normal distributed 43 with CV
non_normal_count24 <- sum(shapiro_df24$p_value < 0.05)
non_normal_count24 #47 metabolites are non-normal distributed 28 with CV


# Part 4 -----
# Distribution plot before and after Imputation and Imputed Values Only 

#change numeric data back
numeric23 <- data23[, 6:ncol(data23)]
numeric24 <- data24[, 6:ncol(data24)]

#function to plot distribution before and after imputation (entire dataset)
plot_whole_distribution <- function(original, imputed, method, year) {
  #numeric columns (starting from column 6)
  numeric_original <- original
  numeric_imputed <- imputed
  
  #count missing values in the original dataset
  missing_count <- sum(is.na(numeric_original))
  
  #convert to long format for plotting
  original_long <- numeric_original %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    mutate(Data = "Original Data")
  
  imputed_long <- numeric_imputed %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    mutate(Data = "Imputed Data")
  
  #identify imputed values (handling NA correctly)
  imputed_values <- is.na(numeric_original) & !is.na(numeric_imputed)
  
  imputed_only_long <- numeric_imputed %>%
    as.data.frame() %>%
    mutate(across(everything(), ~ replace(., !imputed_values[, cur_column()], NA))) %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    filter(!is.na(Value)) %>%
    mutate(Data = "Imputed Values")
  
  #count values in each long format dataset (for ensuring that it is correctly done)
  original_count <- nrow(original_long)
  imputed_count <- nrow(imputed_long)
  imputed_only_count <- nrow(imputed_only_long)
  
  print(paste("Total values in Original Data:", original_count))
  print(paste("Total values in Imputed Data:", imputed_count))
  print(paste("Total imputed values (only new values):", imputed_only_count))
  
  #combine all data
  combined_data <- bind_rows(original_long, imputed_long, imputed_only_long)
  
  #compute mean for each data category
  mean_data <- combined_data %>%
    group_by(Data) %>%
    summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  #plot overall density distribution
  plot <- ggplot(combined_data, aes(x = Value, fill = Data)) +
    geom_density(alpha = 0.5, na.rm = TRUE) + 
    theme_minimal() +
    labs(title = paste0("Overall Density Distribution Before and After Imputation (", method, " ", year, ")"),
         x = "Value",
         y = "Density") +
    geom_vline(data = mean_data %>% filter(Data == "Original Data"),
               aes(xintercept = mean_value), color = "blue", linewidth = 1, linetype = "dashed") + 
    geom_vline(data = mean_data %>% filter(Data == "Imputed Data"),
               aes(xintercept = mean_value), color = "red", linewidth = 1, linetype = "dashed") +
    xlim(-10, 100) + 
    theme(legend.position = "bottom")
  
  print(plot)
  return(combined_data)
}

#example function call
dist_QRILC23 <- plot_whole_distribution(numeric23, QRILC23, "QRILC", "2023")
dist_QRILC24 <- plot_whole_distribution(numeric24, QRILC24, "QRILC", "2024")


# Part 8 ------
# calculate normalized difference of each imputation (before and after) and plot

#2023
#mean before imputation
mean_before23 <- numeric23 %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_Before")

#means after imputation
mean_after23 <- QRILC23 %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_After")

#merge before and after mean values
mean_comparison23 <- left_join(mean_before23, mean_after23, by = "Metabolite")

#compute normalized difference: (Mean_After - Mean_Before) / Mean_Before
mean_comparison23 <- mean_comparison23 %>%
  mutate(Normalized_Difference = (Mean_After - Mean_Before) / Mean_Before)

#2024
#mean before imputation
mean_before24 <- numeric24 %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_Before")

#means after imputation
mean_after24 <- QRILC24 %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_After")

#merge before and after mean values
mean_comparison24 <- left_join(mean_before24, mean_after24, by = "Metabolite")

#compute normalized difference: (Mean_After - Mean_Before) / Mean_Before
mean_comparison24 <- mean_comparison24 %>%
  mutate(Normalized_Difference = (Mean_After - Mean_Before) / Mean_Before)


#plot normalized difference
ggplot(mean_comparison23, aes(x = Metabolite, y = Normalized_Difference, fill = Normalized_Difference)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Normalized Difference in Mean Before and After with QRILC Imputation (2023)",
       x = "Metabolite",
       y = "Normalized Difference") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)


#plot normalized difference
ggplot(mean_comparison24, aes(x = Metabolite, y = Normalized_Difference, fill = Normalized_Difference)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Normalized Difference in Mean Before and After with QRILC Imputation (2024)",
       x = "Metabolite",
       y = "Normalized Difference") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

#plot the density of the normalized difference
ggplot(mean_comparison23, aes(x = Normalized_Difference)) +
  geom_density(fill = "blue", alpha = 0.4, color = "black") +  #density plot
  theme_minimal() +
  labs(title = "Density Plot of Normalized Difference with QRILC Imputation (2023)",
       x = "Normalized Difference",
       y = "Density") +
  xlim(-0.2,0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")  #reference line at 0

#plot the density of the normalized difference
ggplot(mean_comparison24, aes(x = Normalized_Difference)) +
  geom_density(fill = "blue", alpha = 0.4, color = "black") +  #density plot
  theme_minimal() +
  labs(title = "Density Plot of Normalized Difference with QRILC Imputation (2024)",
       x = "Normalized Difference",
       y = "Density") +
  xlim(-0.4,0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")  #reference line at 0


# Part 9 -----
# Kolomogorov-Smirnov test: nonparametric test to test whether two ssamples came from same distirbution

#convert data to long format
#2023
imputed_23_long <- QRILC23 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")

original_23_long <- numeric23 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")

#2024
#convert data to long format
imputed_24_long <- QRILC24 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")

original_24_long <- numeric24 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")


#original 2023 vs imputed 2023
ks.test(original_23_long$Original_Data, imputed_23_long$Imputed_Data)
#original 2024 vs imputed 2024
ks.test(original_24_long$Original_Data, imputed_24_long$Imputed_Data)
#original 2023 vs Original 2024
ks.test(original_23_long$Original_Data, original_24_long$Original_Data)
#imputed 2023 vs imputed 2024
ks.test(imputed_23_long$Imputed_Data, imputed_24_long$Imputed_Data)

#make QQ plot
qqplot(original_23_long$Original_Data, imputed_23_long$Imputed_Data,
       main = "Q-Q Plot: Original vs Imputed Data (2023)",
       xlab = "Original Data Quantiles",
       ylab = "Imputed Data Quantiles",
       col = "blue", pch = 19)
abline(0, 1, col = "red", lwd = 2)  

qqplot(original_24_long$Original_Data, imputed_24_long$Imputed_Data,
       main = "Q-Q Plot: Original vs Imputed Data (2024)",
       xlab = "Original Data Quantiles",
       ylab = "Imputed Data Quantiles",
       col = "blue", pch = 19)
abline(0, 1, col = "red", lwd = 2)  


qqplot(original_23_long$Original_Data, original_24_long$Original_Data,
       main = "Q-Q Plot: 2023 vs 2024 Data (Original)",
       xlab = "Original Data Quantiles",
       ylab = "Imputed Data Quantiles",
       col = "blue", pch = 19)
abline(0, 1, col = "red", lwd = 2)  

qqplot(imputed_23_long$Imputed_Data, imputed_24_long$Imputed_Data,
       main = "Q-Q Plot: 2023 vs 2024 Data (Imputed)",
       xlab = "Original Data Quantiles",
       ylab = "Imputed Data Quantiles",
       col = "blue", pch = 19)
abline(0, 1, col = "red", lwd = 2)  

#save to file 
path23 <- paste0("/Users/marcinebessire/Desktop/project/QRILC_Imputation23.csv")
write_csv(imputed_23_long, path23)

path24 <- paste0("/Users/marcinebessire/Desktop/project/QRILC_Imputation24.csv")
write_csv(imputed_24_long, path24)

