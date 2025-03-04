#load library
library(impute) #for impute.knn() function 
library(corrplot) #for correlation coefficient calculation
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

# Part 1 -------
# Perform KNN imputation of Missing values 

#load data 
data23 <- read.csv("/Users/marcinebessire/Desktop/project/CV_Common_Metabolites23.csv", check.names = FALSE)
data24 <- read.csv("/Users/marcinebessire/Desktop/project/CV_Common_Metabolites24.csv", check.names = FALSE)

#numeric data 
numeric23 <- data23[, 6:ncol(data23)]
numeric24 <- data24[, 6:ncol(data24)]

#need to transpose because imput.knn() has samples in columns and metabolites in rows
#impute.knn() needs matrix format, because of mathematical operations (euclidean distances)

#2023
imputed_KNN23 <- impute.knn(as.matrix(t(numeric23)), rowmax = 0.5, colmax = 1) #transform first
imputed_KNN23 <- as.data.frame(t(imputed_KNN23$data)) #transform back and use transformed data 

#2024
imputed_KNN24 <- impute.knn(as.matrix(t(numeric24)), rowmax = 0.5, colmax = 1)
imputed_KNN24 <- as.data.frame(t(imputed_KNN24$data))

# Part 2 --------
# Wilcoxon rank-sum test (for independent data)

#get metabolite column names
metabolite_cols <- colnames(imputed_KNN23)

#dataframe for the results
results_Wilcoxon_KNN <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)

#loop through each metabolite and run Wilcoxon rank-sum test
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  test_result <- wilcox.test(imputed_KNN23[[metabolite]], imputed_KNN24[[metabolite]], paired = FALSE, exact = FALSE)
  
  #save results to dataframe
  results_Wilcoxon_KNN$p_value[i] <- test_result$p.value
  results_Wilcoxon_KNN$statistic[i] <- test_result$statistic
}


#adjust p-value for multiple testing (Benjamini-Hochberg
results_Wilcoxon_KNN$adj_p_value <- p.adjust(results_Wilcoxon_KNN$p_value, method = "BH")

print(results_Wilcoxon_KNN)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_Wilcoxon_KNN <- results_Wilcoxon_KNN %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_Wilcoxon_KNN) #81 out of 84 were significant and with CV 59/60 significant

# Part 3 -------
# Run unpaired t-test for each metabolite 
#dataframe for the results
results_ttest_KNN <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)

#loop through each metabolite and run t-test for each metabolites
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  #use Welchs test here (not same variance, to test H0 where two groups have equal mean)
  test_result2 <- t.test(imputed_KNN23[[metabolite]], imputed_KNN24[[metabolite]], paired = FALSE, var.equal = FALSE)
  
  #save results to dataframe
  results_ttest_KNN$p_value[i] <- test_result2$p.value
  results_ttest_KNN$statistic[i] <- test_result2$statistic
}

#adjust p-value for multiple testing (Benjamini-Hochberg)
results_ttest_KNN$adj_p_value <- p.adjust(results_ttest_KNN$p_value, method = "BH")

print(results_ttest_KNN)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_ttest_KNN <- results_ttest_KNN %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_ttest_KNN) #6 out of 84 were significant and 3/60 with CV

#grouped bar plot 
#count total number of metabolites
total_metabolites <- length(metabolite_cols) #84

#nr of significant metabolites for each test
significant_wilcoxon <- nrow(significant_results_Wilcoxon_KNN) #81
significant_ttest <- nrow(significant_results_ttest_KNN) #6

#create summary dataframe
test_results <- data.frame(
  Test = rep(c("Wilcoxon", "T-test"), each = 2),
  Category = rep(c("Significant", "Non-Significant"), 2),
  Count = c(significant_wilcoxon, total_metabolites - significant_wilcoxon,
            significant_ttest, total_metabolites - significant_ttest)
)

pdf("/Users/marcinebessire/Desktop/project/KNN_Significance_CV30.pdf", width = 10, height = 6)

#plot grouped bar chart
ggplot(test_results, aes(x = Test, y = Count, fill = Category)) +
  geom_bar(position="dodge", stat="identity") +  #grouped bars
  theme_minimal() +
  labs(title = "Number of Significant vs Non-Significant Metabolites using KNN Imputation",
       x = "Statistical Test",
       y = "Count of Metabolites") +
  scale_fill_manual(values = c("Significant" = "blue", "Non-Significant" ="lightblue")) +
  theme(legend.position = "top") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5)  #add labels

dev.off()

# Part 3.2 -------
# Check normality of data 

#check normality using shapiro.test (becuase t-test assumes normality)
shapiro_results23 <- apply(imputed_KNN23, 2, function(x) shapiro.test(x)$p.value)
shapiro_results24 <- apply(imputed_KNN24, 2, function(x) shapiro.test(x)$p.value)


#convert to a dataframe for easy viewing
shapiro_df23 <- data.frame(Metabolite = names(shapiro_results23), p_value = shapiro_results23) #total 84 metabolites
shapiro_df24 <- data.frame(Metabolite = names(shapiro_results24), p_value = shapiro_results24)

#if p-value < 0.05 then not normal distribution
non_normal_count23 <- sum(shapiro_df23$p_value < 0.05)
non_normal_count23 #61 metabolites are non-normal distributed and 45 with CV
non_normal_count24 <- sum(shapiro_df24$p_value < 0.05)
non_normal_count24 #84 metabolites are non-normal distributed and 60 with CV


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
    labs(title = paste("Distribution of Original Data, Imputed Data, and Imputed Values with CV filtering (", year, ")", sep = ""),
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

density_plot_KNN_23 <- plot_imputation_distribution(numeric23, imputed_KNN23, "2023", "/Users/marcinebessire/Desktop/project/KNN_Distribution_Comparison23_CV30.pdf")
density_plot_KNN_24 <- plot_imputation_distribution(numeric24, imputed_KNN24, "2024", "/Users/marcinebessire/Desktop/project/KNN_Distribution_Comparison24_CV30.pdf")

# #QQ-plot to assess if distribution are the smae 
# plot_qq_comparison <- function(original_data, imputed_data, year, output_file) {
#   #convert data to long format
#   imputed_long <- imputed_data %>%
#     pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")
#   
#   original_long <- original_data %>%
#     pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")
#   
#   #merge original and imputed datasets
#   comparison <- original_long %>%
#     left_join(imputed_long, by = "Metabolite")
#   
#   #save plot
#   pdf(output_file, width = 8, height = 6)
#   
#   #generate the QQ plot
#   p <- ggplot(comparison, aes(sample = Imputed_Data)) +
#     stat_qq(aes(sample = Original_Data)) + 
#     stat_qq_line(aes(sample = Original_Data), color = "red") +
#     labs(title = paste("QQ Plot: Original vs Imputed Data (", year, ")", sep = ""),
#          x = "Original Data Quantiles",
#          y = "Imputed Data Quantiles") +
#     theme_minimal()
#   
#   print(p)
#   dev.off()
# }
# 
# # Generate QQ plots for 2023 and 2024
# qq_plot_23 <- plot_qq_comparison(numeric23, imputed_KNN23, "2023", "/Users/marcinebessire/Desktop/project/KNN_QQ_Plot_Comparison23.pdf")
# qq_plot_24 <- plot_qq_comparison(numeric24, imputed_KNN24, "2024", "/Users/marcinebessire/Desktop/project/KNN_QQ_Plot_Comparison24.pdf")


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
    labs(title = paste("Normalized Difference in Mean Before and After Imputation with CV filtering (", year, ")", sep = ""),
         x = "Metabolite",
         y = "Normalized Difference") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
  
  print(p1)
  
  #plot density of the normalized difference
  p2 <- ggplot(mean_comparison, aes(x = Normalized_Difference)) +
    geom_density(fill = "blue", alpha = 0.4, color = "black") +  # Density plot
    theme_minimal() +
    labs(title = paste("Density Plot of Normalized Difference with KNN Imputation and CV filtering (", year, ")", sep = ""),
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.4, 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")  # Reference line at 0
  
  print(p2)
  
  dev.off()
}

normalized_difference23 <- calculate_normalized_difference(numeric23, imputed_KNN23, "2023", "/Users/marcinebessire/Desktop/project/KNN_normalized_difference23_CV30.pdf")
normalized_difference24 <- calculate_normalized_difference(numeric24, imputed_KNN24, "2024", "/Users/marcinebessire/Desktop/project/KNN_normalized_difference24_CV30.pdf")



# # Extra Part (Correlation)
# # Assess quality of KNN imputation 
# 
# #first check if all values were imputed (no NA left)
# sum(is.na(imputed_KNN23)) #is 0
# sum(is.na(imputed_KNN24)) #is 0 
# 
# #compare data distribution before and after imputation 
# par(mfrow = c(3,3)) #3x3 grid
# par(mar = c(4,4,2,1))
# 
# #Boxplots of 9 metabolites (because 84 too big) for 2023
# for (i in 1:9) {
#   boxplot(imputed_KNN23[[i]], numeric23[[i]],
#           names = c("KNN Imputation", "No Imputation"),
#           main = paste("Metabolite", colnames(imputed_KNN23)[i]),
#           col = c("lightblue", "lightgreen"))
# }
# 
# #Boxplots of 9 metabolites (because 84 too big) for 2024
# for (i in 1:9) {
#   boxplot(imputed_KNN24[[i]], numeric24[[i]],
#           names = c("KNN Imputation", "No Imputation"),
#           main = paste("Metabolite", colnames(imputed_KNN23)[i]),
#           col = c("lightblue", "lightgreen"))
# }


# #Correlation Check, compare correlation matrices
# cor_before_23 <- cor(numeric23, use = "pairwise.complete.obs")
# cor_after_23 <- cor(imputed_KNN23)
# 
# cor_before_24 <- cor(numeric24, use = "pairwise.complete.obs")
# cor_after_24 <- cor(imputed_KNN24)
# 
# # Scatter plot 
# #identify missing vlaues 
# missing_val23 <- is.na(numeric23) #99
# missing_val24 <- is.na(numeric24) #526
# 
# 
# #create logical matrix identifying if a pair of metabolites had missing values
# missing_pairs23 <- (missing_val23 %*% t(missing_val23)) > 0 #True 
# missing_pairs24 <- (missing_val24 %*% t(missing_val24)) > 0
# 
# #extract upper triangular part of matrix 
# #2023
# imputed_indices23 <- missing_pairs23[upper.tri(missing_pairs23, diag = FALSE)]
# cor_values_before23 <- cor_before_23[upper.tri(cor_before_23, diag = FALSE)]
# cor_values_after23 <- cor_after_23[upper.tri(cor_after_23, diag = FALSE)]
# #2024
# imputed_indices24 <- missing_pairs24[upper.tri(missing_pairs24, diag = FALSE)]
# cor_values_before24 <- cor_before_24[upper.tri(cor_before_24, diag = FALSE)]
# cor_values_after24 <- cor_after_24[upper.tri(cor_after_24, diag = FALSE)]
# 
# pdf("/Users/marcinebessire/Desktop/project/KNN_Correlation_Comparison.pdf", width = 8, height = 6)
# 
# #2023
# plot(cor_values_before23, cor_values_after23,
#      xlab = "Before Imputation",
#      ylab = "After Imputation",
#      main = "Correlation Comparison",
#      col = ifelse(imputed_indices23, "red", "blue"), #red for imputed pairs
#      pch = 19) 
# 
# abline(0,1,col = "black", lwd = 2) #ideally points should line on line
# 
# #2024
# plot(cor_values_before24, cor_values_after24,
#      xlab = "Before Imputation",
#      ylab = "After Imputation",
#      main = "Correlation Comparison",
#      col = ifelse(imputed_indices24, "red", "blue"), #red for imputed pairs
#      pch = 19) 
# 
# abline(0,1,col = "black", lwd = 2) #ideally points should line on line
# 
# dev.off()

