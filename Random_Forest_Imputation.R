#load necessary library
library(missForest)
library(dplyr)
library(corrplot)
library(tidyverse)

# Part 1 -------
# Random Forest impputation of MV

#load common metabolite dataframe (2023 and 2024) 
data_23 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites23.csv", check.names = FALSE)
data_24 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites24.csv", check.names = FALSE)

#Only use nuemric columns for imputation
numeric_23 <- data_23[, 6:ncol(data_23)]
numeric_24 <- data_24[, 6:ncol(data_24)]

#make function for RF
RF_imputation <- function(data, ...) {
  #apply miss forest on data to impute MV
  imputed_data <- missForest(data, maxiter = 10, ntree = 100)
  
  #return imputed data
  return(imputed_data)
  
}

#Apply RF function on data 
RF_data23 <- RF_imputation(numeric_23)
RF_data24 <- RF_imputation(numeric_24)

#extract imputed data and save as df
imputed_RF_23 <- as.data.frame(RF_data23$ximp)
imputed_RF_24 <- as.data.frame(RF_data24$ximp)

# Part 2 --------
# Wilcoxon rank-sum test (for independent data)

#get metabolite column names
metabolite_cols <- colnames(imputed_RF_23)

#dataframe for the results
results_Wilcoxon_RF <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)


#loop through each metabolite and run Wilcoxon rank-sum test
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  test_result <- wilcox.test(imputed_RF_23[[metabolite]], imputed_RF_24[[metabolite]], paired = FALSE, exact = FALSE)
  
  #save results to dataframe
  results_Wilcoxon_RF$p_value[i] <- test_result$p.value
  results_Wilcoxon_RF$statistic[i] <- test_result$statistic
}


#adjust p-value for multiple testing (Benjamini-Hochberg
results_Wilcoxon_RF$adj_p_value <- p.adjust(results_Wilcoxon_RF$p_value, method = "BH")

print(results_Wilcoxon_RF)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_Wilcoxon_RF <- results_Wilcoxon_RF %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_Wilcoxon_RF) #82 out of 84 were significant 

# Part 3 -------
# Run unpaired t-test for each metabolite 
#dataframe for the results
results_ttest_RF <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)

#loop through each metabolite and run t-test for each metabolites
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  #use Welchs test here (not same variance, to test H0 where two groups have equal mean)
  test_result2 <- t.test(imputed_RF_23[[metabolite]], imputed_RF_24[[metabolite]], paired = FALSE, var.equal = FALSE)
  
  #save results to dataframe
  results_ttest_RF$p_value[i] <- test_result2$p.value
  results_ttest_RF$statistic[i] <- test_result2$statistic
}

#adjust p-value for multiple testing (Benjamini-Hochberg)
results_ttest_RF$adj_p_value <- p.adjust(results_ttest_RF$p_value, method = "BH")

print(results_ttest_RF)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_ttest_RF <- results_ttest_RF %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_ttest_RF) #80 out of 84 were significant


#grouped bar plot 
#count total number of metabolites
total_metabolites <- length(colnames(imputed_RF_23)) #84

#nr of significant metabolites for each test
significant_wilcoxon <- nrow(significant_results_Wilcoxon_RF) #82
significant_ttest <- nrow(significant_results_ttest_RF) #80


#create summary dataframe
test_results <- data.frame(
  Test = rep(c("Wilcoxon", "T-test"), each = 2),
  Category = rep(c("Significant", "Non-Significant"), 2),
  Count = c(significant_wilcoxon, total_metabolites - significant_wilcoxon,
            significant_ttest, total_metabolites - significant_ttest)
)

pdf("/Users/marcinebessire/Desktop/project/RF_Significance.pdf", width = 10, height = 6)

#plot grouped bar chart
ggplot(test_results, aes(x = Test, y = Count, fill = Category)) +
  geom_bar(position="dodge", stat="identity") +  #grouped bars
  theme_minimal() +
  labs(title = "Number of Significant vs Non-Significant Metabolites using Random Forest Imputation",
       x = "Statistical Test",
       y = "Count of Metabolites") +
  scale_fill_manual(values = c("Significant" = "blue", "Non-Significant" ="lightblue")) +
  theme(legend.position = "top") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5)  #add labels

dev.off()

# Part 3.2 -------
# Check normality of data 

#check normality using shapiro.test (becuase t-test assumes normality)
shapiro_results23 <- apply(imputed_RF_23, 2, function(x) shapiro.test(x)$p.value)
shapiro_results24 <- apply(imputed_RF_24, 2, function(x) shapiro.test(x)$p.value)


#convert to a dataframe for easy viewing
shapiro_df23 <- data.frame(Metabolite = names(shapiro_results23), p_value = shapiro_results23) #total 84 metabolites
shapiro_df24 <- data.frame(Metabolite = names(shapiro_results24), p_value = shapiro_results24)

#if p-value < 0.05 then not normal distribution
non_normal_count23 <- sum(shapiro_df23$p_value < 0.05)
non_normal_count23 #61 metabolites are non-normal distributed 
non_normal_count24 <- sum(shapiro_df24$p_value < 0.05)
non_normal_count24 #64 metabolites are non-normal distributed 


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

density_plot_QRILC_23 <- plot_imputation_distribution(numeric_23, imputed_RF_23, "2023", "/Users/marcinebessire/Desktop/project/RF_Distribution_Comparison23.pdf")
density_plot_QRILC_24 <- plot_imputation_distribution(numeric_24, imputed_RF_24, "2024", "/Users/marcinebessire/Desktop/project/RF_Distribution_Comparison24.pdf")


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
    labs(title = paste("Density Plot of Normalized Difference with RF Imputation (", year, ")", sep = ""),
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.4, 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")  # Reference line at 0
  
  print(p2)
  
  dev.off()
}

normalized_difference23 <- calculate_normalized_difference(numeric_23, imputed_RF_23, "2023", "/Users/marcinebessire/Desktop/project/RF_normalized_difference23.pdf")
normalized_difference24 <- calculate_normalized_difference(numeric_24, imputed_RF_24, "2024", "/Users/marcinebessire/Desktop/project/RF_normalized_difference24.pdf")


# # Extra Part (Correlation) ----- 
# # Check if Impuation method was good 
# 
# #OOB (out of bag imputation error) is the normalized mean squared error (NMSE) => the lower the better
# #NMSE = 0 means perfect prediction
# OOB23 <- RF_data23$OOBerror #0.005697
# OOB24 <- RF_data24$OOBerror #0.40187 higher error, less effecitve imputation 
# 
# # Correlation analysis 
# cor_before_23 <- cor(numeric_23, use = "pairwise.complete.obs")
# cor_after_23 <- cor(imputed_RF_23)
# 
# # Heat map
# par(mfrow = c(1,2)) #for side by side plots
# corrplot(cor_before_23, method = "color", tl.cex = 0.6, title = "Before Imputation (2023)")
# corrplot(cor_after_23, method = "color", tl.cex = 0.6, title = "After Imputation (2023)")
# dev.off()
# 
# # Scatter plot 
# #identify missing vlaues 
# missing_val <- is.na(numeric_23) #total of 99 are missing
# 
# #create logical matrix identifying if a pair of metabolites had missing values
# missing_pairs <- (missing_val %*% t(missing_val)) > 0 #True 
# 
# #extract upper triangular part of matrix 
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

