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
data23 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites23.csv", check.names = FALSE)
data24 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites24.csv", check.names = FALSE)

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
# Distribution plot before and after Imputation and Imputed Values Only 

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
dist_knn23 <- plot_whole_distribution(numeric23, imputed_KNN23, "KNN", "2023")
dist_knn24 <- plot_whole_distribution(numeric24, imputed_KNN24, "KNN", "2024")


# Part 8 ------
# calculate normalized difference of each imputation (before and after) and plot

#mean before imputation
mean_before23 <- numeric23 %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_Before")

sd_before23 <- numeric23 %>%
  summarise(across(everything(), sd, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "SD_Before")


#means after imputation
mean_after23 <- imputed_KNN23 %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_After")

#merge before and after mean values
mean_comparison23 <- left_join(mean_before23, mean_after23, by = "Metabolite")

#compute normalized difference: (Mean_After - Mean_Before) / Mean_Before
mean_comparison23 <- mean_comparison23 %>%
  mutate(Normalized_Difference = (Mean_After - Mean_Before) / Mean_Before)

#mean before imputation
mean_before24 <- numeric24 %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_Before")

#means after imputation
mean_after24 <- imputed_KNN24 %>%
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
  labs(title = "Normalized Difference in Mean Before and After with KNN Imputation (2023)",
       x = "Metabolite",
       y = "Normalized Difference") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)


#plot normalized difference
ggplot(mean_comparison24, aes(x = Metabolite, y = Normalized_Difference, fill = Normalized_Difference)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Normalized Difference in Mean Before and After with KNN Imputation (2024)",
       x = "Metabolite",
       y = "Normalized Difference") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

#plot the density of the normalized difference
ggplot(mean_comparison23, aes(x = Normalized_Difference)) +
  geom_density(fill = "blue", alpha = 0.4, color = "black") +  #density plot
  theme_minimal() +
  labs(title = "Density Plot of Normalized Difference with KNN Imputation (2023)",
       x = "Normalized Difference",
       y = "Density") +
  xlim(-0.2,0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")  #reference line at 0

#plot the density of the normalized difference
ggplot(mean_comparison24, aes(x = Normalized_Difference)) +
  geom_density(fill = "blue", alpha = 0.4, color = "black") +  #density plot
  theme_minimal() +
  labs(title = "Density Plot of Normalized Difference with KNN Imputation (2024)",
       x = "Normalized Difference",
       y = "Density") +
  xlim(-0.4,0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")  #reference line at 0


# Part 9 -----
# Kolomogorov-Smirnov test: nonparametric test to test whether two ssamples came from same distirbution

#convert data to long format
#2023
imputed_23_long <- imputed_KNN23 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")

original_23_long <- numeric23 %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")

#2024
#convert data to long format
imputed_24_long <- imputed_KNN24 %>%
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




