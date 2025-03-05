#Load required library
library(tidyverse)
library(ggplot2)
library(corrplot)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)

# Part 1 ------
# Perform Half-min imputation

#load data with CV threshold 30
final_data_2023 <- read.csv("/Users/marcinebessire/Desktop/project/Final_Data_2023.csv", check.names = FALSE)
final_data_2024 <- read.csv("/Users/marcinebessire/Desktop/project/Final_Data_2024.csv", check.names = FALSE)

#load original common metabolites
original_23 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites23.csv", check.names = FALSE)
original_24 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites24.csv", check.names = FALSE)

#keep only numeric columns
original_23_metabolites <- original_23[, 6:ncol(original_23)]
original_24_metabolites <- original_24[, 6:ncol(original_24)]

#count how many missing values are there
mv23 <- sum(is.na(original_23_metabolites)) #26 MV out of 4980 => 0.522%
mv24 <- sum(is.na(original_24_metabolites)) #333 MV out of 4320 => 7.708%

#imputation for cut 2024
#cut_2024 <- read.csv("/Users/marcinebessire/Desktop/project/Cut_Common_Metabolites24.csv", check.names = FALSE)
#cut_2023 <- read.csv("/Users/marcinebessire/Desktop/project/Cut_Common_Metabolites23.csv", check.names = FALSE)

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

#cut data
#cut_imputed_data_23 <- half_min_imputation(cut_2023, "2023")
#cut_imputed_data_24 <- half_min_imputation(cut_2024, "2024")

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
cat("Number of common columns:", length(common_cols), "\n") #84 and 60 if CV 30%
cat("Common columns:", paste(common_cols, collapse=", ")) 

#make new dataframe for each year with those 84 columns + metadata
final_col <- c(exclude_cols, common_cols)
half_min_23 <- imputed_data_23 %>% select(all_of(final_col))
half_min_24 <- imputed_data_24 %>% select(all_of(final_col))

#save common imputed metabolites to new file
#write_csv(half_min_23, "/Users/marcinebessire/Desktop/project/Common_Half_Min_Imputation23.csv")
#write_csv(half_min_24, "/Users/marcinebessire/Desktop/project/Common_Half_Min_Imputation24.csv")


# Part 3 --------
# Wilcoxon rank-sum test (for independent data)

#select only metabolite columns (starting from column 6)
half_min_23_metabolites <- half_min_23[, 6:ncol(half_min_23)]
half_min_24_metabolites <- half_min_24[, 6:ncol(half_min_24)]

#get metabolite column names
metabolite_cols <- colnames(half_min_23_metabolites)
#cut_metabolite_cols <- colnames(cut_imputed_data_23)

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

print(significant_results_Wilcoxon) #78 out of 84 and if CV 30% then 56/60


# Part 4 -------
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

print(significant_results_ttest) #74 out of 84 and with CV 30% 53/60

#grouped bar plot 
#count total number of metabolites
total_metabolites <- length(metabolite_cols) #84

#nr of significant metabolites for each test
significant_wilcoxon <- nrow(significant_results_Wilcoxon) #78
significant_ttest <- nrow(significant_results_ttest) #74

#create summary dataframe
test_results <- data.frame(
  Test = rep(c("Wilcoxon", "T-test"), each = 2),
  Category = rep(c("Significant", "Non-Significant"), 2),
  Count = c(significant_wilcoxon, total_metabolites - significant_wilcoxon,
            significant_ttest, total_metabolites - significant_ttest)
)

pdf("/Users/marcinebessire/Desktop/project/Halfmin_Significance.pdf", width = 10, height = 6)

#plot grouped bar chart
ggplot(test_results, aes(x = Test, y = Count, fill = Category)) +
  geom_bar(position="dodge", stat="identity") +  #grouped bars
  theme_minimal() +
  labs(title = "Number of Significant vs Non-Significant Metabolites using Half-min Imputation",
       x = "Statistical Test",
       y = "Count of Metabolites") +
  scale_fill_manual(values = c("Significant" = "blue", "Non-Significant" ="lightblue")) +
  theme(legend.position = "top") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5)  #add labels

dev.off()

# Part 4.2 -----
# Check normality of data 

#check normality using shapiro.test (becuase t-test assumes normality)
shapiro_results23 <- apply(half_min_23_metabolites, 2, function(x) shapiro.test(x)$p.value)
shapiro_results24 <- apply(half_min_24_metabolites, 2, function(x) shapiro.test(x)$p.value)
#check before imputation
shapiro_results_original23 <- apply(original_23_metabolites, 2, function(x) shapiro.test(x)$p.value)

#convert to a dataframe for easy viewing
shapiro_df23 <- data.frame(Metabolite = names(shapiro_results23), p_value = shapiro_results23) #total 84 metabolites
shapiro_df24 <- data.frame(Metabolite = names(shapiro_results24), p_value = shapiro_results24)
shapiro_df_original23 <- data.frame(Metabolite = names(shapiro_results_original23), p_value = shapiro_results_original23)

#if p-value < 0.05 then not normal distribution
non_normal_count23 <- sum(shapiro_df23$p_value < 0.05)
non_normal_count23 #64 metabolites are non-normal distributed and with CV 47
non_normal_count24 <- sum(shapiro_df24$p_value < 0.05)
non_normal_count24 #84 metabolites are non-normal distributed and with CV 60
non_normal_count_original23 <- sum(shapiro_df_original23$p_value < 0.05)
non_normal_count_original23 #59 metabolites are non-normal distributed and with CV 60


# Part 5 -------
# Correlation Coefficient for each year (between CV (values before imputation) and Metabolite after Half min Imputation)

#load CV data for both years 
cv_res23 <- read.csv("/Users/marcinebessire/Desktop/project/Common_CV_results23.csv", check.names = FALSE)
cv_res24 <- read.csv("/Users/marcinebessire/Desktop/project/Common_CV_results24.csv", check.names = FALSE)

#t-test results (p-values)
df_ttest <- results_ttest

#rename columns for clarity 
colnames(cv_res23) <- c("Metabolite", "CV [%]")
colnames(cv_res24) <- c("Metabolite", "CV [%]")

#merge dataframe by metabolites
merged_23 <- merge(cv_res23, df_ttest, by = "Metabolite")
merged_24 <- merge(cv_res24, df_ttest, by = "Metabolite")

#calculate correlation coefficient between CV and adjusted p-value 
#2023
cor_23 <- cor(merged_23$`CV [%]`, merged_23$adj_p_value, method = "pearson")
print(paste("Correlation coefficient for 2023:", cor_23)) #-0.0293
#2024
cor_24 <- cor(merged_24$`CV [%]`, merged_24$adj_p_value, method = "pearson")
print(paste("Correlation coefficient for 2024:", cor_24)) #-0.0837

#now plot the data of 2023
ggplot(merged_23, aes(x = `CV [%]`, y = adj_p_value)) +
  geom_point() +
  labs(title = "2023: CV vs Adjusted P-Value",
       x = "CV [%]", 
       y = "Adjusted P-Value") +
  theme_minimal()

#now plot the data of 2024
ggplot(merged_24, aes(x = `CV [%]`, y = adj_p_value)) +
  geom_point() +
  labs(title = "2023: CV vs Adjusted P-Value",
       x = "CV [%]", 
       y = "Adjusted P-Value") +
  xlim(0.0, 50) +
  ylim(0.0, 0.01)
  theme_minimal()

#measure Spearman rank correlation (monotonic relationship)  
cor_23_spearman <- cor(merged_23$CV, merged_23$adj_p_value, method = "spearman")
cor_24_spearman <- cor(merged_24$CV, merged_24$adj_p_value, method = "spearman")
print(paste("Spearman Correlation for 2023:", cor_23_spearman)) #-0.105
print(paste("Spearman Correlation for 2024:", cor_24_spearman)) #-0.05

#check if outliers impact the relationship
threshold_23 <- quantile(merged_23$`CV [%]`, 0.95)  #95th percentile fitlers top 5% of CV values
threshold_24 <- quantile(merged_24$`CV [%]`, 0.95)  

filtered_23 <- merged_23 %>% filter(`CV [%]` < threshold_23)
filtered_24 <- merged_24 %>% filter(`CV [%]` < threshold_24)

cor_23_filtered <- cor(filtered_23$`CV [%]`, filtered_23$adj_p_value, method = "pearson")
cor_24_filtered <- cor(filtered_24$`CV [%]`, filtered_24$adj_p_value, method = "pearson")

print(paste("Filtered Correlation for 2023:", cor_23_filtered)) #-0.2676
print(paste("Filtered Correlation for 2024:", cor_24_filtered)) #-0.2495

dev.off()

# # Part 6 -------
# # check correlation coefficient before and after imputation for each year
# 
# #correlation Check, compare correlation matrices
# cor_before_23 <- cor(original_23_metabolites, use = "pairwise.complete.obs")
# cor_after_23 <- cor(half_min_23_metabolites)
# 
# cor_before_24 <- cor(original_24_metabolites, use = "pairwise.complete.obs")
# cor_after_24 <- cor(half_min_24_metabolites)
# 
# # Scatter plot 
# #identify missing vlaues 
# missing_val23 <- is.na(original_23_metabolites) #99
# missing_val24 <- is.na(original_24_metabolites) #526
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
# 
# #2024
# imputed_indices24 <- missing_pairs24[upper.tri(missing_pairs24, diag = FALSE)]
# cor_values_before24 <- cor_before_24[upper.tri(cor_before_24, diag = FALSE)]
# cor_values_after24 <- cor_after_24[upper.tri(cor_after_24, diag = FALSE)]
# 
# pdf("/Users/marcinebessire/Desktop/project/Halfmin_Correlation_Comparison.pdf", width = 8, height = 6)
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


# Part 7 -----
# Distribution plot of data before and after Imputation 

#2023
#convert data to long format for visualization
half_min_23_long <- half_min_23_metabolites %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")

original_23_long <- original_23_metabolites %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")

#identify imputed values => missing values were replaced by half of the minimum observed value.
imputed_only_23 <- original_23_long %>%
  mutate(Imputed = is.na(Original_Data)) %>%
  filter(Imputed) %>%
  inner_join(half_min_23_long %>%
               distinct(Metabolite, .keep_all = TRUE), by = "Metabolite") %>%
  mutate(Dataset = "Imputed_Values")


#merge original and imputed datasets
comparison_23 <- original_23_long %>%
  left_join(half_min_23_long %>%
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
half_min_24_long <- half_min_24_metabolites %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")

original_24_long <- original_24_metabolites %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")

#identify imputed values => missing values were replaced by half of the minimum observed value.
imputed_only_24 <- original_24_long %>%
  mutate(Imputed = is.na(Original_Data)) %>%
  filter(Imputed) %>%
  inner_join(half_min_24_long %>%
               distinct(Metabolite, .keep_all = TRUE), by = "Metabolite") %>%
  mutate(Dataset = "Imputed_Values")

#merge original and imputed datasets
comparison_24 <- original_24_long %>%
  left_join(half_min_24_long %>%
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
pdf("/Users/marcinebessire/Desktop/project/Halfmin_Distribution_Comparison.pdf", width = 8, height = 6)

#2023 plot 
ggplot(comparison_23, aes(x = Value, fill = Dataset)) +
  geom_density(alpha = 0.5) +  # Transparency for overlapping
  labs(title = "Distribution of Original/Imputed Data and Imputed Values with Half-min Imputation (2023)",
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
  labs(title = "Distribution of Original/Imputed Data, and Imputed Values with Half-min Imputation (2024)",
       x = "Metabolite Value",
       y = "Density") +
  theme_minimal() + 
  scale_fill_manual(values = c("Original_Data" = "lightblue", 
                               "Imputed_Data" = "red", 
                               "Imputed_Values" = "green")) +
  xlim(-10,50)
  
dev.off()

# Part 7.2 -----
# QQ plot comparing distribution

pdf("/Users/marcinebessire/Desktop/project/Halfmin_QQ_Plots.pdf", width = 8, height = 6)

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
pdf("/Users/marcinebessire/Desktop/project/Halfmin_QQ_Plots_NoOutliers.pdf", width = 8, height = 6)

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

# Part 8 ------
# calculate normalized difference of each imputation (before and after) and plot

#mean before imputation
mean_before23 <- original_23_metabolites %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_Before")

sd_before23 <- original_23_metabolites %>%
  summarise(across(everything(), sd, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "SD_Before")


#means after imputation
mean_after23 <- half_min_23_metabolites %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_After")

#merge before and after mean values
mean_comparison23 <- left_join(mean_before23, mean_after23, by = "Metabolite")

#compute normalized difference: (Mean_After - Mean_Before) / Mean_Before
mean_comparison23 <- mean_comparison23 %>%
  mutate(Normalized_Difference = (Mean_After - Mean_Before) / Mean_Before)

#mean before imputation
mean_before24 <- original_24_metabolites %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_Before")

#means after imputation
mean_after24 <- half_min_24_metabolites %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_After")

#merge before and after mean values
mean_comparison24 <- left_join(mean_before24, mean_after24, by = "Metabolite")

#compute normalized difference: (Mean_After - Mean_Before) / Mean_Before
mean_comparison24 <- mean_comparison24 %>%
  mutate(Normalized_Difference = (Mean_After - Mean_Before) / Mean_Before)

pdf("/Users/marcinebessire/Desktop/project/Halfmin_Normalized_Difference_Comparison.pdf", width = 10, height = 6)

#plot normalized difference
ggplot(mean_comparison23, aes(x = Metabolite, y = Normalized_Difference, fill = Normalized_Difference)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Normalized Difference in Mean Before and After Half-min Imputation (2023)",
       x = "Metabolite",
       y = "Normalized Difference") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)


#plot normalized difference
ggplot(mean_comparison24, aes(x = Metabolite, y = Normalized_Difference, fill = Normalized_Difference)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Normalized Difference in Mean Before and After Half-min Imputation (2024)",
       x = "Metabolite",
       y = "Normalized Difference") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

#plot the density of the normalized difference
ggplot(mean_comparison23, aes(x = Normalized_Difference)) +
  geom_density(fill = "blue", alpha = 0.4, color = "black") +  #density plot
  theme_minimal() +
  labs(title = "Density Plot of Normalized Difference with Half-min Imputation (2023)",
       x = "Normalized Difference",
       y = "Density") +
  xlim(-0.2,0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")  #reference line at 0

#plot the density of the normalized difference
ggplot(mean_comparison24, aes(x = Normalized_Difference)) +
  geom_density(fill = "blue", alpha = 0.4, color = "black") +  #density plot
  theme_minimal() +
  labs(title = "Density Plot of Normalized Difference with Half-min Imputation (2024)",
       x = "Normalized Difference",
       y = "Density") +
  xlim(-0.4,0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")  #reference line at 0

dev.off()

# FINALLY: Check row count again 
#2023: should have 99 MV and 6972 total values
dataset_counts <- comparison_23 %>%
  group_by(Dataset) %>%
  summarise(Count = n())

print(dataset_counts)
#2024: should have 526 MV and 6048 total values
dataset_counts <- comparison_24 %>%
  group_by(Dataset) %>%
  summarise(Count = n())

print(dataset_counts)



