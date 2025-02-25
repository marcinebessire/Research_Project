# Load required library
library(tidyverse)
library(ggplot2)
library(corrplot)

# Part 1 ------
# Perform Half-min imputation

#load data 
final_data_2023 <- read.csv("/Users/marcinebessire/Desktop/project/Final_Data_2023.csv", check.names = FALSE)
final_data_2024 <- read.csv("/Users/marcinebessire/Desktop/project/Final_Data_2024.csv", check.names = FALSE)

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
cat("Number of common columns:", length(common_cols), "\n") #84
cat("Common columns:", paste(common_cols, collapse=", ")) 

#make new dataframe for each year with those 84 columns + metadata
final_col <- c(exclude_cols, common_cols)
half_min_23 <- imputed_data_23 %>% select(all_of(final_col))
half_min_24 <- imputed_data_24 %>% select(all_of(final_col))

#save common imputed metabolites to new file
write_csv(half_min_23, "/Users/marcinebessire/Desktop/project/Common_Half_Min_Imputation23.csv")
write_csv(half_min_24, "/Users/marcinebessire/Desktop/project/Common_Half_Min_Imputation24.csv")


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

print(significant_results_Wilcoxon)

# Part 3.2 --------
# Wilcoxon rank-sum test (for cut data)

#select only metabolite columns for cut data
#cut_imputed_data_23 <- cut_imputed_data_23[, 6:ncol(cut_imputed_data_23)]
#cut_imputed_data_24 <- cut_imputed_data_24[, 6:ncol(cut_imputed_data_24)]

# cut_results_Wilcoxon <- data.frame(
#   Metabolite = cut_metabolite_cols,
#   p_value = numeric(length(cut_metabolite_cols)),
#   statistic = numeric(length(cut_metabolite_cols))
# )

#loop through each metabolite and run Wilcoxon rank-sum test for cut dataset
# for (i in seq_along(cut_metabolite_cols)) {
#   cut_metabolite <- cut_metabolite_cols[i]
#   
#   cut_test_result <- wilcox.test(cut_imputed_data_23[[cut_metabolite]], cut_imputed_data_24[[cut_metabolite]], paired = FALSE, exact = FALSE)
#   
#   #save results to dataframe
#   cut_results_Wilcoxon$p_value[i] <- cut_test_result$p.value
#   cut_results_Wilcoxon$statistic[i] <- cut_test_result$statistic
# }

#for cut data 
#cut_results_Wilcoxon$adj_p_value <- p.adjust(cut_results_Wilcoxon$p_value, method = "BH")

#print(cut_results_Wilcoxon)

#view significant metabolites with BH adjusted p-value < 0.05 for cut data
# cut_significant_results_Wilcoxon <- cut_results_Wilcoxon %>% 
#   filter(adj_p_value < 0.05) %>% #usually 0.05 used 
#   arrange(adj_p_value)
# 
# print(cut_significant_results_Wilcoxon)

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

print(significant_results_ttest)

# Part 6 -------
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

#measure Spearman rank correlation (monotinic relationship)  
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


