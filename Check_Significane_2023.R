# Load required library
library(tidyverse)

#load data (after Half-Min imputation) 
data23 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Half_Min_Imputation23.csv", check.names = FALSE)

#take every second row and make 2 separate dataframes
df_1 <- data23[seq(1, nrow(data23), 2), ]
df_2 <- anti_join(data23, df_1)

#select only metabolite columns (starting from column 6)
df_1 <- df_1[, 6:ncol(df_1)]
df_2 <- df_2[, 6:ncol(df_2)]

#get metabolite column names
metabolite_cols <- colnames(df_1)

#Run unpaired t-test for each metabolite 
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
  test_result <- t.test(df_1[[metabolite]], df_2[[metabolite]], paired = FALSE, var.equal = FALSE)
  
  #save results to dataframe
  results_ttest$p_value[i] <- test_result$p.value
  results_ttest$statistic[i] <- test_result$statistic
}

#adjust p-value for multiple testing (Benjamini-Hochberg)
results_ttest$adj_p_value <- p.adjust(results_ttest$p_value, method = "BH")

print(results_ttest)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_ttest <- results_ttest %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_ttest)
