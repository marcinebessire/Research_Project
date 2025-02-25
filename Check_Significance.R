# Load required library
library(tidyverse)

# Signficance Test for 2023
#load data (after Half-Min imputation) 
data23 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Half_Min_Imputation23.csv", check.names = FALSE)
data24 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Half_Min_Imputation24.csv", check.names = FALSE)

#make function to make t-test
t_test_func <- function(data) {
  #take every second row and make 2 separate dataframes
  df_1 <- data[seq(1, nrow(data), 2), ]
  df_2 <- anti_join(data, df_1)
  
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
  
  return(results_ttest)
  
}

#perform t-test
t_test_23 <- t_test_func(data23)
t_test_24 <- t_test_func(data24)

#view significant metabolites with BH adjusted p-value < 0.05
check_significance <- function(data) {
  significant_results_ttest <- data %>% 
    filter(adj_p_value < 0.05) %>% #usually 0.05 used 
    arrange(adj_p_value)
  
  return(significant_results_ttest)
}

#check significance 
res_23 <- check_significance(t_test_23)
res_24 <- check_significance(t_test_24) #also 0 

# Significance test for Before and After New Chip from Data 2024
#cut data of 2024 into two parts (before and after new chip)
before_chip_data24 <- data24 %>% #before new chip
  filter(MonthDay < "09-10") 
after_chip_data24 <- anti_join(data24,before_chip_data24) #after new chip 

func_before_after <- function(data1, data2) {
  #select only metabolite columns (starting from column 6)
  df_1 <- data1[, 6:ncol(data1)]
  df_2 <- data2[, 6:ncol(data2)]
  
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
  
  return(results_ttest)
  
}

#perform t_test
t_test_before_after24 <- func_before_after(before_chip_data24, after_chip_data24)

#check significance 
res_before_after24 <- check_significance(t_test_before_after24) #34 
