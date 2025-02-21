# Load required library
library(tidyverse)

# Part 1 ------
# Perform Half-min imputation

#load data 
data_2023 <- read_csv("/Users/marcinebessire/Desktop/project/Final_Data_2023.csv")
data_2024 <- read_csv("/Users/marcinebessire/Desktop/project/Final_Data_2024.csv")

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
imputed_data_23 <- half_min_imputation(data_2023, "2023")
imputed_data_24 <- half_min_imputation(data_2024, "2024")

# Part 2 ----
# Check Half-Min imputation 
# by comparing the  mean median and standard deviation for and after imputation

#identify columns with missing values for 2023 and 2024
cols_with_na_23 <- colnames(data_2023)[colSums(is.na(data_2023)) > 0]
cols_with_na_24 <- colnames(data_2024)[colSums(is.na(data_2024)) > 0]

#function to compute summary statistics
compute_stats <- function(original_df, imputed_df, imputed_columns){
  #original data
  stats_original <- original_df %>%
    select(all_of(imputed_columns)) %>%
    summarise(across(everything(), list(
      Median = ~median(.x, na.rm = TRUE),
      Mean = ~mean(.x, na.rm = TRUE),
      SD = ~sd(.x, na.rm = TRUE)
    ))) %>%
    mutate(Dataset = "Original") %>%
    relocate(Dataset)
  
  #imputed data
  stats_imputed <- imputed_df %>%
    select(all_of(imputed_columns)) %>%
    summarise(across(everything(), list(
      Median = ~median(.x, na.rm = TRUE),
      Mean = ~mean(.x, na.rm = TRUE),
      SD = ~sd(.x, na.rm = TRUE)
    ))) %>%
    mutate(Dataset = "Imputed") %>%
    relocate(Dataset)
  
  #calculate % change ((imputed-original)/orignal) * 100
  stats_change <- stats_original
  stats_change[,-1] <- ((stats_imputed[,-1] - stats_original[,-1])/stats_original[,-1]) * 100
  stats_change$Dataset <- "% Change"
  
  final_stats <- bind_rows(stats_original, stats_imputed, stats_change)
  
  return(final_stats)
  
}

#compute statistics
summary_comparison23 <- compute_stats(data_2023, imputed_data_23, cols_with_na_23)
summary_comparison24 <- compute_stats(data_2024, imputed_data_24, cols_with_na_24)

#save csv 
output_file23 <- "/Users/marcinebessire/Desktop/project/Summary_Statistics_HalfMin_23.csv"
write_csv(summary_comparison23, output_file23)
output_file24 <- "/Users/marcinebessire/Desktop/project/Summary_Statistics_HalfMin_24.csv"
write_csv(summary_comparison24, output_file24)


