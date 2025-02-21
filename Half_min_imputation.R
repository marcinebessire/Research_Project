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
compute_stats <- function(df, dataset_name,imputed_columns){
  df %>%
    select(all_of(imputed_columns)) %>%
    summarise(
      across(everything(), list(
        Median = ~median(.x, na.rm = TRUE),
        Mean = ~mean(.x, na.rm = TRUE),
        SD = ~sd(.x, na.rm = TRUE)
      ))
    ) %>%
    pivot_longer(cols = everything(), names_to = "Statistic", values_to = "Value") %>%
    mutate(Dataset = dataset_name) %>%
    pivot_wider(names_from = Statistic, values_from = Value) %>%
    relocate(Dataset)
}

#compute statistics
summary_original23 <- compute_stats(data_2023, "Original",cols_with_na_23)
summary_original24 <- compute_stats(data_2024, "Original",cols_with_na_24)
summary_imputed23 <- compute_stats(imputed_data_23, "Imputed",cols_with_na_23)
summary_imputed24 <- compute_stats(imputed_data_24, "Imputed",cols_with_na_24)

#combine results
summary_comparison23 <- bind_rows(summary_original23, summary_imputed23)
summary_comparison24 <- bind_rows(summary_original24, summary_imputed24)

#save csv 
output_file23 <- "/Users/marcinebessire/Desktop/project/Summary_Statistics_HalfMin_23.csv"
write_csv(summary_comparison23, output_file23)
output_file24 <- "/Users/marcinebessire/Desktop/project/Summary_Statistics_HalfMin_24.csv"
write_csv(summary_comparison23, output_file24)


