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

#find imputed column to later select for those 
find_imputed_col <- function(orignal_data, imputed_data){
  imputed_col <- names(orignal_data)[sapply(original_data, function(col) sum(is.na(col)) > 0)]
  return(imputed_col)
}

#function to compare statistics before and after imputation 
compare_stats <- function(original_data, imputed_data){
  stats_original <- original_data %>%
    summarise(across(where(is.numeric), list(mean = ~mean(., na.rm = TRUE),
                                             median = ~median(., na.rm = TRUE),
                                             sd = ~sd(., na.rm = TRUE))))
  
  stats_imputed <- imputed_data %>%
    summarise(across(where(is.numeric), list(mean = ~mean(., na.rm = TRUE),
                                             median = ~median(., na.rm = TRUE),
                                             sd = ~sd(., na.rm = TRUE))))
  
  stats_change <- stats_original - stats_imputed #check difference 

  
  return(stats_change)
}

#identify columns that had missing values 
imputed_col23 <- find_imputed_col(data_2023, imputed_data_23)
#Apply function to 2023
stats_difference_23 <- compare_stats(data_2023, imputed_data_23)

print(stats_difference_23)
write_csv(stats_difference_23, "/Users/marcinebessire/Desktop/project/Stats_Change_HalfMin_2023.csv")



