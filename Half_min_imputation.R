# Load required library
library(tidyverse)

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
