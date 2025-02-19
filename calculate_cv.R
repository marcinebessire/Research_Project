#load libraries required
library(tidyverse)

#load data of 2023 and 2034 (cleaned data)
file_path23 <- "/Users/marcinebessire/Desktop/project/cleaned_data_2023.csv"
file_path24 <- "/Users/marcinebessire/Desktop/project/cleaned_data_2024.csv"

#write function to caluclate cv for each numeric column in each file
calculate_cv <- function(file_path){
  #read csv file
  df <- read_csv(file_path)
  
  #extract filename and date for output file 
  filename <- basename(file_path)
  filedate <- str_extract(filename, "\\d{4}" ) #extract year
  
  #exclude non-data columns
  df_data <- df %>% select(-c(Name, ID, Year, MonthDay, Trial, Whole_Date, Trial_number))
  
  #keep only numeric columns
  df_numeric <- df_data %>% select(where(is.numeric))
  
  #calculate coefficient of variance
  cv_values <- sapply(df_data, function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100)
  
  #convert results into dataframe
  cv_df <- data.frame(Column = names(cv_values), CV_Percentage = cv_values)
  
  #rename column
  colnames(cv_df)[2] <- "CV [%]"
  
  #save results to csv file
  output_path <- paste0("/Users/marcinebessire/Desktop/project/CV_results_", filedate, ".csv")
  write_csv(cv_df, output_path)
  
  #return CV dataframe
  return(cv_df)
}

#run function for both
cv_results_2023 <- calculate_cv(file_path23)
cv_results_2024 <- calculate_cv(file_path24)

#print results 
print(cv_results_2023)
print(cv_results_2024)

