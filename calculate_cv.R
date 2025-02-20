#load libraries required
library(tidyverse)

# Part 1 -----
# Calculate CV

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

# Part 2 -----
# Remove columns with CV > 15-20%

#load CV files an cleaned data files

#2023
cv_file_2023 <- "/Users/marcinebessire/Desktop/project/CV_results_2023.csv"
data_file_2023 <- "/Users/marcinebessire/Desktop/project/cleaned_data_2023.csv"

#2024
cv_file_2024 <- "/Users/marcinebessire/Desktop/project/CV_results_2024.csv"
data_file_2024 <- "/Users/marcinebessire/Desktop/project/cleaned_data_2024.csv"

#function to filter CV 15-20% and return columns names 
filter_cv_columns <- function(cv_file, year){
  #load CV table
  cv_df <- read_csv(cv_file)
  
  #filter CV table wiht threshold 20%
  filtered_cv_df <- cv_df %>%
    filter(cv_df$`CV [%]` < 20) #filters out values greater than 20
  
  #save filter CV table 
  output_path <- paste0("/Users/marcinebessire/Desktop/project/Filtered_CV_", year, ".csv")
  write_csv(filtered_cv_df, output_path)
  
  return(filtered_cv_df$Column)
}

#process data of 2023
filtered_columns_2023 <- filter_cv_columns(cv_file_2023, "2023")

#process data of 2023
filtered_columns_2024 <- filter_cv_columns(cv_file_2024, "2024")

#output of filtering CV 2023: before 314, after fitlering wiht 20% 69 and with 15% 15 => kept 20%
#output of filtering CV 2024: before 660, after fitlering wiht 20% 77 and with 15% 27 => kept 20%

