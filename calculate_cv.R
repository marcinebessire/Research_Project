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
  
  #set inf values to NA (because are like missing values)
  df_numeric[df_numeric == Inf | df_numeric == -Inf] <- NA 
  
  #calculate coefficient of variance
  cv_values <- sapply(df_numeric, function(x) {
    if(all(is.na(x))) {
      return(NA_real_)
    }
    (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
  })
  
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
#cv_results_2024 <- calculate_cv(file_path24) do not do that

#print results 
print(cv_results_2023)
#print(cv_results_2024)

# Part 2 -----
# Remove columns with CV < 20%

#load CV files an cleaned data files

#2023
cv_file_2023 <- "/Users/marcinebessire/Desktop/project/CV_results_2023.csv"
cleaned_data_file_2023 <- "/Users/marcinebessire/Desktop/project/cleaned_data_2023.csv"

#2024
#cv_file_2024 <- "/Users/marcinebessire/Desktop/project/CV_results_2024.csv"
#cleaned_data_file_2024 <- "/Users/marcinebessire/Desktop/project/cleaned_data_2024.csv"

#function to filter CV 15-20% and return columns names 
filter_cv_columns <- function(cv_file, year){
  #load CV table
  cv_df <- read_csv(cv_file)
  
  #filter CV table with threshold 45%
  filtered_cv_df <- cv_df %>%
    filter(cv_df$`CV [%]` < 45) #filters out values greater than 20
  
  #save filter CV table 
  output_path <- paste0("/Users/marcinebessire/Desktop/project/Filtered_CV_", year, ".csv")
  write_csv(filtered_cv_df, output_path)
  
  return(filtered_cv_df$Column)
}

#process data of 2023
filtered_columns_2023 <- filter_cv_columns(cv_file_2023, "2023")

#process data of 2023
#filtered_columns_2024 <- filter_cv_columns(cv_file_2024, "2024")

#output of filtering CV 2023: before 314, after filtering with 20% 68 and with 15% 14 => kept 20%
#output of filtering CV 2024: before 660, after filtering with 20% 76 and with 15% 26 => kept 20%

# Part 3 -----
# Extract columns with CV < 20% in cleaned data files based on name

filter_data_by_columns <- function(cleaned_data_file, filtered_columns, year){
  #load cleaned dataset 
  cleaned_df <- read_csv(cleaned_data_file)
  
  #select again metadata columns 
  metadata_columns <- c("Name", "ID", "Year", "MonthDay", "Trial")
  
  #find names that match between filtered and clenead data
  matching_columns <- intersect (filtered_columns, colnames(cleaned_df))
  
  #for the output keep only metadata (as information) and selected columns from CV filtering
  filtered_data_df <- cleaned_df %>%
    select(all_of(metadata_columns), all_of(matching_columns))
  
  #save output
  output_path <- paste0("/Users/marcinebessire/Desktop/project/Filtered_Data_", year, ".csv")
  write_csv(filtered_data_df, output_path)
  
  #return the filtered dataset
  return(filtered_data_df)
  
}

#process data for 2023
filtered_data_2023 <- filter_data_by_columns(cleaned_data_file_2023, filtered_columns_2023, "2023") #should be 5 + 68 = 73 (correct)
#process data for 2024
filtered_data_2024 <- df_data_cleaned24 #do not filter according to CV ! 


