#load libraries required
library(tidyverse)

#load file path for 2023 and 2024 (fitlered data set)
file_path23 <- "/Users/marcinebessire/Desktop/project/Filtered_Data_2023.csv"
file_path24 <- "/Users/marcinebessire/Desktop/project/Filtered_Data_2024.csv"

#function to remove columns with >80% of NA
remove_MV <- function(filtered_data, year){
  #load filtered data
  data_df <- read.csv(filtered_data)
  
  #select metadata column 
  metadata_columns <- c("Name", "ID", "Year", "MonthDay", "Trial")
  
  #calculate NA in each column
  na_counts <- colSums(is.na(data_df))
  total_rows <- nrow(data_df)
  na_percent <- (na_counts/total_rows)*100 #percentage of NA
  
  #remove columns with more than 80% NA
  columns_to_keep <- names(na_percent[na_percent <= 20]) #keep the one with less or equal 80% NA
  final_cleaned_data <- data_df %>% select(all_of(columns_to_keep))
  
  #save final dataset 
  output_path <- paste0("/Users/marcinebessire/Desktop/project/Final_Data_", year, ".csv")
  write_csv(final_cleaned_data, output_path)
  
  #return final dataset
  return(final_cleaned_data)
}

#process filtered data, call function
#after filtering CV and removing MV and 
final_fitlered_data_2023 <- remove_MV(file_path23, "2023") 
final_fitlered_data_2024 <- remove_MV(file_path24, "2024")
