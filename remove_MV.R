#load libraries required
library(tidyverse)
library(dplyr)


# Part 1 -------
# Remove Missing Values 
# Choose cutoff of 20%

#load file path for 2023 and 2024 (filtered data set) for the one who where filtered
file_path23cv <- "/Users/marcinebessire/Desktop/project/CV_Filtered_Data_2023.csv"
file_path24cv <- "/Users/marcinebessire/Desktop/project/CV_Filtered_Data_2024.csv"

#load data of 2023 and 2034 (cleaned data) for the one w/o fitlering CV
#file_path23 <- "/Users/marcinebessire/Desktop/project/cleaned_data_2023.csv"
#file_path24 <- "/Users/marcinebessire/Desktop/project/cleaned_data_2024.csv"
#cut_file_path24 <- "/Users/marcinebessire/Desktop/project/cut_cleaned_data_2024.csv"

#function to remove columns with >80% of NA
remove_MV <- function(filtered_data, year){
  #load filtered data
  data_df <- read.csv(filtered_data, check.names = FALSE)
  
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
  output_path <- paste0("/Users/marcinebessire/Desktop/project/CV_Final_Data_", year, ".csv") #ajudst if needed
  write_csv(final_cleaned_data, output_path)
  
  #return final dataset
  return(final_cleaned_data)
}

#process filtered data, call function
#final_filtered_data_2023 <- remove_MV(file_path23, "2023") 
#final_filtered_data_2024 <- remove_MV(file_path24, "2024")
#final_cut_filtered_data_2024 <- remove_MV(cut_file_path24, "2024.2")

#process filtered CV data
filtered_data23 <- remove_MV(file_path23cv, "2023")
filtered_data24 <- remove_MV(file_path24cv, "2024")

#Check how many column names are the same in both datasets and which ones

#columns to exclude
exclude_cols <- c('Name', 'ID', 'Year', 'MonthDay', 'Trial')

#identify same columns
#cols_df1 <- setdiff(colnames(final_filtered_data_2023), exclude_cols)
#cols_df2 <- setdiff(colnames(final_filtered_data_2024), exclude_cols)
#cols_df3 <- setdiff(colnames(final_cut_filtered_data_2024), exclude_cols) #for cut data

#identify columns for cv filtered data
cols_df1 <- setdiff(colnames(filtered_data23), exclude_cols)
cols_df2 <- setdiff(colnames(filtered_data24), exclude_cols)

#find common columns
common_cols1 <- intersect(cols_df1, cols_df2)
#common_cols2 <- intersect(cols_df1, cols_df3) #for cut data

#output results
cat("Number of common columns:", length(common_cols1), "\n") #84 for no CV and with CV of 30 it is 60 columns
cat("Common columns:", paste(common_cols1, collapse=", "))
 
#output results for cut data
#cat("Number of common columns:", length(common_cols2), "\n") #19 for cut data
#cat("Common columns:", paste(common_cols2, collapse=", "))

#dataframe 
final_col1 <- c(exclude_cols, common_cols1)
#final_col2 <- c(exclude_cols, common_cols2)

#make new dataframes for each year with those 60 columns + metadata (with CV of 30)
common_metabolites_23cv <- filtered_data23 %>% select(all_of(final_col1))
common_metabolites_24cv <- filtered_data24 %>% select(all_of(final_col1))

#make new dataframes for each year with those 84 columns + metadata
#common_metabolites_23 <- final_filtered_data_2023 %>% select(all_of(final_col1))
#common_metabolites_24 <- final_filtered_data_2024 %>% select(all_of(final_col1))

#for cut data
#common_metabolites_23_2 <- final_filtered_data_2023 %>% select(all_of(final_col2))
#common_metabolites_24_2 <- final_filtered_data_2024 %>% select(all_of(final_col2))

#save common data with CV of 30%
write_csv(common_metabolites_23cv, "/Users/marcinebessire/Desktop/project/CV_Common_Metabolites23.csv")
write_csv(common_metabolites_24cv, "/Users/marcinebessire/Desktop/project/CV_Common_Metabolites24.csv")

#save common data 
#write_csv(common_metabolites_23, "/Users/marcinebessire/Desktop/project/Common_Metabolites23.csv")
#write_csv(common_metabolites_24, "/Users/marcinebessire/Desktop/project/Common_Metabolites24.csv")

#save common data for cut 2024
#write_csv(common_metabolites_23_2, "/Users/marcinebessire/Desktop/project/Cut_Common_Metabolites23.csv")
#write_csv(common_metabolites_24_2, "/Users/marcinebessire/Desktop/project/Cut_Common_Metabolites24.csv")

# Part 2 ----
# Asses normality of data with Shapiro-Wilk test

#select only metabolite columns (starting from column 6)
data_columns23 <- common_metabolites_23[, 6:ncol(common_metabolites_23)]
data_columns24 <- common_metabolites_24[, 6:ncol(common_metabolites_24)]

# Part 3.2 -------
# Check normality of data 

#check normality using shapiro.test (becuase t-test assumes normality)
shapiro_results23 <- apply(data_columns23, 2, function(x) shapiro.test(x)$p.value)
shapiro_results24 <- apply(data_columns24, 2, function(x) shapiro.test(x)$p.value)


#convert to a dataframe for easy viewing
shapiro_df23 <- data.frame(Metabolite = names(shapiro_results23), p_value = shapiro_results23) #total 84 metabolites
shapiro_df24 <- data.frame(Metabolite = names(shapiro_results24), p_value = shapiro_results24)

#if p-value < 0.05 then not normal distribution
non_normal_count23 <- sum(shapiro_df23$p_value < 0.05)
non_normal_count23 #59 metabolites are non-normal distributed 
non_normal_count24 <- sum(shapiro_df24$p_value < 0.05)
non_normal_count24 #56 metabolites are non-normal distributed 


