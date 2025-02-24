# Load necessary library
library(tidyverse) #data visualization, collection of packagess
library(lubridate) #to handle dats better

# Part 1 -------
# Read Data, make unique column names (handle duplications)

# Read the CSV file
file_path <- "/Users/marcinebessire/Desktop/project/Result_out_MS2_filtered_20250217.csv"
df_data23 <- read_csv(file_path, name_repair = "minimal") #keep duplicates

# #change change duplicated names to .2, .3 etc. 
# name_count <- table(colnames(df_data23)) #make table with count of each name
# seen_count <- list() #to check how many names have appeared
# col_names <- colnames(df_data23) #store current column names
# 
# #iterate through column names and rename duplicates
# for (i in seq_along(col_names)){
#   name <- col_names[i]
#   
#   #check if there are duplicates
#   if (name_count[name] > 1){
#     if (!name %in% names(seen_count)){
#       seen_count[[name]] <- 1 #because first occurrence, keep original
#     } else { #if not first occurence update the duplication
#       seen_count[[name]] <- seen_count[[name]] + 1 #increment
#       col_names[i] <- paste0(name,".", seen_count[[name]]) #append the seen count
#     }
#   } 
# }
# 
# #assign the new column names 
# colnames(df_data23) <- col_names

# Part 2 -------
#now merge the column name to get unique names

#make a copy firs
df_data23_copy <- df_data23

#define the values to ignore
ignore_values <- c("Species", "NA", "Trial NA", "MS1")

#new column names
new_colnames <- names(df_data23_copy)

#loop through each column index
for (i in seq_along(names(df_data23_copy))) {
  first_row_val <- as.character(df_data23_copy[1, i])  #get first row value
  column_name <- names(df_data23_copy)[i]  #get column name
  
  #if the value is not in the ignore list, merge it with column name
  if (!(is.na(first_row_val) || first_row_val %in% ignore_values)) {
    new_colnames[i] <- paste0(column_name, "_", first_row_val)
  }
}

#assign the new names
colnames(df_data23_copy) <- new_colnames
#remove first row after merging 
df_data23_copy <- df_data23_copy[-1, ]

# Part 3 ------
# Expand Name column 
# Merge column name and first row of data to get unique names

#extract Year and Date from the Name column
df_data_cleaned23 <- df_data23_copy %>%
  mutate(
    ID = str_remove(Name, "^[^_]+_[^_]+_"), 
    Whole_Date = as.Date(str_extract(Name, "^\\d{8}"), format="%Y%m%d"), #YYYYMMDD
    Year = str_sub(Whole_Date, 1, 4), #YYYY
    MonthDay = format(Whole_Date, "%m-%d"),
    Trial_number = dense_rank(MonthDay), #rank based on MonthDay
    Trial = paste0("Trial ", Trial_number), #Set Trial number for each Trial
  ) %>%
  select(Name, ID, Year, MonthDay, Trial, everything(), -Whole_Date, -Trial_number)

# Part 3 ------
# Convert data columns (keeping metadata unchanged)
info_cols <- c("Name", "ID", "Year", "MonthDay", "Trial") 

df_data_cleaned23 <- df_data_cleaned23 %>%
  mutate(across(-all_of(info_cols), ~ suppressWarnings(as.numeric(.))))  

#replace NA and Inf values with empty strings before writing the CSV
df_data_cleaned23 <- df_data_cleaned23 %>%
  mutate(across(everything(), ~if_else(is.na(.) | is.infinite(.), "", as.character(.))))


print(df_data_cleaned23)

#write new csv
write_csv(df_data_cleaned23, "/Users/marcinebessire/Desktop/project/cleaned_data_2023.csv")


