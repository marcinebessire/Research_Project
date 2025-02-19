# Load necessary library
library(tidyverse) #data visualization, collection of packagess
library(lubridate) #to handle dats better

# Part 1 -------
# Read Data, make unique column names (handle duplications)

# Read the CSV file
file_path <- "/Users/marcinebessire/Desktop/project/Result_out_MS2_filtered_20250218.csv"
df_data <- read_csv(file_path, name_repair = "minimal") #keep duplicates

#change change duplicated names to .2, .3 etc. 
name_count <- table(colnames(df_data)) #make table with count of each name
seen_count <- list() #to check how many names have appeared
col_names <- colnames(df_data) #store current column names

#iterate through column names and rename duplicates
for (i in seq_along(col_names)){
  name <- col_names[i]
  
  #check if there are duplicates
  if (name_count[name] > 1){
    if (!name %in% names(seen_count)){
      seen_count[[name]] <- 1 #because first occurrence, keep original
    } else { #if not first occurence update the duplication
      seen_count[[name]] <- seen_count[[name]] + 1 #increment
      col_names[i] <- paste0(name,".", seen_count[[name]]) #append the seen count
    }
  } 
}

#assign the new column names 
colnames(df_data) <- col_names

#becuase some Names are different identify them and reverse order so all start with number 

# Part 2 ------
# Expand Name column 
# Merge column name and first row of data to get unique names

#extract Year and Date from the Name column
df_data_cleaned <- df_data %>%
  mutate(
    ID = str_remove(Name, "^[^_]+_[^_]+_"), 
    Whole_Date = as.Date(str_extract(Name, "^\\d{8}"), format="%Y%m%d"), #YYYYMMDD
    Year = str_sub(Whole_Date, 1, 4), #YYYY
    MonthDay = format(Whole_Date, "%m-%d"),
    Trial_number = dense_rank(MonthDay), #rank based on MonthDay
    Trial = paste0("Trial ", Trial_number), #Set Trial number for each Trial
  ) %>%
  select(Name, ID, Year, MonthDay, Trial, everything())

#now merge the column name to get unique names
#define the values to ignore
ignore_values <- c("Species", "NA", "Trial NA", "MS1")

#new column names
new_colnames <- names(df_data_cleaned)

#loop through each column index
for (i in seq_along(names(df_data_cleaned))) {
  first_row_val <- as.character(df_data_cleaned[1, i])  #get first row value
  column_name <- names(df_data_cleaned)[i]  #get column name
  
  #if the value is not in the ignore list, merge it with column name
  if (!(is.na(first_row_val) || first_row_val %in% ignore_values)) {
    new_colnames[i] <- paste0(column_name, "_", first_row_val)
  }
}

#assign the new names
colnames(df_data_cleaned) <- new_colnames
#remove first row after merging 
df_data_cleaned <- df_data_cleaned[-1, ]
df_data_cleaned


#write new csv
write_csv(df_data_cleaned, "/Users/marcinebessire/Desktop/project/cleaned_data_2024.csv")
