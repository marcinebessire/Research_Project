# Load necessary library
library(tidyverse) #data visualization, collection of packages
library(lubridate) #to handle dates better

# Part 1 -------
# Read Data, make unique column names (handle duplication)

# Read the CSV file
file_path <- "/Users/marcinebessire/Desktop/project/Result_out_MS2_filtered_20250218.csv"
df_data24 <- read_csv(file_path, name_repair = "minimal") # Keep duplicate names

#change change duplicated names to .2, .3 etc.
name_count <- table(colnames(df_data24)) 
seen_count <- list() 
col_names <- colnames(df_data24)

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
colnames(df_data24) <- col_names

# Part 2 ------
# Expand Name column 
# Merge column name and first row of data to get unique names

df_data_cleaned24 <- df_data24 %>%
  mutate( #mutate to create or edit existing columns in a dataframe
    #extract date from anywhere in the Name column (because not the same)
    Whole_Date = str_extract(Name, "\\d{8}"), #\\d for matching any digits
    Whole_Date = as.Date(Whole_Date, format="%Y%m%d"), 
    Year = year(Whole_Date), 
    MonthDay = format(Whole_Date, "%m-%d"),
    
    #extract only the meaningful part of ID 
    ID = str_remove(Name, ".*?_NIST_?"),  #remove everything before and including "NIST_"
    ID = str_remove(ID, "^DI_"),  #remove "DI_" prefix if present, ^ for the beginning 
    ID = str_remove(ID, "\\d{8}"),  #remove date if it appears again (for the wrongly ordered)
    ID = str_replace(ID, "^_|_$", ""),  #remove underscores
    
    #assign trial numbers based on MonthDay
    Trial_number = dense_rank(MonthDay),  
    Trial = paste0("Trial ", Trial_number)  
  ) %>%
  select(Name, ID, Year, MonthDay, Trial, everything())


# Part 3 -------
#now merge the column name to get unique names

#define the values to ignore
ignore_values <- c("Species", "NA", "Trial NA", "MS1") 
new_colnames <- names(df_data_cleaned24)

#loop through each column index
for (i in seq_along(names(df_data_cleaned24))) {
  first_row_val <- as.character(df_data_cleaned24[1, i])  #get first row value
  column_name <- names(df_data_cleaned24)[i]  #get column name
  
  #if the value is not in the ignore list, merge it with column name
  if (!(is.na(first_row_val) || first_row_val %in% ignore_values)) {
    new_colnames[i] <- paste0(column_name, "_", first_row_val)
  }
}

#assign the new names
colnames(df_data_cleaned24) <- new_colnames
#remove first row after merging
df_data_cleaned24 <- df_data_cleaned24[-1, ] 

# Part 4 ------
# Convert data columns (keeping metadata unchanged)
info_cols <- c("Name", "ID", "Year", "MonthDay", "Trial") 

df_data_cleaned24 <- df_data_cleaned24 %>%
  mutate(across(-all_of(info_cols), ~ suppressWarnings(as.numeric(.))))  

print(df_data_cleaned24)

#write new CSV
write_csv(df_data_cleaned24, "/Users/marcinebessire/Desktop/project/cleaned_data_2024.csv")


