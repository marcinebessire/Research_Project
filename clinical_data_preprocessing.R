#load necessary libraries 
library(readxl)
library(stringr)
library(dplyr)
library(lubridate) #to handle dates better
library(tidyverse)

#load excel sheet with data 
whole_data <- "/Users/marcinebessire/Desktop/project/clinical_data.xlsx"

#Get each sheet
#BAS data
BAS_data <- read_excel(whole_data, sheet = 2, .name_repair = "minimal")
#FAO data
FAO_data <- read_excel(whole_data, sheet = 3, .name_repair = "minimal")
#intact lipids data
intact_lipids_data <- read_excel(whole_data, sheet = 4, .name_repair = "minimal")

#--------------------------------------------
# Part 1: Clean dataset for Sheet 2 and 3 
# -------------------------------------------


#function for first step of cleaning data (remove column and make column names etc)
clean_data <- function(data){
  #remove columns before 3rd columns and remove first and second row
  data_cleaned <- data[, 3:ncol(data)]
  data_cleaned <- data_cleaned[3:nrow(data_cleaned),]
  
  #make second row the column name 
  colnames(data_cleaned) <- as.character(data_cleaned[1, ])  #assign second row as column names
  data_cleaned <- data_cleaned[-1, ]  #and remove the row that is now column names
  
  #remove sample ID row
  data_cleaned <- data_cleaned[2:nrow(data_cleaned),]
  
  return(data_cleaned)
  
}

#call function to clean data
BAS_data_cleaned <- clean_data(BAS_data)
FAO_data_cleaned <- clean_data(FAO_data)

#count how many E0 are in both datasets
E0_count_BAS <- nrow(BAS_data_cleaned[grep("E0", BAS_data_cleaned[[1]]), ]) #20
E0_count_FAO <- nrow(FAO_data_cleaned[grep("E0", FAO_data_cleaned[[1]]), ]) #20

#function to extract E0
e0_extraction <- function(data){
  #extract E0
  data_filtered <- data[grep("E0", data[[1]]), ]
  
  #rename column 1 to ID for further processing
  names(data_filtered)[1]<- paste("ID")
  
  return(data_filtered)
}

#call function to extract E0
BAS_data_cleaned <- e0_extraction(BAS_data_cleaned)
FAO_data_cleaned <- e0_extraction(FAO_data_cleaned)


#function to extract participant number "PX", date, time and visit 1 or visit 2
make_new_columns <- function(data){
  data <- data %>%
    mutate(
      Participant = str_extract(as.character(data[[1]]), "P[0-9]+"), #extract "P" followed by digits
      Date = str_extract(data[[1]], "[0-9]{6}"), #6 digit date
      Time = str_extract(data[[1]], "E[0-9]+") %>% str_replace("E", "") #extract number after E to get time
    ) %>%
    mutate(
      Day = str_sub(Date,1,2), #extract first two digits for day
      Month = str_sub(Date,3,4), #extract next two digits for month
      Year = paste0("20", str_sub(Date,5,6)), #extract year in format e.g. 2020
      MonthDaycalc = paste0(Month, Day), #create MonthDay column for ordering properly
      MonthDay = paste0(Month, "/", Day) #format for visualization
    ) %>%
    mutate(
      MonthDaycalc = as.numeric(MonthDaycalc), #convert to numeric for ordering
      Participant_Number = as.numeric(str_extract(Participant, "[0-9]+")) #extract numeric part of Participant for sorting
    ) %>%
    arrange(Participant_Number, MonthDaycalc) %>% #sort by date and participant
    group_by(Participant) %>%
    mutate(Visit = ifelse(row_number() == 1, "Visit 1", "Visit 2")) %>% #assign visit depending on date  %>%
    select(ID,Participant, MonthDay, Year, Time, Visit, everything()) %>%
    select(-Date, -Day, -Month, -MonthDaycalc, -Participant_Number) %>%
  
  return(data)
}

#call function 
BAS_data_extended <- make_new_columns(BAS_data_cleaned)
FAO_data_extended <- make_new_columns(FAO_data_cleaned)

#function to convert data to numeric and remove whole columns with all same vlaue
convert_columns_to_numeric <- function(data) {
  data[, 7:ncol(data)] <- lapply(data[, 7:ncol(data)], function(x) {
    suppressWarnings(as.numeric(x)) #convert to numeric and replace non-numeric values with NA
  })
  
  #remove columns with dame values 
  data <- data[, sapply(data, function(col) lenght(unique(na.omit(col))) > 1)]
  return(data)
}

convert_columns_to_numeric <- function(data) {
  #convert columns from 7th onwards to numeric
  data[, 7:ncol(data)] <- lapply(data[, 7:ncol(data)], function(x) {
    suppressWarnings(as.numeric(x)) # Convert to numeric and replace non-numeric values with NA
  })
  
  # Remove columns where all values are the same (including all NA)
  data <- data[, sapply(data, function(col) length(unique(na.omit(col))) > 1)]
  
  return(data)
}

#apply function to convert to numeric
BAS_data_final <- convert_columns_to_numeric(BAS_data_extended)
FAO_data_final <- convert_columns_to_numeric(FAO_data_extended)

#save to csv file 
write_csv(BAS_data_final, "/Users/marcinebessire/Desktop/project/BAS_data.csv")
write_csv(FAO_data_final, "/Users/marcinebessire/Desktop/project/FAO_data.csv")


#--------------------------------------------
# Part 2: Clean dataset for Sheet 4
# -------------------------------------------
