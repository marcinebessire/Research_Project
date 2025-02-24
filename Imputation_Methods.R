# Wrapper functions for Different Imputation Methods

#install packages required for imputation
#install.packages("missForest")
#install.packages("magrittr")
#install.packages("imputeLCMD")
#BiocManager::install("impute")
#BiocManager::install("pcaMethods")

#load necessary libraries 
library(missForest)
library(impute)
library(magrittr)
library(imputeLCMD)
library(tidyverse)

#load datasets 
final_data_2023 <- read_csv("/Users/marcinebessire/Desktop/project/Final_Data_2023.csv")

#Random Forest imputation
RF_wrapper <- function(data, ...) {
  numeric_data <- data %>% select(-c(Name, ID, Year, MonthDay, Trial)) #select only data (not metadata)
  imputed_data <- missForest(numeric_data, ...)[[1]]
  result <- bind_cols(data %>% select(Name, ID, Year, MonthDay, Trial), imputed_data)
  return(result)
}

RF_data <- RF_wrapper(final_data_2023)

final_data_2023
