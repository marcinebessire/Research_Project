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
  numeric_or_factor <- data %>% select(where(~ is.numeric(.x) || is.factor(.x)))
  result <- missForest(numeric_or_factor, ...)[[1]]
  return (result)
}

RF_data <- RF_wrapper(final_data_2023)

final_data_2023
