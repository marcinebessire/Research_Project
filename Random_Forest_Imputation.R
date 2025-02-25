#load necessary library
library(missForest)
library(dplyr)

#load common metabolite dataframe (2023 and 2024) 
data_23 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites23.csv", check.names = FALSE)
data_24 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites24.csv", check.names = FALSE)

#Only use nuemric columns for imputation
numeric_23 <- data_23[, 6:ncol(data_23)]
numeric_24 <- data_24[, 6:ncol(data_24)]

#make function for RF
RF_imputation <- function(data, ...) {
  #apply miss forest on data to impute MV
  imputed_data <- missForest(data, maxiter = 10, ntree = 100)
  
  #return imputed data
  return(imputed_data)
  
}

#Apply RF function on data 
RF_data23 <- RF_imputation(numeric_23)
RF_data24 <- RF_imputation(numeric_24)

#extract imputed data and save as df
imputed_RF_23 <- as.data.frame(RF_data23$ximp)
imputed_RF_24 <- as.data.frame(RF_data24$ximp)
