#load necessary libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(missForest)
library(imputeLCMD)

#load data with missing values 
FAO_5pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_5pct.csv", check.names = FALSE)
FAO_10pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_10pct.csv", check.names = FALSE)
FAO_20pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_20pct.csv", check.names = FALSE)
FAO_25pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_25pct.csv", check.names = FALSE)
FAO_30pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_30pct.csv", check.names = FALSE)
FAO_40pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_40pct.csv", check.names = FALSE)

# ------------------------------------
# Part 1: Half-min Imputation
# ------------------------------------

#function for Half-minimum imputation
half_min_imputation <- function(data){
  #create copy
  data_imputed <- data 

  #metadata
  meta_data <- data_imputed[, 1:5]
  #select only numeric data
  numeric <- data_imputed[, 6:ncol(data_imputed)]
  
  #loop through column
  for (col in names(numeric)) { 
    min_val <- min(numeric[[col]], na.rm = TRUE) #find min value 
    numeric[[col]][is.na(numeric[[col]])] <- 0.5 * min_val #calculate half min and replace NA
  }
  
  final_df <- cbind(meta_data, numeric)
  
  return(final_df)
}

#call function for Half-min imputation
Halfmin_5pct <- half_min_imputation(FAO_5pct)
Halfmin_10pct <- half_min_imputation(FAO_10pct)
Halfmin_20pct <- half_min_imputation(FAO_20pct)
Halfmin_25pct <- half_min_imputation(FAO_25pct)
Halfmin_30pct <- half_min_imputation(FAO_30pct)
Halfmin_40pct <- half_min_imputation(FAO_40pct)

# ------------------------------------
# Part 2: KNN Imputation
# ------------------------------------

KNN_imputation <- function(data) {
  #make copy 
  data_imputed <- data
  
  #select only numeric data
  numeric <- data_imputed[, 6:ncol(data_imputed)]
  
  #metadata
  meta_data <- data_imputed[, 1:5]
  
  #transform first into matrix
  imputed_data <- impute.knn(as.matrix(t(numeric)), rowmax = 0.5, colmax = 1) 
  
  #transform back into dataframe
  imputed_data <- as.data.frame(t(imputed_data$data))
  
  final_df <- cbind(meta_data, imputed_data)
  
  return(final_df)
}

#call function for KNN
KNN_5pct <- KNN_imputation(FAO_5pct)
KNN_10pct <- KNN_imputation(FAO_10pct)
KNN_20pct <- KNN_imputation(FAO_20pct)
KNN_25pct <- KNN_imputation(FAO_25pct)
KNN_30pct <- KNN_imputation(FAO_30pct)
KNN_40pct <- KNN_imputation(FAO_40pct)

# ------------------------------------
# Part 3: RF Imputation
# ------------------------------------

#make function for RF
RF_imputation <- function(data) {
  #make a copy 
  data_copy <- data 

  #select only numeric data
  numeric <- data_copy[, 6:ncol(data_copy)]
  
  #metadata
  meta_data <- data_copy[, 1:5]
  
  #apply miss forest on data to impute MV
  imputed_data <- missForest(numeric, maxiter = 10, ntree = 100)
  imputed_df <- as.data.frame(imputed_data$ximp)
  
  #return final dataframe
  final_df <- cbind(meta_data, imputed_df)
}

#call function for RF
RF_5pct <- RF_imputation(FAO_5pct)
RF_10pct <- RF_imputation(FAO_10pct)
RF_20pct <- RF_imputation(FAO_20pct)
RF_25pct <- RF_imputation(FAO_25pct)
RF_30pct <- RF_imputation(FAO_30pct)
RF_40pct <- RF_imputation(FAO_40pct)

# ------------------------------------
# Part 4: QRILC Imputation
# ------------------------------------

#QRILC imputation 
QRILC_impuation <- function(data) {
  #copy data
  data_copy <- data
  
  #select only numeric data
  numeric <- data_copy[, 6:ncol(data_copy)]
  
  #metadata
  meta_data <- data_copy[, 1:5]
  
  #transfer to log data
  log_data <- log(numeric + 1e-6)  #add small constant to avoid log(0)

  #make imputation
  imputed_data <- impute.QRILC(log_data)[[1]] #returns list, extract imputed matrix
  
  #convert back to exp 
  exp_imputed_data <- exp(imputed_data) - 1e-6
  
  #save as datarame
  imputed_df <- as.data.frame(exp_imputed_data)
  
  #return final dataframe
  final_df <- cbind(meta_data, imputed_df)
  return(final_df)
}

#call function for QRILC
QRILC_5pct <- QRILC_impuation(FAO_5pct)
QRILC_10pct <- QRILC_impuation(FAO_10pct)
QRILC_20pct <- QRILC_impuation(FAO_20pct)
QRILC_25pct <- QRILC_impuation(FAO_25pct)
QRILC_30pct <- QRILC_impuation(FAO_30pct)
QRILC_40pct <- QRILC_impuation(FAO_40pct)


