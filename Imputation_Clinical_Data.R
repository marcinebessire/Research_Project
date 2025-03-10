#load necessary libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(missForest)
library(imputeLCMD)

#load original data
FAO_original <- read.csv("/Users/marcinebessire/Desktop/project/FAO_data.csv", check.names = FALSE)

#load data with missing values 
FAO_5pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_5pct.csv", check.names = FALSE)
FAO_10pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_10pct.csv", check.names = FALSE)
FAO_20pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_20pct.csv", check.names = FALSE)
FAO_25pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_25pct.csv", check.names = FALSE)
FAO_30pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_30pct.csv", check.names = FALSE)
FAO_40pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_40pct.csv", check.names = FALSE)

# ------------------------------------
# Part 1: Imputation methods 
# ------------------------------------

# ------------------------------------
# Part 1.1: Half-min Imputation 
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
# Part 1.2: KNN Imputation
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
# Part 1.: RF Imputation
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
# Part 1.4: QRILC Imputation
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

# ------------------------------------
# Part 2: Statistical Tests  
# ------------------------------------

# ------------------------------------
# Part 2.1: T-test
# ------------------------------------

t_test_func <- function(original, imputed) {
  #numeric columns
  numeric <- original[, 6:ncol(original)]
  
  #metabolte columns
  metabolite_cols <- colnames(numeric)
  
  #save resutls in dataframe
  results <- data.frame(
    Metabolite = metabolite_cols,
    p_value = numeric(length(metabolite_cols)),
    statistic = numeric(length(metabolite_cols))
  )
  
  #loop throuh each metabolite and run t-test 
  for (i in seq_along(metabolite_cols)) {
    metabolite <- metabolite_cols[i]
    
    #perform T-test
    t_test_result <- t.test(original[[metabolite]], imputed[[metabolite]], paired = FALSE, var.equal = FALSE)
    
    #save to dataframe
    results$p_value[i] <- t_test_result$p.value
    results$statistic[i] <- t_test_result$statistic
  }
  
  #adjust p-value
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  
  return(results)

}

#call t-test function
#Half-min
t_test_half_min_5pct <- t_test_func(FAO_original, Halfmin_5pct)
t_test_half_min_20pct <- t_test_func(FAO_original, Halfmin_20pct)
t_test_half_min_25pct <- t_test_func(FAO_original, Halfmin_25pct)
t_test_half_min_30pct <- t_test_func(FAO_original, Halfmin_30pct)
t_test_half_min_40pct <- t_test_func(FAO_original, Halfmin_40pct)
#KNN
t_test_KNN_5pct <- t_test_func(FAO_original, KNN_5pct)
t_test_KNN_10pct <- t_test_func(FAO_original, KNN_10pct)
t_test_KNN_20pct <- t_test_func(FAO_original, KNN_20pct)
t_test_KNN_25pct <- t_test_func(FAO_original, KNN_25pct)
t_test_KNN_30pct <- t_test_func(FAO_original, KNN_30pct)
t_test_KNN_40pct <- t_test_func(FAO_original, KNN_40pct)
#RF
t_test_RF_20pct <- t_test_func(FAO_original, RF_20pct)
t_test_RF_25pct <- t_test_func(FAO_original, RF_25pct)
t_test_RF_40pct <- t_test_func(FAO_original, RF_40pct)
#QRILC
t_test_QRILC_20pct <- t_test_func(FAO_original, QRILC_20pct)
t_test_QRILC_25pct <- t_test_func(FAO_original, QRILC_25pct)
t_test_QRILC_40pct <- t_test_func(FAO_original, QRILC_40pct)

# ------------------------------------
# Part 2.2: Wilcoxon rank-sum Test
# ------------------------------------

wilcoxon_func <- function(original, imputed) {
  #numeric columns
  numeric <- original[, 6:ncol(original)]
  
  #metabolte columns
  metabolite_cols <- colnames(numeric)
  
  #save resutls in dataframe
  results <- data.frame(
    Metabolite = metabolite_cols,
    p_value = numeric(length(metabolite_cols)),
    statistic = numeric(length(metabolite_cols))
  )
  
  #loop throuh each metabolite and run t-test 
  for (i in seq_along(metabolite_cols)) {
    metabolite <- metabolite_cols[i]
    
    #perform T-test
    wilcox_test_result <- wilcox.test(original[[metabolite]], imputed[[metabolite]], paired = FALSE, exact = FALSE)
    
    #save to dataframe
    results$p_value[i] <- wilcox_test_result$p.value
    results$statistic[i] <- wilcox_test_result$statistic
  }
  
  #adjust p-value
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  
  return(results)
  
}

#call wilcoxon function 
#Halfmin
wilcox_half_min_5pct <- wilcoxon_func(FAO_original, Halfmin_5pct) 
wilcox_half_min_20pct <- t_test_func(FAO_original, Halfmin_20pct)
wilcox_half_min_25pct <- t_test_func(FAO_original, Halfmin_25pct)
wilcox_half_min_30pct <- t_test_func(FAO_original, Halfmin_30pct)
wilcox_half_min_40pct <- wilcoxon_func(FAO_original, Halfmin_40pct)
#KNN
wilcox_KNN_5pct <- wilcoxon_func(FAO_original, KNN_5pct) 
wilcox_KNN_20pct <- t_test_func(FAO_original, KNN_20pct)
wilcox_KNN_25pct <- t_test_func(FAO_original, KNN_25pct)
wilcox_KNN_30pct <- t_test_func(FAO_original, KNN_30pct)
wilcox_KNN_40pct <- wilcoxon_func(FAO_original, KNN_40pct)
#RF
wilcox_RF_20pct <- wilcoxon_func(FAO_original, RF_20pct) 
wilcox_RF_25pct <- wilcoxon_func(FAO_original, RF_25pct) 
wilcox_RF_40pct <- wilcoxon_func(FAO_original, RF_40pct) 
#QRILC
wilcox_QRILC_20pct <- wilcoxon_func(FAO_original, QRILC_20pct) 
wilcox_QRILC_25pct <- wilcoxon_func(FAO_original, QRILC_25pct) 
wilcox_QRILC_40pct <- wilcoxon_func(FAO_original, QRILC_40pct) 

# ------------------------------------
# Part 2.3: Check Signficance of Tests
# ------------------------------------

#check for significance (p value < 0.05) 
significance <- function(data) {
  #view significant metabolites with BH adjusted p-value < 0.05
  data_significant <- data %>% 
    filter(adj_p_value < 0.05) %>% #usually 0.05 used 
    arrange(adj_p_value)
  
  return(data_significant)
}

#call significance function on T-test results
#Halfmin
signif_halfmin_5pct_t_test <- significance(t_test_half_min_5pct) #0/34
signif_halfmin_20pct_t_test <- significance(t_test_half_min_20pct) #0/34
signif_halfmin_25pct_t_test <- significance(t_test_half_min_25pct) #0/34
signif_halfmin_30pct_t_test <- significance(t_test_half_min_30pct) #0/34
signif_halfmin_40pct_t_test <- significance(t_test_half_min_40pct) #14/34
#KNN
signif_KNN_5pct_t_test <- significance(t_test_KNN_5pct) #0/34
signif_KNN_10pct_t_test <- significance(t_test_KNN_10pct) #0/34
signif_KNN_20pct_t_test <- significance(t_test_KNN_20pct) #0/34
signif_KNN_25pct_t_test <- significance(t_test_KNN_25pct) #0/34
signif_KNN_40pct_t_test <- significance(t_test_KNN_40pct) #7/34
#RF
signif_RF_20pct_t_test <- significance(t_test_RF_20pct) #0/34
signif_RF_25pct_t_test <- significance(t_test_RF_25pct) #0/34
signif_RF_40pct_t_test <- significance(t_test_RF_40pct) #0/34
#QRILC
signif_QRILC_20pct_t_test <- significance(t_test_QRILC_20pct) #0/34
signif_QRILC_25pct_t_test <- significance(t_test_QRILC_25pct) #0/34
signif_QRILC_40pct_t_test <- significance(t_test_QRILC_40pct) #0/34


#call significance function on wilcoxon results
#Halfmin
signif_halfmin_5pct_wilcox <- significance(wilcox_half_min_5pct) #0/34
signif_halfmin_20pct_wilcox <- significance(wilcox_half_min_20pct) #0/34
signif_halfmin_25pct_wilcox <- significance(wilcox_half_min_25pct) #0/34
signif_halfmin_30pct_wilcox <- significance(wilcox_half_min_30pct) #0/34
signif_halfmin_40pct_wilcox <- significance(wilcox_half_min_40pct) #21/34
#KNN
signif_KNN_20pct_wilcox <- significance(wilcox_half_min_20pct) #0/34
signif_KNN_25pct_wilcox <- significance(wilcox_half_min_25pct) #0/34
signif_KNN_40pct_wilcox <- significance(wilcox_half_min_40pct) #21/34
#RF
signif_RF_20pct_wilcox <- significance(wilcox_RF_20pct) #0/34
signif_RF_25pct_wilcox <- significance(wilcox_RF_25pct) #0/34
signif_RF_40pct_wilcox <- significance(wilcox_RF_40pct) #0/34
#QRILC
signif_QRILC_20pct_wilcox <- significance(wilcox_QRILC_20pct) #0/34
signif_QRILC_25pct_wilcox <- significance(wilcox_QRILC_25pct) #0/34
signif_QRILC_40pct_wilcox <- significance(wilcox_QRILC_40pct) #0/34

