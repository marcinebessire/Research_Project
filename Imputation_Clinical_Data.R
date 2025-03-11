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
# Part 1.3: RF Imputation
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
# Part 2.1: Shapiro-Wilk test
# ------------------------------------

#check normality using shapiro.test 
shapiro_test <- function(data){
  #copy data
  data_copy <- data
  
  #select only numeric data
  numeric <- data_copy[, 6:ncol(data_copy)]
  
  #do shapiro test
  shapiro_results <- apply(numeric, 2, function(x) shapiro.test(x)$p.value)
  
  #concert to dataframe
  shapiro_df <- data.frame(Metabolite = names(shapiro_results), p_value = shapiro_results) 
  
  #print how many are significant = if p-value < 0.05 then not normal distribution
  significant_count <- sum(shapiro_df$p_value < 0.05)
  cat("Number of non-normally distributed metabolites:", significant_count, "\n")
  
  #return dataframe
  return(shapiro_df)
}

#call shaprio function
shapiro_original <- shapiro_test(FAO_original) #16/34 are significant 


# ------------------------------------
# Part 2.2: T-test
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
t_test_half_min_10pct <- t_test_func(FAO_original, Halfmin_10pct)
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
t_test_RF_5pct <- t_test_func(FAO_original, RF_5pct)
t_test_RF_10pct <- t_test_func(FAO_original, RF_10pct)
t_test_RF_20pct <- t_test_func(FAO_original, RF_20pct)
t_test_RF_25pct <- t_test_func(FAO_original, RF_25pct)
t_test_RF_30pct <- t_test_func(FAO_original, RF_30pct)
t_test_RF_40pct <- t_test_func(FAO_original, RF_40pct)
#QRILC
t_test_QRILC_5pct <- t_test_func(FAO_original, QRILC_5pct)
t_test_QRILC_10pct <- t_test_func(FAO_original, QRILC_10pct)
t_test_QRILC_20pct <- t_test_func(FAO_original, QRILC_20pct)
t_test_QRILC_25pct <- t_test_func(FAO_original, QRILC_25pct)
t_test_QRILC_30pct <- t_test_func(FAO_original, QRILC_30pct)
t_test_QRILC_40pct <- t_test_func(FAO_original, QRILC_40pct)

#Store results in dataframe
df_t_test <- bind_rows(
  t_test_half_min_5pct %>% mutate(Method = "Halfmin", Missing = "5%"),
  t_test_half_min_10pct %>% mutate(Method = "Halfmin", Missing = "10%"),
  t_test_half_min_20pct %>% mutate(Method = "Halfmin", Missing = "20%"),
  t_test_half_min_25pct %>% mutate(Method = "Halfmin", Missing = "25%"),
  t_test_half_min_30pct %>% mutate(Method = "Halfmin", Missing = "30%"),
  t_test_half_min_40pct %>% mutate(Method = "Halfmin", Missing = "40%"),
  t_test_KNN_5pct %>% mutate(Method = "KNN", Missing = "5%"),
  t_test_KNN_10pct %>% mutate(Method = "KNN", Missing = "10%"),
  t_test_KNN_20pct %>% mutate(Method = "KNN", Missing = "20%"),
  t_test_KNN_25pct %>% mutate(Method = "KNN", Missing = "25%"),
  t_test_KNN_30pct %>% mutate(Method = "KNN", Missing = "30%"),
  t_test_KNN_40pct %>% mutate(Method = "KNN", Missing = "40%"),
  t_test_RF_5pct %>% mutate(Method = "RF", Missing = "5%"),
  t_test_RF_10pct %>% mutate(Method = "RF", Missing = "10%"),
  t_test_RF_20pct %>% mutate(Method = "RF", Missing = "20%"),
  t_test_RF_25pct %>% mutate(Method = "RF", Missing = "25%"),
  t_test_RF_30pct %>% mutate(Method = "RF", Missing = "30%"),
  t_test_RF_40pct %>% mutate(Method = "RF", Missing = "40%"),
  t_test_QRILC_5pct %>% mutate(Method = "QRILC", Missing = "5%"),
  t_test_QRILC_10pct %>% mutate(Method = "QRILC", Missing = "10%"),
  t_test_QRILC_20pct %>% mutate(Method = "QRILC", Missing = "20%"),
  t_test_QRILC_25pct %>% mutate(Method = "QRILC", Missing = "25%"),
  t_test_QRILC_30pct %>% mutate(Method = "QRILC", Missing = "30%"),
  t_test_QRILC_40pct %>% mutate(Method = "QRILC", Missing = "40%")
)


# ------------------------------------
# Part 2.3: Wilcoxon rank-sum Test
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
wilcox_half_min_10pct <- wilcoxon_func(FAO_original, Halfmin_10pct) 
wilcox_half_min_20pct <- t_test_func(FAO_original, Halfmin_20pct)
wilcox_half_min_25pct <- t_test_func(FAO_original, Halfmin_25pct)
wilcox_half_min_30pct <- t_test_func(FAO_original, Halfmin_30pct)
wilcox_half_min_40pct <- wilcoxon_func(FAO_original, Halfmin_40pct)
#KNN
wilcox_KNN_5pct <- wilcoxon_func(FAO_original, KNN_5pct) 
wilcox_KNN_10pct <- wilcoxon_func(FAO_original, KNN_10pct) 
wilcox_KNN_20pct <- t_test_func(FAO_original, KNN_20pct)
wilcox_KNN_25pct <- t_test_func(FAO_original, KNN_25pct)
wilcox_KNN_30pct <- t_test_func(FAO_original, KNN_30pct)
wilcox_KNN_40pct <- wilcoxon_func(FAO_original, KNN_40pct)
#RF
wilcox_RF_5pct <- wilcoxon_func(FAO_original, RF_5pct) 
wilcox_RF_10pct <- wilcoxon_func(FAO_original, RF_10pct) 
wilcox_RF_20pct <- wilcoxon_func(FAO_original, RF_20pct) 
wilcox_RF_25pct <- wilcoxon_func(FAO_original, RF_25pct) 
wilcox_RF_30pct <- wilcoxon_func(FAO_original, RF_30pct) 
wilcox_RF_40pct <- wilcoxon_func(FAO_original, RF_40pct) 
#QRILC
wilcox_QRILC_5pct <- wilcoxon_func(FAO_original, QRILC_5pct) 
wilcox_QRILC_10pct <- wilcoxon_func(FAO_original, QRILC_10pct) 
wilcox_QRILC_20pct <- wilcoxon_func(FAO_original, QRILC_20pct) 
wilcox_QRILC_25pct <- wilcoxon_func(FAO_original, QRILC_25pct) 
wilcox_QRILC_30pct <- wilcoxon_func(FAO_original, QRILC_30pct) 
wilcox_QRILC_40pct <- wilcoxon_func(FAO_original, QRILC_40pct) 

# ------------------------------------
# Part 2.4: Wilcoxon singed-rank Test (Visit 1 vs 2)
# ------------------------------------

#function to compare significance between Visit 1 and Visit 2
compare_visits_wilcoxon <- function(data) {
  #numeric columns
  numeric <- data[, 6:ncol(data)]
  
  #metabolte columns
  metabolite_cols <- colnames(numeric)
  
  #save resutls in dataframe
  results <- data.frame(
    Metabolite = metabolite_cols,
    p_value = numeric(length(metabolite_cols)),
    statistic = numeric(length(metabolite_cols))
  )
  
  for (i in seq_along(metabolite_cols)) {
    metabolite <- metabolite_cols[i]
    
    #select particiaptn who have both visit 1 and 2 
    paired_data <- data %>%
      filter(Visit %in% c("Visit 1", "Visit 2")) %>%
      select(Participant, Visit, all_of(metabolite)) %>%
      spread(Visit, metabolite) #same as pivot_wider 
  
    #only compete pairs (no NA)
    paired_data <- paired_data %>% filter(!is.na(`Visit 1`) & !is.na(`Visit 2`))
    
    #perform wilcoxon singed-rank test
    if (nrow(paired_data) > 2) {
      #use paired here because same patient
      wilcox_test_result <- wilcox.test(paired_data$`Visit 1`, paired_data$`Visit 2`, paired = TRUE, exact = FALSE)
      
      #store resuls
      results$p_value[i] <- wilcox_test_result$p.value
      results$statistic[i] <- wilcox_test_result$statistic
    } else {
      #if not enough data
      results$p_value[i] <- NA
      results$statistic[i] <- NA
    }
  }
  #adjust p-value using BH
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  
  return(results)
  
}

#call function to compare Visit 1 vs. Visit 2
#original 
visit_original_res <- compare_visits_wilcoxon(FAO_original)
visit_origininal_signif <- significance(visit_original_res) #16/34 significant
#Halfmin
visit_Halfmin5pct_res <- compare_visits_wilcoxon(Halfmin_5pct)
visit_Halfmin5pct_signif <- significance(visit_Halfmin5pct_res) #13
visit_Halfmin10pct_res <- compare_visits_wilcoxon(Halfmin_10pct)
visit_Halfmin10pct_signif <- significance(visit_Halfmin10pct_res) #0/34
visit_Halfmin20pct_res <- compare_visits_wilcoxon(Halfmin_20pct)
visit_Halfmin20pct_signif <- significance(visit_Halfmin20pct_res) #0/34
visit_Halfmin25pct_res <- compare_visits_wilcoxon(Halfmin_25pct)
visit_Halfmin25pct_signif <- significance(visit_Halfmin25pct_res) #0/34
visit_Halfmin30pct_res <- compare_visits_wilcoxon(Halfmin_30pct)
visit_Halfmin30pct_signif <- significance(visit_Halfmin30pct_res) #0/34
visit_Halfmin40pct_res <- compare_visits_wilcoxon(Halfmin_40pct)
visit_Halfmin40pct_signif <- significance(visit_Halfmin40pct_res) #0/34
#KNN
visit_KNN5pct_res <- compare_visits_wilcoxon(KNN_5pct)
visit_KNN5pct_signif <- significance(visit_KNN5pct_res) #14/34
visit_KNN10pct_res <- compare_visits_wilcoxon(KNN_10pct)
visit_KNN10pct_signif <- significance(visit_KNN10pct_res) #14/34
visit_KNN20pct_res <- compare_visits_wilcoxon(KNN_20pct)
visit_KNN20pct_signif <- significance(visit_KNN20pct_res) #10/34
visit_KNN25pct_res <- compare_visits_wilcoxon(KNN_25pct)
visit_KNN25pct_signif <- significance(visit_KNN25pct_res) #6/34
visit_KNN30pct_res <- compare_visits_wilcoxon(KNN_30pct)
visit_KNN30pct_signif <- significance(visit_KNN30pct_res) #8/34
visit_KNN40pct_res <- compare_visits_wilcoxon(KNN_40pct)
visit_KNN40pct_signif <- significance(visit_KNN40pct_res) #0/34
#RF
visit_RF5pct_res <- compare_visits_wilcoxon(RF_5pct)
visit_RF5pct_signif <- significance(visit_RF5pct_res) #15/34
visit_RF10pct_res <- compare_visits_wilcoxon(RF_10pct)
visit_RF10pct_signif <- significance(visit_RF10pct_res) #17/34
visit_RF20pct_res <- compare_visits_wilcoxon(RF_20pct)
visit_RF20pct_signif <- significance(visit_RF20pct_res) #13/34
visit_RF25pct_res <- compare_visits_wilcoxon(RF_25pct)
visit_RF25pct_signif <- significance(visit_RF25pct_res) #15/34
visit_RF30pct_res <- compare_visits_wilcoxon(RF_30pct)
visit_RF30pct_signif <- significance(visit_RF30pct_res) #16/34
visit_RF40pct_res <- compare_visits_wilcoxon(RF_40pct)
visit_RF40pct_signif <- significance(visit_RF40pct_res) #18/34
#QRILC
visit_QRILC5pct_res <- compare_visits_wilcoxon(QRILC_5pct)
visit_QRILC5pct_signif <- significance(visit_QRILC5pct_res) #16/34
visit_QRILC10pct_res <- compare_visits_wilcoxon(QRILC_10pct)
visit_QRILC10pct_signif <- significance(visit_QRILC10pct_res) #0/34
visit_QRILC20pct_res <- compare_visits_wilcoxon(QRILC_20pct)
visit_QRILC20pct_signif <- significance(visit_QRILC20pct_res) #0/34
visit_QRILC25pct_res <- compare_visits_wilcoxon(QRILC_25pct)
visit_QRILC25pct_signif <- significance(visit_QRILC25pct_res) #0/34
visit_QRILC30pct_res <- compare_visits_wilcoxon(QRILC_30pct)
visit_QRILC30pct_signif <- significance(visit_QRILC30pct_res) #0/34
visit_QRILC40pct_res <- compare_visits_wilcoxon(QRILC_40pct)
visit_QRILC40pct_signif <- significance(visit_QRILC40pct_res) #0/34

#store results in dataframe
df_significant_visit_tests <- data.frame(
  Method = rep(c("Halfmin", "KNN", "RF", "QRILC"), each = 6),  
  Missing_Percentage = rep(c("5%", "10%", "20%", "25%", "30%", "40%"), times = 4), 
  Signif_Wilcoxon_Visit = c(
    nrow(visit_Halfmin5pct_signif),  
    nrow(visit_Halfmin10pct_signif),  
    nrow(visit_Halfmin20pct_signif),  
    nrow(visit_Halfmin25pct_signif),  
    nrow(visit_Halfmin30pct_signif),  
    nrow(visit_Halfmin40pct_signif),  
    nrow(visit_KNN5pct_signif),  
    nrow(visit_KNN10pct_signif),  
    nrow(visit_KNN20pct_signif),  
    nrow(visit_KNN25pct_signif),  
    nrow(visit_KNN30pct_signif),  
    nrow(visit_KNN40pct_signif),  
    nrow(visit_RF5pct_signif),  
    nrow(visit_RF10pct_signif),  
    nrow(visit_RF20pct_signif),  
    nrow(visit_RF25pct_signif),  
    nrow(visit_RF30pct_signif),  
    nrow(visit_RF40pct_signif),  
    nrow(visit_QRILC5pct_signif),  
    nrow(visit_QRILC10pct_signif),  
    nrow(visit_QRILC20pct_signif),  
    nrow(visit_QRILC25pct_signif),  
    nrow(visit_QRILC30pct_signif),  
    nrow(visit_QRILC40pct_signif)
  )
)

#save dataframe 
write_csv(df_significant_visit_tests, "/Users/marcinebessire/Desktop/project/Result_Wilcoxon_Visits_FAO.csv")

# ------------------------------------
# Part 2.5: Check Significance of Tests
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
signif_halfmin_10pct_t_test <- significance(t_test_half_min_10pct) #0/34
signif_halfmin_20pct_t_test <- significance(t_test_half_min_20pct) #0/34
signif_halfmin_25pct_t_test <- significance(t_test_half_min_25pct) #0/34
signif_halfmin_30pct_t_test <- significance(t_test_half_min_30pct) #0/34
signif_halfmin_40pct_t_test <- significance(t_test_half_min_40pct) #14/34
#KNN
signif_KNN_5pct_t_test <- significance(t_test_KNN_5pct) #0/34
signif_KNN_10pct_t_test <- significance(t_test_KNN_10pct) #0/34
signif_KNN_20pct_t_test <- significance(t_test_KNN_20pct) #0/34
signif_KNN_25pct_t_test <- significance(t_test_KNN_25pct) #0/34
signif_KNN_30pct_t_test <- significance(t_test_KNN_30pct) #0/34
signif_KNN_40pct_t_test <- significance(t_test_KNN_40pct) #7/34
#RF
signif_RF_5pct_t_test <- significance(t_test_RF_5pct) #0/34
signif_RF_10pct_t_test <- significance(t_test_RF_10pct) #0/34
signif_RF_20pct_t_test <- significance(t_test_RF_20pct) #0/34
signif_RF_25pct_t_test <- significance(t_test_RF_25pct) #0/34
signif_RF_30pct_t_test <- significance(t_test_RF_30pct) #0/34
signif_RF_40pct_t_test <- significance(t_test_RF_40pct) #0/34
#QRILC
signif_QRILC_5pct_t_test <- significance(t_test_QRILC_5pct) #0/34
signif_QRILC_10pct_t_test <- significance(t_test_QRILC_10pct) #0/34
signif_QRILC_20pct_t_test <- significance(t_test_QRILC_20pct) #0/34
signif_QRILC_25pct_t_test <- significance(t_test_QRILC_25pct) #0/34
signif_QRILC_30pct_t_test <- significance(t_test_QRILC_30pct) #0/34
signif_QRILC_40pct_t_test <- significance(t_test_QRILC_40pct) #0/34


#call significance function on wilcoxon results
#Halfmin
signif_halfmin_5pct_wilcox <- significance(wilcox_half_min_5pct) #0/34
signif_halfmin_10pct_wilcox <- significance(wilcox_half_min_10pct) #0/34
signif_halfmin_20pct_wilcox <- significance(wilcox_half_min_20pct) #0/34
signif_halfmin_25pct_wilcox <- significance(wilcox_half_min_25pct) #0/34
signif_halfmin_30pct_wilcox <- significance(wilcox_half_min_30pct) #0/34
signif_halfmin_40pct_wilcox <- significance(wilcox_half_min_40pct) #21/34
#KNN
signif_KNN_5pct_wilcox <- significance(wilcox_half_min_5pct) #0/34
signif_KNN_10pct_wilcox <- significance(wilcox_half_min_10pct) #0/34
signif_KNN_20pct_wilcox <- significance(wilcox_half_min_20pct) #0/34
signif_KNN_25pct_wilcox <- significance(wilcox_half_min_25pct) #0/34
signif_KNN_30pct_wilcox <- significance(wilcox_half_min_30pct) #0/34
signif_KNN_40pct_wilcox <- significance(wilcox_half_min_40pct) #21/34
#RF
signif_RF_5pct_wilcox <- significance(wilcox_RF_5pct) #0/34
signif_RF_10pct_wilcox <- significance(wilcox_RF_10pct) #0/34
signif_RF_20pct_wilcox <- significance(wilcox_RF_20pct) #0/34
signif_RF_25pct_wilcox <- significance(wilcox_RF_25pct) #0/34
signif_RF_30pct_wilcox <- significance(wilcox_RF_30pct) #0/34
signif_RF_40pct_wilcox <- significance(wilcox_RF_40pct) #0/34
#QRILC
signif_QRILC_5pct_wilcox <- significance(wilcox_QRILC_5pct) #0/34
signif_QRILC_10pct_wilcox <- significance(wilcox_QRILC_10pct) #0/34
signif_QRILC_20pct_wilcox <- significance(wilcox_QRILC_20pct) #0/34
signif_QRILC_25pct_wilcox <- significance(wilcox_QRILC_25pct) #0/34
signif_QRILC_30pct_wilcox <- significance(wilcox_QRILC_30pct) #0/34
signif_QRILC_40pct_wilcox <- significance(wilcox_QRILC_40pct) #5/34


#Store significant results
df_significant_tests <- data.frame(
  Method = rep(c("Halfmin", "KNN", "RF", "QRILC"), each = 6),  
  Missing_Percentage = rep(c("5%", "10%", "20%", "25%", "30%", "40%"), times = 4),  
  Signif_T_Test = c(
    nrow(signif_halfmin_5pct_t_test),  
    nrow(signif_halfmin_10pct_t_test),  
    nrow(signif_halfmin_20pct_t_test),  
    nrow(signif_halfmin_25pct_t_test),  
    nrow(signif_halfmin_30pct_t_test),  
    nrow(signif_halfmin_40pct_t_test),  
    nrow(signif_KNN_5pct_t_test),  
    nrow(signif_KNN_10pct_t_test),  
    nrow(signif_KNN_20pct_t_test),  
    nrow(signif_KNN_25pct_t_test),  
    nrow(signif_KNN_30pct_t_test),  
    nrow(signif_KNN_40pct_t_test),  
    nrow(signif_RF_5pct_t_test),  
    nrow(signif_RF_10pct_t_test),  
    nrow(signif_RF_20pct_t_test),  
    nrow(signif_RF_25pct_t_test),  
    nrow(signif_RF_30pct_t_test),  
    nrow(signif_RF_40pct_t_test),  
    nrow(signif_QRILC_5pct_t_test),  
    nrow(signif_QRILC_10pct_t_test),  
    nrow(signif_QRILC_20pct_t_test),  
    nrow(signif_QRILC_25pct_t_test),  
    nrow(signif_QRILC_30pct_t_test),  
    nrow(signif_QRILC_40pct_t_test)
  ),
  Signif_Wilcoxon = c(
    nrow(signif_halfmin_5pct_wilcox),  
    nrow(signif_halfmin_10pct_wilcox),  
    nrow(signif_halfmin_20pct_wilcox),  
    nrow(signif_halfmin_25pct_wilcox),  
    nrow(signif_halfmin_30pct_wilcox),  
    nrow(signif_halfmin_40pct_wilcox),  
    nrow(signif_KNN_5pct_wilcox),  
    nrow(signif_KNN_10pct_wilcox),  
    nrow(signif_KNN_20pct_wilcox),  
    nrow(signif_KNN_25pct_wilcox),  
    nrow(signif_KNN_30pct_wilcox),  
    nrow(signif_KNN_40pct_wilcox),  
    nrow(signif_RF_5pct_wilcox),  
    nrow(signif_RF_10pct_wilcox),  
    nrow(signif_RF_20pct_wilcox),  
    nrow(signif_RF_25pct_wilcox),  
    nrow(signif_RF_30pct_wilcox),  
    nrow(signif_RF_40pct_wilcox),  
    nrow(signif_QRILC_5pct_wilcox),  
    nrow(signif_QRILC_10pct_wilcox),  
    nrow(signif_QRILC_20pct_wilcox),  
    nrow(signif_QRILC_25pct_wilcox),  
    nrow(signif_QRILC_30pct_wilcox),  
    nrow(signif_QRILC_40pct_wilcox)
  )
)

write_csv(df_significant_tests, "/Users/marcinebessire/Desktop/project/Result_StatTests_FAO.csv")


# ------------------------------------
# Part 3: NRMSE
# ------------------------------------

#function to calculate NRMSE
calculate_nrmse <- function(original, imputed) {
  #take numeric columns 
  numeric_col_names <- names(original)[6:ncol(original)]
  
  #calcualte for each column nrmse
  nrmse_values <- sapply(numeric_col_names, function(col){
    actual_val <- original[[col]]
    predicted_val <- imputed[[col]]
    
    #ensure no missing values for comparison
    valid_indices <- !is.na(actual_val) & !is.na(predicted_val)
    
    if (sum(valid_indices) > 2) {  #if enouhg datapoints
      mse <- mean((actual_val[valid_indices] - predicted_val[valid_indices])^2) #mean squared error
      rmse <- sqrt(mse) #root mean squared error
      norm_factor <- max(actual_val, na.rm = TRUE) - min(actual_val, na.rm = TRUE) #denominator 
      
      if (norm_factor > 0) {
        return(rmse/norm_factor) #NRMSE
      } else {
        return(NA) #avoid dividing by 0
      }
    } else {
      return(NA)  #return NA else (if not > 2)
    }
  })

  return(data.frame(Variable = numeric_col_names, NRMSE = nrmse_values))
  
}

#call function to compute NRMSE
#Halfmin
nrmse_res_halfmin_5pct <- calculate_nrmse(FAO_original, Halfmin_5pct)
nrmse_res_halfmin_10pct <- calculate_nrmse(FAO_original, Halfmin_10pct)
nrmse_res_halfmin_20pct <- calculate_nrmse(FAO_original, Halfmin_20pct)
nrmse_res_halfmin_25pct <- calculate_nrmse(FAO_original, Halfmin_25pct)
nrmse_res_halfmin_30pct <- calculate_nrmse(FAO_original, Halfmin_30pct)
nrmse_res_halfmin_40pct <- calculate_nrmse(FAO_original, Halfmin_40pct)
#KNN
nrmse_res_KNN_5pct <- calculate_nrmse(FAO_original, KNN_5pct)
nrmse_res_KNN_10pct <- calculate_nrmse(FAO_original, KNN_10pct)
nrmse_res_KNN_20pct <- calculate_nrmse(FAO_original, KNN_20pct)
nrmse_res_KNN_25pct <- calculate_nrmse(FAO_original, KNN_25pct)
nrmse_res_KNN_30pct <- calculate_nrmse(FAO_original, KNN_30pct)
nrmse_res_KNN_40pct <- calculate_nrmse(FAO_original, KNN_40pct)
#RF
nrmse_res_RF_5pct <- calculate_nrmse(FAO_original, RF_5pct)
nrmse_res_RF_10pct <- calculate_nrmse(FAO_original, RF_10pct)
nrmse_res_RF_20pct <- calculate_nrmse(FAO_original, RF_20pct)
nrmse_res_RF_25pct <- calculate_nrmse(FAO_original, RF_25pct)
nrmse_res_RF_30pct <- calculate_nrmse(FAO_original, RF_30pct)
nrmse_res_RF_40pct <- calculate_nrmse(FAO_original, RF_40pct)
#QRILC
nrmse_res_QRILC_5pct <- calculate_nrmse(FAO_original, QRILC_5pct)
nrmse_res_QRILC_10pct <- calculate_nrmse(FAO_original, QRILC_10pct)
nrmse_res_QRILC_20pct <- calculate_nrmse(FAO_original, QRILC_20pct)
nrmse_res_QRILC_25pct <- calculate_nrmse(FAO_original, QRILC_25pct)
nrmse_res_QRILC_30pct <- calculate_nrmse(FAO_original, QRILC_30pct)
nrmse_res_QRILC_40pct <- calculate_nrmse(FAO_original, QRILC_40pct)

# ------------------------------------
# Part 3.1: NRMSE Plot 
# ------------------------------------



# ------------------------------------
# Part 4: Normlaized Mean Difference
# ------------------------------------

#function to calcualte normalized mean difference
norm_mean_diff <- function(original, imputed, method, percentage) {
  #numeric columns
  original_numeric <- original[, 6:ncol(original)]
  imputed_numeric <- imputed[, 6:ncol(imputed)]
  
  #mean before imputation
  mean_before <- original_numeric %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_Before")
  
  #means after imputation
  mean_after <- imputed_numeric %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_After")
  
  #merge before and after mean values
  mean_comparison <- left_join(mean_before, mean_after, by = "Metabolite")
  
  #compute normalized difference: (Mean_After - Mean_Before) / Mean_Before
  mean_comparison <- mean_comparison %>%
    mutate(Normalized_Difference = (Mean_After - Mean_Before) / Mean_Before)
  
  #title 
  plot_title <- paste0("Normalized Difference with ", percentage, "% Missing Values using ", method, " Imputation")
  
  #plot the density of the normalized difference
  plot <- ggplot(mean_comparison, aes(x = Normalized_Difference)) +
    geom_density(fill = "blue", alpha = 0.4, color = "black") + 
    theme_minimal() +
    labs(title = plot_title,  
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.2, 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")  
  
  print(plot)  
  return(mean_comparison)  
}

#call normalized difference function
#Halfmin
norm_diff_Halfmin_5pct <- norm_mean_diff(FAO_original, Halfmin_5pct, "Half-min", 5)
norm_diff_Halfmin_30pct <- norm_mean_diff(FAO_original, Halfmin_30pct, "Half-min", 30)
norm_diff_Halfmin_40pct <- norm_mean_diff(FAO_original, Halfmin_40pct, "Half-min", 40)
#KNN
norm_diff_KNN_5pct <- norm_mean_diff(FAO_original, KNN_5pct, "KNN", 5)
norm_diff_KNN_40pct <- norm_mean_diff(FAO_original, KNN_40pct, "KNN", 40)
#RF
norm_diff_RF_5pct <- norm_mean_diff(FAO_original, RF_5pct, "RF", 5)
norm_diff_RF_40pct <- norm_mean_diff(FAO_original, RF_40pct, "RF", 40)
#QRILC
norm_diff_QRILC_5pct <- norm_mean_diff(FAO_original, QRILC_5pct, "QRILC", 5)
norm_diff_QRILC_40pct <- norm_mean_diff(FAO_original, QRILC_40pct, "QRILC", 40)


# ------------------------------------
# Part 5: Plot Distirbution
# ------------------------------------

#function to plot the distirbution before and after imputation
distirbution_plot <- function(original, imputed, method, percentage){
  #numeric columns
  original_numeric <- original[, 6:ncol(original)]
  imputed_numeric <- imputed[, 6:ncol(imputed)]
  
  #convert data to long format for visualization
  imputed_long <- imputed_numeric %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Imputed_Data")
  
  original_long <- original_numeric %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Original_Data")
  
  #identify imputed values => missing values were replaced by half of the minimum observed value.
  imputed_only <- original_long %>%
    inner_join(imputed_long, by = "Metabolite") %>%
    filter(Original_Data != Imputed_Data) %>%
    mutate(Dataset = "Imputed_Values")
  
  #merge original and imputed datasets
  comparison <- original_long %>%
    left_join(imputed_long %>%
                distinct(Metabolite, .keep_all = TRUE), by = "Metabolite") %>%
    pivot_longer(cols = c("Original_Data", "Imputed_Data"), 
                 names_to = "Dataset", values_to = "Value")
  
  #add Imputed_Only as a separate dataset
  imputed_only <- imputed_only %>%
    mutate(Value = Imputed_Data, Dataset = "Imputed_Values") %>%
    select(Metabolite, Dataset, Value)
  
  #combine both datasets
  comparison <- bind_rows(comparison, imputed_only)
  
  #title 
  plot_title <- paste0("Distirbution of Original Data and Data with ", percentage, "% Missing Values using ", method, " Imputation")

  #plot 
  plot <- ggplot(comparison, aes(x = Value, fill = Dataset)) +
    geom_density(alpha = 0.5) +  #transparency for overlapping
    labs(title = plot_title,
         x = "Metabolite Value",
         y = "Density") +
    theme_minimal() + 
    scale_fill_manual(values = c("Original_Data" = "lightblue", 
                                 "Imputed_Data" = "red", 
                                 "Imputed_Values" = "lightyellow")) +
    xlim(-100,500)
    
    print(plot)
    return(comparison)
} 

#call function for Distirbution plot
#Halfmin
dist_Halfmin_5pct <- distirbution_plot(FAO_original, Halfmin_5pct, "Half-min", 5)
dist_Halfmin_40pct <- distirbution_plot(FAO_original, Halfmin_40pct, "Half-min", 40)
#KNN
dist_KNN_5pct <- distirbution_plot(FAO_original, KNN_5pct, "KNN", 5)
dist_KNN_40pct <- distirbution_plot(FAO_original, KNN_40pct, "KNN", 40)
#RF
dist_RF_5pct <- distirbution_plot(FAO_original, RF_5pct, "RF", 5)
dist_RF_40pct <- distirbution_plot(FAO_original, RF_40pct, "RF", 40)
#QRILC
dist_QRILC_5pct <- distirbution_plot(FAO_original, QRILC_5pct, "QRILC", 5)
dist_QRILC_40pct <- distirbution_plot(FAO_original, QRILC_40pct, "QRILC", 40)




