#load necessary libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(missForest)
library(imputeLCMD)
library(reshape2) #for melt
library(FSA) #for Dunns test


#load original data
BAS_original <- read.csv("/Users/marcinebessire/Desktop/project/BAS_data.csv", check.names = FALSE)

# ------------------------------------
# Mean, Median and CV of FAO original
# ------------------------------------

#numeric data
BAS_numeric <- BAS_original[, 6:ncol(BAS_original)]

#initialize empty dataset
summary_stats <- data.frame(Metabolite = colnames(BAS_numeric),
                            Mean = NA,
                            Median = NA, 
                            CV = NA)

#loop through each column and compute statistics
for (i in seq_along(BAS_numeric)) {
  col_data <- BAS_numeric[[i]]
  
  #compte mean, median and CV 
  mean_val <- mean(col_data, na.rm = TRUE)
  median_val <- median(col_data, na.rm = TRUE)
  cv_val <- (sd(col_data, na.rm = TRUE) / mean_val) * 100 #in percent
  
  #store in dataframe
  summary_stats[i, "Mean"] <- mean_val
  summary_stats[i, "Median"] <- median_val
  summary_stats[i, "CV"] <- cv_val
}

#change to long format for plotting 
summary_long <- melt(summary_stats, id.vars = "Metabolite", measure.vars = c("Mean", "Median"))

#bar plot mean median
ggplot(summary_long, aes(x = reorder(Metabolite, value), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Mean and Median of Metabolites",
       x = "Metabolite", y = "Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("steelblue", "orange"))

#bar plot CV per metabolite
ggplot(summary_stats, aes(x = reorder(Metabolite, -CV), y = CV)) +
  geom_bar(stat = "identity", fill = "magenta", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Coefficient of Variation (CV) per Metabolite",
       x = "Metabolite", y = "CV (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ---------------------------------------------------------------
# Remove values with CV > 20% and then remove col with MV > 20%
# ---------------------------------------------------------------

#function to caluclate cv for each numeric column in each file
calculate_cv <- function(data) {
  #exclude non-data columns
  df_data <- data[, 1:5]
  
  #keep only numeric columns
  df_numeric <- data[, 6:ncol(data)]
  
  #set inf values to NA (because are like missing values)
  df_numeric[df_numeric == Inf | df_numeric == -Inf] <- NA 
  
  #calculate coefficient of variance
  cv_values <- sapply(df_numeric, function(x) {
    if(all(is.na(x))) {
      return(NA_real_)
    }
    (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
  })
  
  #convert results into data frame
  cv_df <- data.frame(Column = names(cv_values), CV_Percentage = cv_values)
  
  #rename column
  colnames(cv_df)[2] <- "CV [%]"
  
  #return CV dataframe
  return(cv_df)
}

#calculate cv for BAS data
cv_data_BAS <- calculate_cv(BAS_original) #keep all because most of them are above 20% 

# ------------------------------
# Remove columns with MV > 20%
# ------------------------------

#remove col with >20% MV
remove_mv <- function(df, threshold = 0.2) {
  #compute the number of missing values per column
  na_counts <- colSums(is.na(df))
  
  #compute the threshold count (20% of total rows)
  na_threshold <- threshold * nrow(df)
  
  #identify columns to keep (those with <= 20% missing values)
  cols_to_keep <- names(df)[na_counts <= na_threshold]
  
  #return the filtered dataframe
  return(df[, cols_to_keep, drop = FALSE])
}

#call function to remove based on MV 
BAS_final <- remove_mv(BAS_original) #17/18 Metabolites remained 

# ------------------------------------
# Part 0: Visit 1 vs Visit 2 Plot 
# ------------------------------------

#before removing mv
#reshape data into long format (use melt to convert dataframe from wide into logn format)
BAS_long_before <- melt(BAS_original, id.vars = c("ID", "Participant", "MonthDay", "Year", "Visit"),
                 variable.name = "Metabolite", value.name = "Value")

#generate boxplot for each metabolite 
ggplot(BAS_long_before, aes(x = Visit, y = Value, fill = Visit)) +
  geom_boxplot() +
  facet_wrap(~Metabolite, scales = "free") +
  theme_minimal() +
  labs(title = "Comparison of Metabolites Between Visits", 
       x = "Visit",
       y = "Metabolite Measurement") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

#after removing one column
#reshape data into long format (use melt to convert dataframe from wide into logn format)
BAS_long_after <- melt(BAS_final, id.vars = c("ID", "Participant", "MonthDay", "Year", "Visit"),
                        variable.name = "Metabolite", value.name = "Value")

#generate boxplot for each metabolite 
ggplot(BAS_long_after, aes(x = Visit, y = Value, fill = Visit)) +
  geom_boxplot() +
  facet_wrap(~Metabolite, scales = "free") +
  theme_minimal() +
  labs(title = "Comparison of Metabolites Between Visits", 
       x = "Visit",
       y = "Metabolite Measurement") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))


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
  numeric_data <- data_imputed[, 6:ncol(data_imputed)]
  
  #loop through column
  for (col in names(numeric_data)) { 
    min_val <- min(numeric_data[[col]], na.rm = TRUE) #find min value 
    numeric_data[[col]][is.na(numeric_data[[col]])] <- 0.5 * min_val #calculate half min and replace NA
  }
  
  final_df <- cbind(meta_data, numeric_data)
  
  return(final_df)
}

#call function for Half-min imputation
Halfmin_BAS <- half_min_imputation(BAS_final)

# ------------------------------------
# Part 1.2: KNN Imputation
# ------------------------------------

KNN_imputation <- function(data) {
  #make copy 
  data_imputed <- data
  
  #select only numeric data
  numeric_data <- data_imputed[, 6:ncol(data_imputed)]
  
  #metadata
  meta_data <- data_imputed[, 1:5]
  
  #transform first into matrix
  imputed_data <- impute.knn(as.matrix(t(numeric_data)), rowmax = 0.5, colmax = 1) 
  
  #transform back into dataframe
  imputed_data <- as.data.frame(t(imputed_data$data))
  
  final_df <- cbind(meta_data, imputed_data)
  
  return(final_df)
}

#call function for KNN
KNN_BAS <- KNN_imputation(BAS_final)

# ------------------------------------
# Part 1.3: RF Imputation
# ------------------------------------

#make function for RF
RF_imputation <- function(data) {
  #make a copy 
  data_copy <- data 
  
  #select only numeric data
  numeric_data <- data_copy[, 6:ncol(data_copy)]
  
  #metadata
  meta_data <- data_copy[, 1:5]
  
  #apply miss forest on data to impute MV
  imputed_data <- missForest(numeric_data, maxiter = 10, ntree = 100)
  imputed_df <- as.data.frame(imputed_data$ximp)
  
  #return final dataframe
  final_df <- cbind(meta_data, imputed_df)
}

#call function for RF
RF_BAS <- RF_imputation(BAS_final)

# ------------------------------------
# Part 1.4: QRILC Imputation
# ------------------------------------

#QRILC imputation 
QRILC_impuation <- function(data) {
  #copy data
  data_copy <- data
  
  #select only numeric data
  numeric_data <- data_copy[, 6:ncol(data_copy)]
  
  #metadata
  meta_data <- data_copy[, 1:5]
  
  #transfer to log data
  log_data <- log(numeric_data + 1e-6)  #add small constant to avoid log(0)
  
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
QRILC_BAS <- QRILC_impuation(BAS_final)


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
  numeric_data <- data_copy[, 6:ncol(data_copy)]
  
  #shapiro test
  shapiro_results <- apply(numeric_data, 2, function(x) shapiro.test(x)$p.value)
  
  #adjust p-value 
  adjusted_p_values <- p.adjust(shapiro_results, method = "BH")
  
  #concert to dataframe
  shapiro_df <- data.frame(
    Metabolite = names(shapiro_results), 
    p_value = shapiro_results,
    Adjusted_P_Value = adjusted_p_values
  ) 
  
  #print how many are significant = if p-value < 0.05 then not normal distribution
  significant_count <- sum(shapiro_df$Adjusted_P_Value < 0.05, na.rm = TRUE)
  cat("\n-------------------------------\n")
  cat("Number of non-normally distributed metabolites:", significant_count, "\n")
  cat("-------------------------------\n")
  
  
  #return dataframe
  return(shapiro_df)
}

#call shapiro function
#original
shapiro_original <- shapiro_test(BAS_original) #15/34 are significant
shapiro_final <- shapiro_test(BAS_final) #14/34 are significant
#halfmin
shapiro_halfmin <- shapiro_test(Halfmin_BAS) #14/34 are significant
#knn
shapiro_KNN <- shapiro_test(KNN_BAS) #14/34 are significant
#rf
shapiro_RF <- shapiro_test(RF_BAS) #14/34 are significant
#QRILC
shapiro_QRILC <- shapiro_test(QRILC_BAS) #14/34 are significant

#store the results in a data frame
shapiro_results_df <- data.frame(
  Method = c("Original", "Final", "Halfmin", "KNN", "RF", "QRILC"),
  Non_Normal_Count = c(
    sum(shapiro_original$Adjusted_P_Value < 0.05, na.rm = TRUE),
    sum(shapiro_final$Adjusted_P_Value < 0.05, na.rm = TRUE),
    sum(shapiro_halfmin$Adjusted_P_Value < 0.05, na.rm = TRUE),
    sum(shapiro_KNN$Adjusted_P_Value < 0.05, na.rm = TRUE),
    sum(shapiro_RF$Adjusted_P_Value < 0.05, na.rm = TRUE),
    sum(shapiro_QRILC$Adjusted_P_Value < 0.05, na.rm = TRUE)
  )
)

#bar Plot
ggplot(shapiro_results_df, aes(x = Method, y = Non_Normal_Count, fill = Method)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Non-Normal Metabolites", 
       x = "Method", 
       y = "Count of Non-Normal Metabolites") +
  theme_minimal()


#CONTINUE HERE -----------------


# ------------------------------------
# Part 2.2: T-test
# ------------------------------------

t_test_func <- function(original, imputed) {
  #numeric columns
  numeric_data <- original[, 6:ncol(original)]
  
  #metabolte columns
  metabolite_cols <- colnames(numeric_data)
  
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
  
  #count significant results
  significant_count <- sum(results$adj_p_value < 0.05)
  
  # Print summary
  cat("\n-------------------------------\n")
  cat("Number of significantly different metabolites:", significant_count, "\n")
  cat("-------------------------------\n")
  
  
  return(results)
  
}

#call t-test function
#Half-min
t_test_half_min_5pct <- t_test_func(FAO_original, Halfmin_5pct) #0/34
t_test_half_min_10pct <- t_test_func(FAO_original, Halfmin_10pct) #0/34
t_test_half_min_20pct <- t_test_func(FAO_original, Halfmin_20pct) #0/34
t_test_half_min_25pct <- t_test_func(FAO_original, Halfmin_25pct) #0/34
t_test_half_min_30pct <- t_test_func(FAO_original, Halfmin_30pct) #0/34
t_test_half_min_40pct <- t_test_func(FAO_original, Halfmin_40pct) #13/34
#KNN
t_test_KNN_5pct <- t_test_func(FAO_original, KNN_5pct) #0/34
t_test_KNN_10pct <- t_test_func(FAO_original, KNN_10pct) #0/34
t_test_KNN_20pct <- t_test_func(FAO_original, KNN_20pct) #0/34
t_test_KNN_25pct <- t_test_func(FAO_original, KNN_25pct) #0/34
t_test_KNN_30pct <- t_test_func(FAO_original, KNN_30pct) #0/34
t_test_KNN_40pct <- t_test_func(FAO_original, KNN_40pct) #7/34
#RF
t_test_RF_5pct <- t_test_func(FAO_original, RF_5pct) #0/34
t_test_RF_10pct <- t_test_func(FAO_original, RF_10pct) #0/34
t_test_RF_20pct <- t_test_func(FAO_original, RF_20pct) #0/34
t_test_RF_25pct <- t_test_func(FAO_original, RF_25pct) #0/34
t_test_RF_30pct <- t_test_func(FAO_original, RF_30pct) #0/34
t_test_RF_40pct <- t_test_func(FAO_original, RF_40pct) #0/34
#QRILC
t_test_QRILC_5pct <- t_test_func(FAO_original, QRILC_5pct) #0/34
t_test_QRILC_10pct <- t_test_func(FAO_original, QRILC_10pct) #0/34
t_test_QRILC_20pct <- t_test_func(FAO_original, QRILC_20pct) #0/34
t_test_QRILC_25pct <- t_test_func(FAO_original, QRILC_25pct) #0/34
t_test_QRILC_30pct <- t_test_func(FAO_original, QRILC_30pct) #0/34
t_test_QRILC_40pct <- t_test_func(FAO_original, QRILC_40pct) #0/34


# ------------------------------------
# Part 2.3: Wilcoxon rank-sum Test
# ------------------------------------

wilcoxon_func <- function(original, imputed) {
  #numeric columns
  numeric_data <- original[, 6:ncol(original)]
  
  #metabolte columns
  metabolite_cols <- colnames(numeric_data)
  
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
  
  #count significant results
  significant_count <- sum(results$adj_p_value < 0.05)
  
  #print summary
  cat("\n-------------------------------\n")
  cat("Number of significantly different metabolites (Wilcoxon test):", significant_count, "\n")
  cat("-------------------------------\n")
  
  
  return(results)
}

#call wilcoxon function 
#Halfmin
wilcox_half_min_5pct <- wilcoxon_func(FAO_original, Halfmin_5pct) #0/34
wilcox_half_min_10pct <- wilcoxon_func(FAO_original, Halfmin_10pct) #0/34
wilcox_half_min_20pct <- wilcoxon_func(FAO_original, Halfmin_20pct) #0/34
wilcox_half_min_25pct <- wilcoxon_func(FAO_original, Halfmin_25pct) #0/34
wilcox_half_min_30pct <- wilcoxon_func(FAO_original, Halfmin_30pct) #0/34
wilcox_half_min_40pct <- wilcoxon_func(FAO_original, Halfmin_40pct) #15/34
#KNN
wilcox_KNN_5pct <- wilcoxon_func(FAO_original, KNN_5pct) #0/34
wilcox_KNN_10pct <- wilcoxon_func(FAO_original, KNN_10pct) #0/34
wilcox_KNN_20pct <- wilcoxon_func(FAO_original, KNN_20pct) #0/34
wilcox_KNN_25pct <- wilcoxon_func(FAO_original, KNN_25pct) #0/34
wilcox_KNN_30pct <- wilcoxon_func(FAO_original, KNN_30pct) #0/34
wilcox_KNN_40pct <- wilcoxon_func(FAO_original, KNN_40pct) #0/34
#RF
wilcox_RF_5pct <- wilcoxon_func(FAO_original, RF_5pct) #0/34
wilcox_RF_10pct <- wilcoxon_func(FAO_original, RF_10pct) #0/34
wilcox_RF_20pct <- wilcoxon_func(FAO_original, RF_20pct) #0/34
wilcox_RF_25pct <- wilcoxon_func(FAO_original, RF_25pct) #0/34
wilcox_RF_30pct <- wilcoxon_func(FAO_original, RF_30pct) #0/34
wilcox_RF_40pct <- wilcoxon_func(FAO_original, RF_40pct) #0/34
#QRILC
wilcox_QRILC_5pct <- wilcoxon_func(FAO_original, QRILC_5pct) #0/34
wilcox_QRILC_10pct <- wilcoxon_func(FAO_original, QRILC_10pct) #0/34
wilcox_QRILC_20pct <- wilcoxon_func(FAO_original, QRILC_20pct) #0/34
wilcox_QRILC_25pct <- wilcoxon_func(FAO_original, QRILC_25pct) #0/34
wilcox_QRILC_30pct <- wilcoxon_func(FAO_original, QRILC_30pct) #0/34
wilcox_QRILC_40pct <- wilcoxon_func(FAO_original, QRILC_40pct) #0/34

# ------------------------------------
# Part 2.4: Wilcoxon and T-test for Visit 1 vs 2
# ------------------------------------

#comapte visit 1 and visit 2 for each metabolite
visit_statistical_tests <- function(data) {
  #numeric metabolite columns (from column 6 onward)
  metabolites <- colnames(data)[6:ncol(data)]
  
  #separate Visit 1 and Visit 2
  visit1 <- data %>% filter(Visit == "Visit 1") %>% select(Participant, all_of(metabolites))
  visit2 <- data %>% filter(Visit == "Visit 2") %>% select(Participant, all_of(metabolites))
  
  #merge Visit 1 and Visit 2 data by participant
  paired_data <- merge(visit1, visit2, by = "Participant", suffixes = c("_V1", "_V2"))
  
  #convert all columns to numeric 
  data_before <- as.data.frame(lapply(paired_data[, paste0(metabolites, "_V1")], as.numeric))
  data_after <- as.data.frame(lapply(paired_data[, paste0(metabolites, "_V2")], as.numeric))
  
  #vectors to store p-values
  p_values_wilcox <- numeric(length(metabolites))
  p_values_ttest <- numeric(length(metabolites))
  
  #perform tests for each metabolite
  for (i in seq_along(metabolites)) {
    if (any(is.na(data_before[, i])) || any(is.na(data_after[, i]))) {
      warning(paste("Non-numeric values detected in column", metabolites[i], ". Skipping this column."))
      next
    }
    
    #wilcoxon signed-rank test
    test_result_wilcox <- wilcox.test(data_before[, i], data_after[, i], paired = TRUE, exact = FALSE)
    p_values_wilcox[i] <- test_result_wilcox$p.value
    
    #paired t-test
    test_result_ttest <- t.test(data_before[, i], data_after[, i], paired = TRUE)
    p_values_ttest[i] <- test_result_ttest$p.value
  }
  
  #adjust p-values using Benjamini-Hochberg (FDR) correction
  adjusted_p_values_wilcox <- p.adjust(p_values_wilcox, method = "BH")
  adjusted_p_values_ttest <- p.adjust(p_values_ttest, method = "BH")
  
  #store results in a dataframe
  results <- data.frame(
    Metabolite = metabolites,
    Wilcoxon_P_Value = p_values_wilcox,
    Wilcoxon_Adjusted_P = adjusted_p_values_wilcox,
    TTest_P_Value = p_values_ttest,
    TTest_Adjusted_P = adjusted_p_values_ttest
  )
  
  #define significance threshold
  alpha <- 0.05
  
  #count significant metabolites
  significant_wilcox <- sum(results$Wilcoxon_Adjusted_P < alpha, na.rm = TRUE)
  significant_ttest <- sum(results$TTest_Adjusted_P < alpha, na.rm = TRUE)
  
  #print results properly
  cat("\n-------------------------------\n")
  cat("Number of significant metabolites:\n")
  cat("- Wilcoxon test:", significant_wilcox, "\n")
  cat("- Paired t-test:", significant_ttest, "\n")
  cat("-------------------------------\n")
  
  
  #return results dataframe
  return(results)
}

#call function to compare Visit 1 vs. Visit 2
#original 
visit_original_res <- visit_statistical_tests(FAO_original)  #18/34 W and 16/34 T
#Halfmin
visit_Halfmin5pct_res <- visit_statistical_tests(Halfmin_5pct) #13/34 W and 15/34 T
visit_Halfmin10pct_res <- visit_statistical_tests(Halfmin_10pct) #0/34 W and 8/34 T
visit_Halfmin20pct_res <- visit_statistical_tests(Halfmin_20pct) #0/34 W and 6/34 T
visit_Halfmin25pct_res <- visit_statistical_tests(Halfmin_25pct) #0/34 W and 5/34 T
visit_Halfmin30pct_res <- visit_statistical_tests(Halfmin_30pct) #0/34 W and 0/34 T
visit_Halfmin40pct_res <- visit_statistical_tests(Halfmin_40pct) #0/34 W and 0/34 T
#KNN
visit_KNN5pct_res <- visit_statistical_tests(KNN_5pct) #14/34 W 15/34 T
visit_KNN10pct_res <- visit_statistical_tests(KNN_10pct) #11/34 W 9/34 T
visit_KNN20pct_res <- visit_statistical_tests(KNN_20pct) #10/34 W 10/34 T
visit_KNN25pct_res <- visit_statistical_tests(KNN_25pct) #8/34 W 10/34 T
visit_KNN30pct_res <- visit_statistical_tests(KNN_30pct) #8/34 W 7/34 T
visit_KNN40pct_res <- visit_statistical_tests(KNN_40pct) #0/34 W 10/34 T
#RF
visit_RF5pct_res <- visit_statistical_tests(RF_5pct) #17/34 W 17/34 T
visit_RF10pct_res <- visit_statistical_tests(RF_10pct) #18/34 W 17/34 T
visit_RF20pct_res <- visit_statistical_tests(RF_20pct) #15/34 W 17/34 T
visit_RF25pct_res <- visit_statistical_tests(RF_25pct) #17/34 W 17/34 T
visit_RF30pct_res <- visit_statistical_tests(RF_30pct) #20/34 W 21/34 T
visit_RF40pct_res <- visit_statistical_tests(RF_40pct) #20/34 W 18/34 T
#QRILC
visit_QRILC5pct_res <- visit_statistical_tests(QRILC_5pct) #17/34 W 16/34 T
visit_QRILC10pct_res <- visit_statistical_tests(QRILC_10pct) #11/34 W 17/34 T
visit_QRILC20pct_res <- visit_statistical_tests(QRILC_20pct) #0/34 W 7/34 T
visit_QRILC25pct_res <- visit_statistical_tests(QRILC_25pct) #0/34 W 5/34 T
visit_QRILC30pct_res <- visit_statistical_tests(QRILC_30pct) #0/34 W 1/34 T
visit_QRILC40pct_res <- visit_statistical_tests(QRILC_40pct) #0/34 W 0/34 T

#Create data frame for number of significant metabolites
data <- data.frame(
  Method = rep(c("Halfmin", "Halfmin", "Halfmin", "Halfmin", "Halfmin", "Halfmin", "Halfmin",
                 "KNN", "KNN", "KNN", "KNN", "KNN", "KNN", "KNN",
                 "RF", "RF", "RF", "RF", "RF", "RF", "RF",
                 "QRILC", "QRILC", "QRILC", "QRILC", "QRILC", "QRILC", "QRILC"), each = 1),
  Percentage = c("Original", "5%", "10%", "20%", "25%", "30%", "40%",
                 "Original", "5%", "10%", "20%", "25%", "30%", "40%",
                 "Original", "5%", "10%", "20%", "25%", "30%", "40%",
                 "Original", "5%", "10%", "20%", "25%", "30%", "40%"),
  Wilcoxon = c(18, 13, 0, 0, 0, 0, 0,   
               18, 14, 11, 10, 8, 8, 0,  
               18, 17, 18, 15, 17, 20, 20,  
               18, 17, 11, 0, 0, 0, 0),  
  TTest = c(16, 15, 8, 6, 5, 0, 0,  
            16, 15, 9, 10, 10, 7, 10,  
            16, 17, 17, 17, 17, 21, 18,  
            16, 16, 17, 7, 5, 1, 0)  
)

# Convert Percentage to factor to maintain ordering
data$Percentage <- factor(data$Percentage, levels = c("Original", "5%", "10%", "20%", "25%", "30%", "40%"))

# Reshape data for plotting
data_long <- data %>% 
  pivot_longer(cols = c("Wilcoxon", "TTest"), names_to = "Test", values_to = "Significant_Metabolites")

# Create grouped bar plot ensuring "Original" appears first in each facet
ggplot(data_long, aes(x = Percentage, y = Significant_Metabolites, fill = Test)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Significant_Metabolites), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 5) +  # Add exact numbers above bars
  facet_wrap(~ Method, scales = "free_x", nrow = 2) +  # Ensure each imputation method has its own facet
  labs(title = "Number of Significant Metabolites (Visit 1 vs Visit 2)",
       x = "Missingness Percentage",
       y = "Number of Significant Metabolites",
       fill = "Statistical Test") +
  theme_minimal(base_size = 14) +  # Increase font size
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 14), # Increase facet label size
        panel.spacing = unit(2, "lines")) # More spacing between facets


# ------------------------------------
# Part 3: NRMSE
# ------------------------------------

#function to calculate NRMSE
calculate_weighted_nrmse <- function(original, imputed, method, percentage) {
  #take numeric columns 
  numeric_col_names <- names(original)[6:ncol(original)]
  
  #calcualte for each column nrmse
  nrmse_values <- sapply(numeric_col_names, function(col){
    actual_val <- original[[col]]
    predicted_val <- imputed[[col]]
    
    #ensure no missing values for comparison
    valid_indices <- !is.na(actual_val) & !is.na(predicted_val)
    
    if (sum(valid_indices) > 2) {  #if enough data points
      mse <- mean((actual_val[valid_indices] - predicted_val[valid_indices])^2) #mean squared error
      rmse <- sqrt(mse) #root mean squared error
      norm_factor <- max(actual_val[valid_indices], na.rm = TRUE) - min(actual_val[valid_indices], na.rm = TRUE) #denominator 
      
      if (norm_factor > 0) {
        nrmse <- rmse / norm_factor #NRMSE
        weighted_nrmse <- nrmse * percentage
        return(weighted_nrmse)
      } else {
        return(NA) #avoid dividing by 0
      }
    } else {
      return(NA)  #return NA else (if not > 2)
    }
  })
  
  return(data.frame(
    Metabolite = numeric_col_names, 
    Imputation_Method = method,
    MCAR_Proportion = (percentage*100),
    Weighted_NRMSE = nrmse_values
  ))
}

#call function to compute NRMSE
#Halfmin
nrmse_res_halfmin_5pct <- calculate_weighted_nrmse(FAO_original, Halfmin_5pct, "Halfmin", 0.05)
nrmse_res_halfmin_10pct <- calculate_weighted_nrmse(FAO_original, Halfmin_10pct, "Halfmin", 0.1)
nrmse_res_halfmin_20pct <- calculate_weighted_nrmse(FAO_original, Halfmin_20pct,"Halfmin", 0.2)
nrmse_res_halfmin_25pct <- calculate_weighted_nrmse(FAO_original, Halfmin_25pct, "Halfmin", 0.25)
nrmse_res_halfmin_30pct <- calculate_weighted_nrmse(FAO_original, Halfmin_30pct, "Halfmin", 0.3)
nrmse_res_halfmin_40pct <- calculate_weighted_nrmse(FAO_original, Halfmin_40pct, "Halfmin", 0.4)
#KNN
nrmse_res_KNN_5pct <- calculate_weighted_nrmse(FAO_original, KNN_5pct, "KNN", 0.05)
nrmse_res_KNN_10pct <- calculate_weighted_nrmse(FAO_original, KNN_10pct, "KNN", 0.1)
nrmse_res_KNN_20pct <- calculate_weighted_nrmse(FAO_original, KNN_20pct, "KNN", 0.2)
nrmse_res_KNN_25pct <- calculate_weighted_nrmse(FAO_original, KNN_25pct, "KNN", 0.25)
nrmse_res_KNN_30pct <- calculate_weighted_nrmse(FAO_original, KNN_30pct, "KNN", 0.3)
nrmse_res_KNN_40pct <- calculate_weighted_nrmse(FAO_original, KNN_40pct, "KNN", 0.4)
#RF
nrmse_res_RF_5pct <- calculate_weighted_nrmse(FAO_original, RF_5pct, "RF", 0.05)
nrmse_res_RF_10pct <- calculate_weighted_nrmse(FAO_original, RF_10pct, "RF", 0.1)
nrmse_res_RF_20pct <- calculate_weighted_nrmse(FAO_original, RF_20pct, "RF", 0.2)
nrmse_res_RF_25pct <- calculate_weighted_nrmse(FAO_original, RF_25pct, "RF", 0.25)
nrmse_res_RF_30pct <- calculate_weighted_nrmse(FAO_original, RF_30pct, "RF", 0.3)
nrmse_res_RF_40pct <- calculate_weighted_nrmse(FAO_original, RF_40pct, "RF", 0.4)
#QRILC
nrmse_res_QRILC_5pct <- calculate_weighted_nrmse(FAO_original, QRILC_5pct, "QRILC", 0.05)
nrmse_res_QRILC_10pct <- calculate_weighted_nrmse(FAO_original, QRILC_10pct, "QRILC", 0.1)
nrmse_res_QRILC_20pct <- calculate_weighted_nrmse(FAO_original, QRILC_20pct, "QRILC", 0.2)
nrmse_res_QRILC_25pct <- calculate_weighted_nrmse(FAO_original, QRILC_25pct, "QRILC", 0.25)
nrmse_res_QRILC_30pct <- calculate_weighted_nrmse(FAO_original, QRILC_30pct, "QRILC", 0.3)
nrmse_res_QRILC_40pct <- calculate_weighted_nrmse(FAO_original, QRILC_40pct, "QRILC", 0.4)

# ------------------------------------
# Part 3.1: NRMSE Plot 
# ------------------------------------

#for plotting combine all nrmse results into one dataframe
nrmse_data <- bind_rows(
  nrmse_res_halfmin_5pct, nrmse_res_halfmin_10pct, nrmse_res_halfmin_20pct, 
  nrmse_res_halfmin_25pct, nrmse_res_halfmin_30pct, nrmse_res_halfmin_40pct,
  
  nrmse_res_KNN_5pct, nrmse_res_KNN_10pct, nrmse_res_KNN_20pct, 
  nrmse_res_KNN_25pct, nrmse_res_KNN_30pct, nrmse_res_KNN_40pct,
  
  nrmse_res_RF_5pct, nrmse_res_RF_10pct, nrmse_res_RF_20pct, 
  nrmse_res_RF_25pct, nrmse_res_RF_30pct, nrmse_res_RF_40pct,
  
  nrmse_res_QRILC_5pct, nrmse_res_QRILC_10pct, nrmse_res_QRILC_20pct, 
  nrmse_res_QRILC_25pct, nrmse_res_QRILC_30pct, nrmse_res_QRILC_40pct
)

#convert MCAR_Proportion to a factor for correct ordering in the plot
nrmse_data$MCAR_Proportion <- as.factor(nrmse_data$MCAR_Proportion)

#Generate Boxplot 
ggplot(nrmse_data, aes(x = MCAR_Proportion, y = Weighted_NRMSE, fill = Imputation_Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + #boxplot wiht transparency and w/o outliers
  scale_fill_manual(values = c("lightblue", "orange", "blue", "magenta")) +
  labs(
    title = "Weighted NRMSE across Imputation Method and MCAR Proportions",
    x = "MCAR Proportion (%)",
    y = "Weigthed NRMSE", 
    fill = "Imputation Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Keep legend for clarity
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "gray85"),  # Light gray grid
    panel.grid.minor = element_blank()  # Remove minor grid
  ) +
  ylim(0,0.4)


# ------------------------------------
# Part 4: Normalized Mean Difference
# ------------------------------------

#function to calculate normalized mean difference
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
norm_diff_Halfmin_10pct <- norm_mean_diff(FAO_original, Halfmin_10pct, "Half-min", 10)
norm_diff_Halfmin_20pct <- norm_mean_diff(FAO_original, Halfmin_20pct, "Half-min", 20)
norm_diff_Halfmin_25pct <- norm_mean_diff(FAO_original, Halfmin_25pct, "Half-min", 25)
norm_diff_Halfmin_30pct <- norm_mean_diff(FAO_original, Halfmin_30pct, "Half-min", 30)
norm_diff_Halfmin_40pct <- norm_mean_diff(FAO_original, Halfmin_40pct, "Half-min", 40)
#KNN
norm_diff_KNN_5pct <- norm_mean_diff(FAO_original, KNN_5pct, "KNN", 5)
norm_diff_KNN_10pct <- norm_mean_diff(FAO_original, KNN_10pct, "KNN", 10)
norm_diff_KNN_20pct <- norm_mean_diff(FAO_original, KNN_20pct, "KNN", 20)
norm_diff_KNN_25pct <- norm_mean_diff(FAO_original, KNN_25pct, "KNN", 25)
norm_diff_KNN_30pct <- norm_mean_diff(FAO_original, KNN_30pct, "KNN", 30)
norm_diff_KNN_40pct <- norm_mean_diff(FAO_original, KNN_40pct, "KNN", 40)
#RF
norm_diff_RF_5pct <- norm_mean_diff(FAO_original, RF_5pct, "RF", 5)
norm_diff_RF_10pct <- norm_mean_diff(FAO_original, RF_10pct, "RF", 10)
norm_diff_RF_20pct <- norm_mean_diff(FAO_original, RF_20pct, "RF", 20)
norm_diff_RF_25pct <- norm_mean_diff(FAO_original, RF_25pct, "RF", 25)
norm_diff_RF_30pct <- norm_mean_diff(FAO_original, RF_30pct, "RF", 30)
norm_diff_RF_40pct <- norm_mean_diff(FAO_original, RF_40pct, "RF", 40)
#QRILC
norm_diff_QRILC_5pct <- norm_mean_diff(FAO_original, QRILC_5pct, "QRILC", 5)
norm_diff_QRILC_10pct <- norm_mean_diff(FAO_original, QRILC_10pct, "QRILC", 10)
norm_diff_QRILC_20pct <- norm_mean_diff(FAO_original, QRILC_20pct, "QRILC", 20)
norm_diff_QRILC_25pct <- norm_mean_diff(FAO_original, QRILC_25pct, "QRILC", 25)
norm_diff_QRILC_30pct <- norm_mean_diff(FAO_original, QRILC_30pct, "QRILC", 30)
norm_diff_QRILC_40pct <- norm_mean_diff(FAO_original, QRILC_40pct, "QRILC", 40)

# ------------------------------------
# Part 4.1: Plot of all NMD per Imputation Method
# ------------------------------------

#function to compute and return normalized mean difference data
norm_mean_diff_data <- function(original, imputed, method, percentage) {
  #numeric columns
  original_numeric <- original[, 6:ncol(original)]
  imputed_numeric <- imputed[, 6:ncol(imputed)]
  
  #compute mean before imputation
  mean_before <- original_numeric %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_Before")
  
  #compute mean after imputation
  mean_after <- imputed_numeric %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Mean_After")
  
  #merge before and after mean values
  mean_comparison <- left_join(mean_before, mean_after, by = "Metabolite")
  
  #compute normalized difference: (Mean_After - Mean_Before) / Mean_Before
  mean_comparison <- mean_comparison %>%
    mutate(Normalized_Difference = (Mean_After - Mean_Before) / Mean_Before,
           Method = method,
           Percentage = paste0(percentage, "%"))
  
  return(mean_comparison)
}

#combine all data for each imputation method
halfmin_data <- bind_rows(
  norm_mean_diff_data(FAO_original, Halfmin_5pct, "Half-min", 5),
  norm_mean_diff_data(FAO_original, Halfmin_10pct, "Half-min", 10),
  norm_mean_diff_data(FAO_original, Halfmin_20pct, "Half-min", 20),
  norm_mean_diff_data(FAO_original, Halfmin_25pct, "Half-min", 25),
  norm_mean_diff_data(FAO_original, Halfmin_30pct, "Half-min", 30),
  norm_mean_diff_data(FAO_original, Halfmin_40pct, "Half-min", 40)
)

knn_data <- bind_rows(
  norm_mean_diff_data(FAO_original, KNN_5pct, "KNN", 5),
  norm_mean_diff_data(FAO_original, KNN_10pct, "KNN", 10),
  norm_mean_diff_data(FAO_original, KNN_20pct, "KNN", 20),
  norm_mean_diff_data(FAO_original, KNN_25pct, "KNN", 25),
  norm_mean_diff_data(FAO_original, KNN_30pct, "KNN", 30),
  norm_mean_diff_data(FAO_original, KNN_40pct, "KNN", 40)
)

rf_data <- bind_rows(
  norm_mean_diff_data(FAO_original, RF_5pct, "RF", 5),
  norm_mean_diff_data(FAO_original, RF_10pct, "RF", 10),
  norm_mean_diff_data(FAO_original, RF_20pct, "RF", 20),
  norm_mean_diff_data(FAO_original, RF_25pct, "RF", 25),
  norm_mean_diff_data(FAO_original, RF_30pct, "RF", 30),
  norm_mean_diff_data(FAO_original, RF_40pct, "RF", 40)
)

qrilc_data <- bind_rows(
  norm_mean_diff_data(FAO_original, QRILC_5pct, "QRILC", 5),
  norm_mean_diff_data(FAO_original, QRILC_10pct, "QRILC", 10),
  norm_mean_diff_data(FAO_original, QRILC_20pct, "QRILC", 20),
  norm_mean_diff_data(FAO_original, QRILC_25pct, "QRILC", 25),
  norm_mean_diff_data(FAO_original, QRILC_30pct, "QRILC", 30),
  norm_mean_diff_data(FAO_original, QRILC_40pct, "QRILC", 40)
)

#function to plot density for each method
plot_density <- function(data, method) {
  ggplot(data, aes(x = Normalized_Difference, fill = Percentage, color = Percentage)) +
    geom_density(alpha = 0.3) +
    theme_minimal() +
    labs(title = paste("Normalized Difference for", method, "Imputation"),
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.2, 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    theme(legend.title = element_blank(), legend.position = "right")
}

#generate plots
plot_halfmin <- plot_density(halfmin_data, "Half-min")
plot_knn <- plot_density(knn_data, "KNN")
plot_rf <- plot_density(rf_data, "RF")
plot_qrilc <- plot_density(qrilc_data, "QRILC")

#display plots
print(plot_halfmin)
print(plot_knn)
print(plot_rf)
print(plot_qrilc)


# ------------------------------------
# Part 5: Plot Distribution
# ------------------------------------

#function to plot the distribution before and after imputation and only imputed values 
plot_distribution <- function(original, imputed, method, percentage) {
  #numeric columns 
  numeric_original <- original[, 6:ncol(original)]
  numeric_imputed <- imputed[, 6:ncol(imputed)]
  
  #convert to long format for plotting
  original_long <- numeric_original %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    mutate(Data = "Original")
  
  imputed_long <- numeric_imputed %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    mutate(Data = "Imputed")
  
  #combine data
  combined_data <- bind_rows(original_long, imputed_long)
  
  #caclualte mean 
  mean_data <- combined_data %>%
    group_by(Metabolite, Data) %>%
    summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  #plot 
  plot <- ggplot(combined_data, aes(x = Value, fill = Data)) +
    geom_density(alpha = 0.5, size = 1) + #transparency
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal() +
    labs(title = paste0("Density Distribution Before and After Imputation (", method, ", ", percentage, "% Missing)"),
         x = "Value",
         y = "Density") +
    #add mean lines
    geom_vline(data = mean_data %>% filter(Data == "Original"),
               aes(xintercept = mean_value), color = "blue", linewidth = 0.8, linetype = "dashed") +
    geom_vline(data = mean_data %>% filter(Data == "Imputed"),
               aes(xintercept = mean_value), color = "red", linewidth = 0.8, linetype = "dashed") +
    theme(legend.position = "bottom")
  
  
  print(plot)
  return(combined_data)
  
}

#call function to plot distirbution
#Halfmin
dist_Halfmin_5pct <- plot_distribution(FAO_original, Halfmin_5pct, "Half-min", 5)
dist_Halfmin_10pct <- plot_distribution(FAO_original, Halfmin_10pct, "Half-min", 10)
dist_Halfmin_20pct <- plot_distribution(FAO_original, Halfmin_20pct, "Half-min", 20)
dist_Halfmin_25pct <- plot_distribution(FAO_original, Halfmin_25pct, "Half-min", 25)
dist_Halfmin_30pct <- plot_distribution(FAO_original, Halfmin_30pct, "Half-min", 30)
dist_Halfmin_40pct <- plot_distribution(FAO_original, Halfmin_40pct, "Half-min", 40)
#KNN
dist_KNN_5pct <- plot_distribution(FAO_original, KNN_5pct, "KNN", 5)
dist_KNN_10pct <- plot_distribution(FAO_original, KNN_10pct, "KNN", 10)
dist_KNN_20pct <- plot_distribution(FAO_original, KNN_20pct, "KNN", 20)
dist_KNN_25pct <- plot_distribution(FAO_original, KNN_25pct, "KNN", 25)
dist_KNN_30pct <- plot_distribution(FAO_original, KNN_30pct, "KNN", 30)
dist_KNN_40pct <- plot_distribution(FAO_original, KNN_40pct, "KNN", 40)
#RF
dist_RF_5pct <- plot_distribution(FAO_original, RF_5pct, "RF", 5)
dist_RF_10pct <- plot_distribution(FAO_original, RF_10pct, "RF", 10)
dist_RF_20pct <- plot_distribution(FAO_original, RF_20pct, "RF", 20)
dist_RF_25pct <- plot_distribution(FAO_original, RF_25pct, "RF", 25)
dist_RF_30pct <- plot_distribution(FAO_original, RF_30pct, "RF", 30)
dist_RF_40pct <- plot_distribution(FAO_original, RF_40pct, "RF", 40)
#QRILC
dist_QRILC_5pct <- plot_distribution(FAO_original, QRILC_5pct, "QRILC", 5)
dist_QRILC_10pct <- plot_distribution(FAO_original, QRILC_10pct, "QRILC", 10)
dist_QRILC_20pct <- plot_distribution(FAO_original, QRILC_20pct, "QRILC", 20)
dist_QRILC_25pct <- plot_distribution(FAO_original, QRILC_25pct, "QRILC", 25)
dist_QRILC_30pct <- plot_distribution(FAO_original, QRILC_30pct, "QRILC", 30)
dist_QRILC_40pct <- plot_distribution(FAO_original, QRILC_40pct, "QRILC", 40)

# ------------------------------------
# Part 5.1: Distribution Plot (Whole Dataset)
# ------------------------------------

#function to plot distribution before and after imputation (entire dataset)
plot_whole_distribution <- function(original, imputed, method, percentage) {
  #numeric columns 
  numeric_original <- original[, 6:ncol(original)]
  numeric_imputed <- imputed[, 6:ncol(imputed)]
  
  #convert to long format for plotting
  original_long <- numeric_original %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    mutate(Data = "Original Data")
  
  imputed_long <- numeric_imputed %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    mutate(Data = "Imputed Data")
  
  #identify imputed values
  imputed_values <- numeric_original != numeric_imputed
  
  imputed_only_long <- numeric_imputed %>%
    as.data.frame() %>%
    replace(!imputed_values, NA) %>%  #keep only changed (imputed) values
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    filter(!is.na(Value)) %>%
    mutate(Data = "Imputed Values")  
  
  #combine both data
  combined_data <- bind_rows(original_long, imputed_long, imputed_only_long)
  
  #compute mean
  mean_data <- combined_data %>%
    group_by(Data) %>%
    summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  
  #plot overall density distribution
  plot <- ggplot(combined_data, aes(x = Value, fill = Data)) +
    geom_density(alpha = 0.5) + 
    theme_minimal() +
    labs(title = paste0("Overall Density Distribution Before and After Imputation (", method, ", ", percentage, "% Missing)"),
         x = "Value",
         y = "Density") +
    geom_vline(data = mean_data %>% filter(Data == "Original Data"),
               aes(xintercept = mean_value), color = "blue", linewidth = 1, linetype = "dashed") + 
    geom_vline(data = mean_data %>% filter(Data == "Imputed Data"),
               aes(xintercept = mean_value), color = "red", linewidth = 1, linetype = "dashed") +
    xlim(-100, 400) + 
    theme(legend.position = "bottom")
  
  print(plot)
  return(combined_data)
}

#call function to plot whole distirbution
#Halfmin
whole_dist_Halfmin_5pct <- plot_whole_distribution(FAO_original, Halfmin_5pct, "Half-min", 5)
whole_dist_Halfmin_40pct <- plot_whole_distribution(FAO_original, Halfmin_40pct, "Half-min", 40)
#KNN
whole_dist_KNN_5pct <- plot_whole_distribution(FAO_original, KNN_5pct, "KNN", 5)
whole_dist_KNN_40pct <- plot_whole_distribution(FAO_original, KNN_40pct, "KNN", 40)
#RF
whole_dist_RF_5pct <- plot_whole_distribution(FAO_original, RF_5pct, "RF", 5)
whole_dist_RF_40pct <- plot_whole_distribution(FAO_original, RF_40pct, "RF", 40)
#QRILC
whole_dist_QRILC_5pct <- plot_whole_distribution(FAO_original, QRILC_5pct, "QRILC", 5)
whole_dist_QRILC_40pct <- plot_whole_distribution(FAO_original, QRILC_40pct, "QRILC", 40)

# ------------------------------------
# Part 6: ANOVA
# ------------------------------------

#fit an ANOVA model with an interaction term between imputation method and missingness level 
#apply log to normalize data
anova_model_mcar <- aov(log(Weighted_NRMSE) ~ Imputation_Method * MCAR_Proportion, data = nrmse_data)

#show ANOVA summary
summary(anova_model_mcar)

# For Imputation_Method
# --------
#3 df (degree of freedom) means 4 imputation method - 1 
#Sum Sq = 201.8 (variation explained by impuation method)
#Mean Sq = 67.26 (variation per degree of freedom)
#F-value = 89.048 (string effect if high value)
#P-value highly significant < 2e-16
#Significant effect of imputation method on NRMSE

# For MCAR_Proportion
# --------
#5 df (degree of freedom) means 6 levels of missing data 
#Sum Sq = 1159.5 (variation explained by impuation method)
#Mean Sq = 231.90 (variation per degree of freedom)
#F-value = 307.036 (strong effect, very high value)
#P-value highly significant < 2e-16
#Proportion of missing data has a strong signficnant effect on NRMSE

# For interaction of both
# -------
#P-value is 1, so no signification i.e. there is no interaction effect meaning that the effect of imputation methond on NRMSE does not depend on the missingness level 
#The effect of imputation method is consistens across different missingness levels

# --------------------------------------------------------
# Part 6.1: Check Residuals and Normality for ANOVA result
# -------------------------------------------------------

#check residuals for normality
#histogram of residuals (extracts results from anova model)
#residul look symmetry and a bit bell shaped then it suggests normalizy 
ggplot(data.frame(residuals = residuals(anova_model_mcar)), aes(x = residuals)) +
  geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency") +
  theme_minimal()

#Q-Q plot of residuals
#plot residula against theoretical normal distirbutions
#red line shows perfect normal distirbution
#S-shaped pattern = possibles skewness
qqnorm(residuals(anova_model_mcar), col = "blue")
qqline(residuals(anova_model_mcar), col = "red")

#Tukey-Anscombe plot to check residuals vs fitted values
#x-axis = fitted values (predicted) and y-axis = residueals (error)
#residuals should be evenly spread, if the points fan out or forma pattern the assumption of homoscedascity is violated 
plot(fitted(anova_model_mcar), resid(anova_model_mcar), 
     main = "Tukey-Anscombe Plot", 
     col = "blue", 
     xlab = "Fitted Values (Predicted by ANOVA Model)", 
     ylab = "Residuals (Errors)")

# ------------------------------------
# Part 7: Kruskal-Wallis
# ------------------------------------

#perform kruskal-wallis test
kruskal_test <- kruskal.test(Weighted_NRMSE ~ Imputation_Method, data = nrmse_data)
kruskal_test2 <- kruskal.test(Weighted_NRMSE ~ MCAR_Proportion, data = nrmse_data)
#show results
print(kruskal_test)
print(kruskal_test2) 
#p-value storngly signficant 
#chi-square = 98.713 (higher value means larger difference between groups)
#at least one imputation method significantly differs from the others in terms of NRMSE

# ------------------------------------
# Part 8: Dunn's Test for each imputation
# ------------------------------------

#perform Dunn's Test for pairwise comparison (BH correction for multiple testing)
dunn_test <- dunnTest(Weighted_NRMSE ~ Imputation_Method, data = nrmse_data, method = "bh")

#print the results
print(dunn_test)

#plot Dunns test results
ggplot(nrmse_data, aes(x = Imputation_Method, y = Weighted_NRMSE, fill = Imputation_Method)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Pairwise Comparisons of Imputation Methods (Dunn's Test)",
       x = "Imputation Method",
       y = "Weighted NRMSE") +
  ylim(0,1)

