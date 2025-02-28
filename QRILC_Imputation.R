#install packages required for QRILC
#install.packages("imputeLCMD")

#load necessary libraries
library(imputeLCMD)
library(tidyverse)
library(corrplot)

# Part 1 -------
# QRILC Imputation

#whole data 
data23 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites23.csv", check.names = FALSE)
data24 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites23.csv", check.names = FALSE)

#remove metadata from whole data to get numeric data 
numeric23 <- data23[, 6:ncol(data23)]
numeric24 <- data24[, 6:ncol(data24)]

#needs transformation to matrix 
numeric23 <- as.matrix(numeric23)
numeric24 <- as.matrix(numeric24)

#Log transformation (helps for very small values, prevents negative imputed values)
log_data23 <- log(numeric23 + 1e-6) #add small constant to avoid log(0)
log_data24 <- log(numeric24 + 1e-6) 

#QRILC imputation 
QRILC_impuation <- function(data, ...) {
  imputed_data <- impute.QRILC(data)[[1]] #returns list, extract imputed matrix
  return(imputed_data)
}

#call function for imputation 
imputed_QRILC23 <- QRILC_impuation(log_data23)
imputed_QRILC24 <- QRILC_impuation(log_data24)

#convert back to data frame and exponentiate back 
exp_QRILC23 <- exp(imputed_QRILC23) - 1e-6
exp_QRILC24 <- exp(imputed_QRILC24) - 1e-6
QRILC23 <- as.data.frame(exp_QRILC23)
QRILC24 <- as.data.frame(exp_QRILC24)


# Part 2 --------
# Wilcoxon rank-sum test (for independent data)

#get metabolite column names
metabolite_cols <- colnames(QRILC23)

#dataframe for the results
results_Wilcoxon_QRILC <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)


#loop through each metabolite and run Wilcoxon rank-sum test
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  test_result <- wilcox.test(QRILC23[[metabolite]], QRILC24[[metabolite]], paired = FALSE, exact = FALSE)
  
  #save results to dataframe
  results_Wilcoxon_QRILC$p_value[i] <- test_result$p.value
  results_Wilcoxon_QRILC$statistic[i] <- test_result$statistic
}


#adjust p-value for multiple testing (Benjamini-Hochberg
results_Wilcoxon_QRILC$adj_p_value <- p.adjust(results_Wilcoxon_QRILC$p_value, method = "BH")

print(results_Wilcoxon_QRILC)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_Wilcoxon_QRILC <- results_Wilcoxon_QRILC %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_Wilcoxon_QRILC) #0 out of 84 were significant 

# Part 3 -------
# Run unpaired t-test for each metabolite 
#dataframe for the results
results_ttest_QRILC <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)

#loop through each metabolite and run t-test for each metabolites
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  #use Welchs test here (not same variance, to test H0 where two groups have equal mean)
  test_result2 <- t.test(QRILC23[[metabolite]], QRILC24[[metabolite]], paired = FALSE, var.equal = FALSE)
  
  #save results to dataframe
  results_ttest_QRILC$p_value[i] <- test_result2$p.value
  results_ttest_QRILC$statistic[i] <- test_result2$statistic
}

#adjust p-value for multiple testing (Benjamini-Hochberg)
results_ttest_QRILC$adj_p_value <- p.adjust(results_ttest_QRILC$p_value, method = "BH")

print(results_ttest_QRILC)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_ttest_QRILC <- results_ttest_QRILC %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_ttest_QRILC) #0 out of 84 were significant


# Part 4 ----- 
# Check if Impuation method was good 

# Correlation analysis 
cor_before_23 <- cor(numeric23, use = "pairwise.complete.obs")
cor_after_23 <- cor(QRILC23)
cor_before_24 <- cor(numeric24, use = "pairwise.complete.obs")
cor_after_24 <- cor(QRILC24)

# Heat map
par(mfrow = c(1,2)) #for side by side plots
corrplot(cor_before_23, method = "color", tl.cex = 0.6, title = "Before Imputation (2023)")
corrplot(cor_after_23, method = "color", tl.cex = 0.6, title = "After Imputation (2023)")
dev.off()

# Scatter plot 
#identify missing vlaues 
missing_val <- is.na(numeric23) #total of 99 are missing
missing_val2 <- is.na(numeric24)

#create logical matrix identifying if a pair of metabolites had missing values
missing_pairs <- (missing_val %*% t(missing_val)) > 0 #True 
missing_pairs2 <- (missing_val2 %*% t(missing_val2)) > 0 #True

#extract upper triangular part of matrix 2023
imputed_indices <- missing_pairs[upper.tri(missing_pairs, diag = FALSE)]
cor_values_before23 <- cor_before_23[upper.tri(cor_before_23, diag = FALSE)]
cor_values_after23 <- cor_after_23[upper.tri(cor_after_23, diag = FALSE)]

plot(cor_values_before23, cor_values_after23,
     xlab = "Before Imputation",
     ylab = "After Imputation",
     main = "Correlation Comparison",
     col = ifelse(imputed_indices, "red", "blue"), #red for imputed pairs
     pch = 19) 

abline(0,1,col = "black", lwd = 2) #ideally points should line on line

#extract upper triangular part of matrix 2024
imputed_indices2 <- missing_pairs2[upper.tri(missing_pairs2, diag = FALSE)]
cor_values_before24 <- cor_before_24[upper.tri(cor_before_24, diag = FALSE)]
cor_values_after24 <- cor_after_24[upper.tri(cor_after_24, diag = FALSE)]

plot(cor_values_before24, cor_values_after24,
     xlab = "Before Imputation",
     ylab = "After Imputation",
     main = "Correlation Comparison",
     col = ifelse(imputed_indices2, "red", "blue"), #red for imputed pairs
     pch = 19) 

abline(0,1,col = "black", lwd = 2) #ideally points should line on line

