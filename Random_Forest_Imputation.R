#load necessary library
library(missForest)
library(dplyr)
library(corrplot)

# Part 1 -------
# Random Forest impputation of MV

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

# Part 2 --------
# Wilcoxon rank-sum test (for independent data)

#get metabolite column names
metabolite_cols <- colnames(imputed_RF_23)

#dataframe for the results
results_Wilcoxon_RF <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)


#loop through each metabolite and run Wilcoxon rank-sum test
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  test_result <- wilcox.test(imputed_RF_23[[metabolite]], imputed_RF_24[[metabolite]], paired = FALSE, exact = FALSE)
  
  #save results to dataframe
  results_Wilcoxon_RF$p_value[i] <- test_result$p.value
  results_Wilcoxon_RF$statistic[i] <- test_result$statistic
}


#adjust p-value for multiple testing (Benjamini-Hochberg
results_Wilcoxon_RF$adj_p_value <- p.adjust(results_Wilcoxon_RF$p_value, method = "BH")

print(results_Wilcoxon_RF)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_Wilcoxon_RF <- results_Wilcoxon_RF %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_Wilcoxon_RF) #82 out of 84 were significant 

# Part 3 -------
# Run unpaired t-test for each metabolite 
#dataframe for the results
results_ttest_RF <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)

#loop through each metabolite and run t-test for each metabolites
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  #use Welchs test here (not same variance, to test H0 where two groups have equal mean)
  test_result2 <- t.test(imputed_RF_23[[metabolite]], imputed_RF_24[[metabolite]], paired = FALSE, var.equal = FALSE)
  
  #save results to dataframe
  results_ttest_RF$p_value[i] <- test_result2$p.value
  results_ttest_RF$statistic[i] <- test_result2$statistic
}

#adjust p-value for multiple testing (Benjamini-Hochberg)
results_ttest_RF$adj_p_value <- p.adjust(results_ttest_RF$p_value, method = "BH")

print(results_ttest_RF)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_ttest_RF <- results_ttest_RF %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_ttest_RF) #80 out of 84 were significant

# Part 4 ----- 
# Check if Impuation method was good 

#OOB (out of bag imputation error) is the normalized mean squared error (NMSE) => the lower the better
#NMSE = 0 means perfect prediction
OOB23 <- RF_data23$OOBerror #0.005697
OOB24 <- RF_data24$OOBerror #0.40187 higher error, less effecitve imputation 

# Correlation analysis 
cor_before_23 <- cor(numeric_23, use = "pairwise.complete.obs")
cor_after_23 <- cor(imputed_RF_23)

# Heat map
par(mfrow = c(1,2)) #for side by side plots
corrplot(cor_before_23, method = "color", tl.cex = 0.6, title = "Before Imputation (2023)")
corrplot(cor_after_23, method = "color", tl.cex = 0.6, title = "After Imputation (2023)")
dev.off()

# Scatter plot 
#identify missing vlaues 
missing_val <- is.na(numeric_23) #total of 99 are missing

#create logical matrix identifying if a pair of metabolites had missing values
missing_pairs <- (missing_val %*% t(missing_val)) > 0 #True 

#extract upper triangular part of matrix 
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

