#load library
library(impute) #for impute.knn() function 
library(corrplot) #for correlaion coefficient calculation

# Part 1 -------
# Perform KNN imputation of Missing values 

#load data 
data23 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites23.csv", check.names = FALSE)
data24 <- read.csv("/Users/marcinebessire/Desktop/project/Common_Metabolites24.csv", check.names = FALSE)

#numeric data 
numeric23 <- data23[, 6:ncol(data23)]
numeric24 <- data24[, 6:ncol(data24)]

#need to transpose because imput.knn() has samples in columns and metabolites in rows
#impute.knn() needs matrix format, because of mathematical operations (euclidean distances)

#2023
imputed_KNN23 <- impute.knn(as.matrix(t(numeric23)), rowmax = 0.5, colmax = 1) #transform first
imputed_KNN23 <- as.data.frame(t(imputed_KNN23$data)) #transform back and use transformed data 

#2024
imputed_KNN24 <- impute.knn(as.matrix(t(numeric24)), rowmax = 0.5, colmax = 1)
imputed_KNN24 <- as.data.frame(t(imputed_KNN24$data))

# Part 2 --------
# Wilcoxon rank-sum test (for independent data)

#get metabolite column names
metabolite_cols <- colnames(imputed_KNN23)

#dataframe for the results
results_Wilcoxon_KNN <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)


#loop through each metabolite and run Wilcoxon rank-sum test
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  test_result <- wilcox.test(imputed_KNN23[[metabolite]], imputed_KNN24[[metabolite]], paired = FALSE, exact = FALSE)
  
  #save results to dataframe
  results_Wilcoxon_KNN$p_value[i] <- test_result$p.value
  results_Wilcoxon_KNN$statistic[i] <- test_result$statistic
}


#adjust p-value for multiple testing (Benjamini-Hochberg
results_Wilcoxon_KNN$adj_p_value <- p.adjust(results_Wilcoxon_KNN$p_value, method = "BH")

print(results_Wilcoxon_KNN)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_Wilcoxon_KNN <- results_Wilcoxon_KNN %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_Wilcoxon_KNN) #81 out of 84 were significant 

# Part 3 -------
# Run unpaired t-test for each metabolite 
#dataframe for the results
results_ttest_KNN <- data.frame(
  Metabolite = metabolite_cols,
  p_value = numeric(length(metabolite_cols)),
  statistic = numeric(length(metabolite_cols))
)

#loop through each metabolite and run t-test for each metabolites
for (i in seq_along(metabolite_cols)) {
  metabolite <- metabolite_cols[i]
  
  #use Welchs test here (not same variance, to test H0 where two groups have equal mean)
  test_result2 <- t.test(imputed_KNN23[[metabolite]], imputed_KNN24[[metabolite]], paired = FALSE, var.equal = FALSE)
  
  #save results to dataframe
  results_ttest_KNN$p_value[i] <- test_result2$p.value
  results_ttest_KNN$statistic[i] <- test_result2$statistic
}

#adjust p-value for multiple testing (Benjamini-Hochberg)
results_ttest_KNN$adj_p_value <- p.adjust(results_ttest_KNN$p_value, method = "BH")

print(results_ttest_KNN)

#view significant metabolites with BH adjusted p-value < 0.05
significant_results_ttest_KNN <- results_ttest_KNN %>% 
  filter(adj_p_value < 0.05) %>% #usually 0.05 used 
  arrange(adj_p_value)

print(significant_results_ttest_KNN) #6 out of 84 were significant


# Part 4 ----
# Assess quality of KNN imputation 

#first check if all values were imputed (no NA left)
sum(is.na(imputed_KNN23)) #is 0
sum(is.na(imputed_KNN24)) #is 0 

#compare data distribution before and after imputation 
par(mfrow = c(3,3)) #3x3 grid
par(mar = c(4,4,2,1))

#Boxplots of 9 metabolites (because 84 too big) for 2023
for (i in 1:9) {
  boxplot(imputed_KNN23[[i]], numeric23[[i]],
          names = c("KNN Imputation", "No Imputation"),
          main = paste("Metabolite", colnames(imputed_KNN23)[i]),
          col = c("lightblue", "lightgreen"))
}

#Boxplots of 9 metabolites (because 84 too big) for 2024
for (i in 1:9) {
  boxplot(imputed_KNN24[[i]], numeric24[[i]],
          names = c("KNN Imputation", "No Imputation"),
          main = paste("Metabolite", colnames(imputed_KNN23)[i]),
          col = c("lightblue", "lightgreen"))
}

#Correlation Check, compare correlation matrices
cor_before_23 <- cor(numeric23, use = "pairwise.complete.obs")
cor_after_23 <- cor(imputed_KNN23)

cor_before_24 <- cor(numeric24, use = "pairwise.complete.obs")
cor_after_24 <- cor(imputed_KNN24)

# Scatter plot 
#identify missing vlaues 
missing_val23 <- is.na(numeric23) #99
missing_val24 <- is.na(numeric24) #526


#create logical matrix identifying if a pair of metabolites had missing values
missing_pairs23 <- (missing_val23 %*% t(missing_val23)) > 0 #True 
missing_pairs24 <- (missing_val24 %*% t(missing_val24)) > 0

#extract upper triangular part of matrix 
#2023
imputed_indices23 <- missing_pairs23[upper.tri(missing_pairs23, diag = FALSE)]
cor_values_before23 <- cor_before_23[upper.tri(cor_before_23, diag = FALSE)]
cor_values_after23 <- cor_after_23[upper.tri(cor_after_23, diag = FALSE)]
#2024
imputed_indices24 <- missing_pairs24[upper.tri(missing_pairs24, diag = FALSE)]
cor_values_before24 <- cor_before_24[upper.tri(cor_before_24, diag = FALSE)]
cor_values_after24 <- cor_after_24[upper.tri(cor_after_24, diag = FALSE)]


dev.off()
#2023
plot(cor_values_before23, cor_values_after23,
     xlab = "Before Imputation",
     ylab = "After Imputation",
     main = "Correlation Comparison",
     col = ifelse(imputed_indices23, "red", "blue"), #red for imputed pairs
     pch = 19) 

abline(0,1,col = "black", lwd = 2) #ideally points should line on line

#2024
plot(cor_values_before24, cor_values_after24,
     xlab = "Before Imputation",
     ylab = "After Imputation",
     main = "Correlation Comparison",
     col = ifelse(imputed_indices24, "red", "blue"), #red for imputed pairs
     pch = 19) 

abline(0,1,col = "black", lwd = 2) #ideally points should line on line

