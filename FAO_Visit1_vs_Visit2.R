#load necessary libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(missForest)
library(imputeLCMD)
library(reshape2) #for melt
library(FSA) #for Dunns test
library(gtools) #for mixedsort
library(forcats) #for factor ordering
library(broom)
library(patchwork)


#load original data
FAO_original <- read.csv("/Users/marcinebessire/Desktop/project/FAO_data.csv", check.names = FALSE)

#load data with missing values 
FAO_10pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_10pct.csv", check.names = FALSE)
FAO_20pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_20pct.csv", check.names = FALSE)
FAO_30pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_30pct.csv", check.names = FALSE)
FAO_40pct <- read.csv("/Users/marcinebessire/Desktop/project/FAO_40pct.csv", check.names = FALSE)

#separate data by Visit
FAO_10pct_1 <- FAO_10pct %>% filter(Visit == "Visit 1")
FAO_10pct_2 <- FAO_10pct %>% filter(Visit == "Visit 2")
FAO_20pct_1 <- FAO_20pct %>% filter(Visit == "Visit 1")
FAO_20pct_2 <- FAO_20pct %>% filter(Visit == "Visit 2")
FAO_30pct_1 <- FAO_30pct %>% filter(Visit == "Visit 1")
FAO_30pct_2 <- FAO_30pct %>% filter(Visit == "Visit 2")
FAO_40pct_1 <- FAO_40pct %>% filter(Visit == "Visit 1")
FAO_40pct_2 <- FAO_40pct %>% filter(Visit == "Visit 2")

# ------------------------------------
# Mean, Median and CV of FAO original
# ------------------------------------

#separate data by Visit
visit1_data <- FAO_original %>% filter(Visit == "Visit 1")
visit2_data <- FAO_original %>% filter(Visit == "Visit 2")

#numeric data
FAO_numeric_visit1 <- visit1_data[, 6:ncol(visit1_data)]
FAO_numeric_visit2 <- visit2_data[, 6:ncol(visit2_data)]

#Visit 1
#initialize empty dataset
summary_stats1 <- data.frame(Metabolite = colnames(FAO_numeric_visit1),
                            Mean = NA,
                            Median = NA, 
                            CV = NA)

#loop through each column and compute statistics
for (i in seq_along(FAO_numeric_visit1)) {
  col_data <- FAO_numeric_visit1[[i]]
  
  #compte mean, median and CV 
  mean_val <- mean(col_data, na.rm = TRUE)
  median_val <- median(col_data, na.rm = TRUE)
  cv_val <- (sd(col_data, na.rm = TRUE) / mean_val) * 100 #in percent
  
  #store in dataframe
  summary_stats1[i, "Mean"] <- mean_val
  summary_stats1[i, "Median"] <- median_val
  summary_stats1[i, "CV"] <- cv_val
}

#change to long format for plotting 
summary_long1 <- melt(summary_stats1, id.vars = "Metabolite", measure.vars = c("Mean", "Median"))

#bar plot mean median
ggplot(summary_long1, aes(x = fct_reorder(Metabolite, Metabolite), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Mean and Median of Metabolites",
       x = "Metabolite", y = "Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("steelblue", "orange"))

#bar plot CV per metabolite
ggplot(summary_stats1, aes(x = fct_reorder(Metabolite, -CV), y = CV)) +
  geom_bar(stat = "identity", fill = "magenta", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Coefficient of Variation (CV) per Metabolite",
       x = "Metabolite", y = "CV (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Visit 2
#initialize empty dataset
summary_stats2 <- data.frame(Metabolite = colnames(FAO_numeric_visit2),
                             Mean = NA,
                             Median = NA, 
                             CV = NA)

#loop through each column and compute statistics
for (i in seq_along(FAO_numeric_visit2)) {
  col_data <- FAO_numeric_visit2[[i]]
  
  #compte mean, median and CV 
  mean_val <- mean(col_data, na.rm = TRUE)
  median_val <- median(col_data, na.rm = TRUE)
  cv_val <- (sd(col_data, na.rm = TRUE) / mean_val) * 100 #in percent
  
  #store in dataframe
  summary_stats2[i, "Mean"] <- mean_val
  summary_stats2[i, "Median"] <- median_val
  summary_stats2[i, "CV"] <- cv_val
}

#change to long format for plotting 
summary_long2 <- melt(summary_stats2, id.vars = "Metabolite", measure.vars = c("Mean", "Median"))

#bar plot mean median
ggplot(summary_long2, aes(x = fct_reorder(Metabolite, Metabolite), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Mean and Median of Metabolites",
       x = "Metabolite", y = "Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("steelblue", "orange"))

#bar plot CV per metabolite
ggplot(summary_stats2, aes(x = fct_reorder(Metabolite, -CV), y = CV)) +
  geom_bar(stat = "identity", fill = "magenta", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Coefficient of Variation (CV) per Metabolite",
       x = "Metabolite", y = "CV (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ------------------------------------
# Part 0: Visit 1 vs Visit 2 Plot 
# ------------------------------------

#reshape data into long format (use melt to convert dataframe from wide into long format)
FAO_long <- melt(FAO_original, id.vars = c("ID", "Participant", "MonthDay", "Year", "Visit"),
                 variable.name = "Metabolite", value.name = "Value")

#generate boxplot for each metabolite 
ggplot(FAO_long, aes(x = Visit, y = Value, fill = Visit)) +
  geom_boxplot() +
  facet_wrap(~Metabolite, scales = "free") +
  theme_minimal(base_size = 14) +  # Adjusts the base font size
  labs(title = "Comparison of Metabolites Between Visits", 
       x = "Visit",
       y = "Metabolite Measurement") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # Increase x-axis text size
        axis.text.y = element_text(size = 16),  # Increase y-axis text size
        axis.title = element_text(size = 16),  # Increase axis title size
        strip.text = element_text(size = 16),  # Increase facet label size
        plot.title = element_text(size = 18, face = "bold"))  # Increase title size


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
#Visit1 and 2
Halfmin_10pct_1 <- half_min_imputation(FAO_10pct_1)
Halfmin_10pct_2 <- half_min_imputation(FAO_10pct_2)
Halfmin_20pct_1 <- half_min_imputation(FAO_20pct_1)
Halfmin_20pct_2 <- half_min_imputation(FAO_20pct_2)
Halfmin_30pct_1 <- half_min_imputation(FAO_30pct_1)
Halfmin_30pct_2 <- half_min_imputation(FAO_30pct_2)
Halfmin_40pct_1 <- half_min_imputation(FAO_40pct_1)
Halfmin_40pct_2 <- half_min_imputation(FAO_40pct_2)

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
#Visit 1
KNN_10pct_1 <- KNN_imputation(FAO_10pct_1)
KNN_20pct_1 <- KNN_imputation(FAO_20pct_1)
KNN_30pct_1 <- KNN_imputation(FAO_30pct_1)
KNN_40pct_1 <- KNN_imputation(FAO_40pct_1)
#Visit 2
KNN_10pct_2 <- KNN_imputation(FAO_10pct_2)
KNN_20pct_2 <- KNN_imputation(FAO_20pct_2)
KNN_30pct_2 <- KNN_imputation(FAO_30pct_2)
KNN_40pct_2 <- KNN_imputation(FAO_40pct_2)

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
#Visit 1
RF_10pct_1 <- RF_imputation(FAO_10pct_1)
RF_20pct_1 <- RF_imputation(FAO_20pct_1)
RF_30pct_1 <- RF_imputation(FAO_30pct_1)
RF_40pct_1 <- RF_imputation(FAO_40pct_1)
#Visit 2
RF_10pct_2 <- RF_imputation(FAO_10pct_2)
RF_20pct_2 <- RF_imputation(FAO_20pct_2)
RF_30pct_2 <- RF_imputation(FAO_30pct_2)
RF_40pct_2 <- RF_imputation(FAO_40pct_2)

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
#Visit 1
QRILC_10pct_1 <- QRILC_impuation(FAO_10pct_1)
QRILC_20pct_1 <- QRILC_impuation(FAO_20pct_1)
QRILC_30pct_1 <- QRILC_impuation(FAO_30pct_1)
QRILC_40pct_1 <- QRILC_impuation(FAO_40pct_1)
#Visit 2
QRILC_10pct_2 <- QRILC_impuation(FAO_10pct_2)
QRILC_20pct_2 <- QRILC_impuation(FAO_20pct_2)
QRILC_30pct_2 <- QRILC_impuation(FAO_30pct_2)
QRILC_40pct_2 <- QRILC_impuation(FAO_40pct_2)

# ---------------------------------
# Part 1.5: Combine datasets into 1
# ----------------------------------

#make funciton to merge and order dataframes for each Imputation
merge_and_order <- function(data1, data2) {
  # Combine datasets
  merged_data <- rbind(data1, data2)
  
  # Convert Participant to a properly sorted factor
  merged_data$Participant <- factor(merged_data$Participant,
                                    levels = mixedsort(unique(merged_data$Participant)))
  
  # Create a proper date column for ordering
  merged_data$FullDate <- as.Date(paste(merged_data$Year, merged_data$MonthDay, sep = "-"), format = "%Y-%m/%d")
  
  # Order by Participant and FullDate (which correctly sorts by Year first, then MonthDay)
  merged_data <- merged_data[order(merged_data$Participant, merged_data$FullDate), ]
  
  # Remove temporary FullDate column (optional)
  merged_data$FullDate <- NULL
  
  # Reset row numbering to ensure sequential indices
  row.names(merged_data) <- NULL  
  
  return(merged_data)
}

#call function
#Halfmin
Halfmin_10pct_tot <- merge_and_order(Halfmin_10pct_1, Halfmin_10pct_2)
Halfmin_20pct_tot <- merge_and_order(Halfmin_20pct_1, Halfmin_20pct_2)
Halfmin_30pct_tot <- merge_and_order(Halfmin_30pct_1, Halfmin_30pct_2)
Halfmin_40pct_tot <- merge_and_order(Halfmin_10pct_1, Halfmin_40pct_2)
#KNN
KNN_10pct_tot <- merge_and_order(KNN_10pct_1, KNN_10pct_2)
KNN_20pct_tot <- merge_and_order(KNN_20pct_1, KNN_20pct_2)
KNN_30pct_tot <- merge_and_order(KNN_30pct_1, KNN_30pct_2)
KNN_40pct_tot <- merge_and_order(KNN_40pct_1, KNN_40pct_2)
#RF
RF_10pct_tot <- merge_and_order(RF_10pct_1, RF_10pct_2)
RF_20pct_tot <- merge_and_order(RF_20pct_1, RF_20pct_2)
RF_30pct_tot <- merge_and_order(RF_30pct_1, RF_30pct_2)
RF_40pct_tot <- merge_and_order(RF_40pct_1, RF_40pct_2)
#QRILC
QRILC_10pct_tot <- merge_and_order(QRILC_10pct_1, QRILC_10pct_2)
QRILC_20pct_tot <- merge_and_order(QRILC_20pct_1, QRILC_20pct_2)
QRILC_30pct_tot <- merge_and_order(QRILC_30pct_1, QRILC_30pct_2)
QRILC_40pct_tot <- merge_and_order(QRILC_40pct_1, QRILC_40pct_2)

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
shapiro_original <- shapiro_test(FAO_original) #11/34 are significant
shapiro_original1 <- shapiro_test(visit1_data) #2/34
shapiro_original2 <- shapiro_test(visit2_data) #1/34

#knn
#visit 1
shapiro1_KNN10 <- shapiro_test(KNN_10pct_1) #8 are significant
shapiro1_KNN20 <- shapiro_test(KNN_20pct_1) #9 are significant
shapiro1_KNN30 <- shapiro_test(KNN_30pct_1) #7 are significant
shapiro1_KNN40 <- shapiro_test(KNN_40pct_1) #7 are significant
#visit 2
shapiro2_KNN10 <- shapiro_test(KNN_10pct_2) #6 are significant
shapiro2_KNN20 <- shapiro_test(KNN_20pct_2) #8 are significant
shapiro2_KNN30 <- shapiro_test(KNN_30pct_2) #10 are significant
shapiro2_KNN40 <- shapiro_test(KNN_40pct_2) #3 are significant

#Half-min
#visit 1
shapiro1_Halfmin10 <- shapiro_test(Halfmin_10pct_1) #1 are significant
shapiro1_Halfmin20 <- shapiro_test(Halfmin_20pct_1) #2 are significant
shapiro1_Halfmin30 <- shapiro_test(Halfmin_30pct_1) #1 are significant
shapiro1_Halfmin40 <- shapiro_test(Halfmin_40pct_1) #0 are significant
#visit 2
shapiro2_Halfmin10 <- shapiro_test(Halfmin_10pct_2) #1 are significant
shapiro2_Halfmin20 <- shapiro_test(Halfmin_20pct_2) #0 are significant
shapiro2_Halfmin30 <- shapiro_test(Halfmin_30pct_2) #1 are significant
shapiro2_Halfmin40 <- shapiro_test(Halfmin_40pct_2) #3 are significant

#RF
#visit 1
shaprio1_RF10 <- shapiro_test(RF_10pct_1) #2 are significant 
shaprio1_RF20 <- shapiro_test(RF_20pct_1) #1 are significant 
shaprio1_RF30 <- shapiro_test(RF_30pct_1) #2 are significant 
shaprio1_RF40 <- shapiro_test(RF_40pct_1) #1 are significant 
#visit 2
shaprio2_RF10 <- shapiro_test(RF_10pct_2) #1 are significant 
shaprio2_RF20 <- shapiro_test(RF_20pct_2) #2 are significant 
shaprio2_RF30 <- shapiro_test(RF_30pct_2) #4 are significant 
shaprio2_RF40 <- shapiro_test(RF_40pct_2) #2 are significant

#QRILC
#visit1
shapiro1_QRILC10 <- shapiro_test(QRILC_10pct_1) #1 are significant
shapiro1_QRILC20 <- shapiro_test(QRILC_20pct_1) #1 are significant
shapiro1_QRILC30 <- shapiro_test(QRILC_30pct_1) #1 are significant
shapiro1_QRILC40 <- shapiro_test(QRILC_40pct_1) #0 are significant
#visit2
shapiro2_QRILC10 <- shapiro_test(QRILC_10pct_2) #1 are significant
shapiro2_QRILC20 <- shapiro_test(QRILC_20pct_2) #2 are significant
shapiro2_QRILC30 <- shapiro_test(QRILC_30pct_2) #5 are significant
shapiro2_QRILC40 <- shapiro_test(QRILC_40pct_2) #3 are significant

#Visit 1
#create dataframe with results
shapiro_summary1 <- data.frame(
  Method = c(
    rep("Original", 1),
    rep("KNN", 4),
    rep("Half-min", 4),
    rep("RF", 4),
    rep("QRILC", 4)
  ),
  Missingness = c(
    0, # Original has only one value
    
    10, 20, 30, 40,  # KNN
    10, 20, 30, 40,  # Half-min
    10, 20, 30, 40,  # RF
    10, 20, 30, 40   # QRILC
  ),
  Non_Normal_Count = c(
    2, # Original dataset
    
    8, 9, 7, 7,  # KNN
    1, 2, 1, 0,     # Half-min
    2, 1, 2, 1,     # RF
    1, 1, 1, 0     # QRILC
  )
)

ggplot(shapiro_summary1, aes(x = factor(Missingness), y = Non_Normal_Count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Number of Non-Normally Distributed Metabolites (Visit 1)",
       x = "Missingness Percentage (%)",
       y = "Count of Non-Normal Metabolites",
       fill = "Imputation Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#extract the Non_Normal_Count for the "Original" method
original_non_normal1 <- shapiro_summary1 %>%
  filter(Method == "Original") %>%
  pull(Non_Normal_Count)

ggplot(shapiro_summary1, aes(x = Missingness, y = Non_Normal_Count, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = original_non_normal1, linetype = "dashed", color = "lightgreen", size = 1) + # Dashed reference line
  theme_minimal(base_size = 16) +  # Increase overall font size
  labs(
    title = "Effect of Imputation on Normality of Metabolites (Visit 1)",
    x = "Missingness Percentage (%)",
    y = "Count of Non-Normal Metabolites",
    color = "Imputation Method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),  # Improve x-axis readability
    axis.text.y = element_text(size = 14, face = "bold"),  # Improve y-axis readability
    axis.title.x = element_text(size = 16, face = "bold"),  # Make x-axis title bigger and bold
    axis.title.y = element_text(size = 16, face = "bold"),  # Make y-axis title bigger and bold
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Center and enlarge title
    legend.title = element_text(size = 16, face = "bold"),  # Increase legend title size
    legend.text = element_text(size = 14)  # Increase legend text size
  )


#Visit 2
#create dataframe with results
shapiro_summary2 <- data.frame(
  Method = c(
    rep("Original", 1),
    rep("KNN", 4),
    rep("Half-min", 4),
    rep("RF", 4),
    rep("QRILC", 4)
  ),
  Missingness = c(
    0, # Original has only one value
    
    10, 20, 30, 40,  # KNN
    10, 20, 30, 40,  # Half-min
    10, 20, 30, 40,  # RF
    10, 20, 30, 40   # QRILC
  ),
  Non_Normal_Count = c(
    1, # Original dataset
    
    6, 8, 10, 3,  # KNN
    1, 0, 1, 3,     # Half-min
    1, 2, 4, 2,     # RF
    1, 2, 5, 3     # QRILC
  )
)

ggplot(shapiro_summary2, aes(x = factor(Missingness), y = Non_Normal_Count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Number of Non-Normally Distributed Metabolites (Visit 2)",
       x = "Missingness Percentage (%)",
       y = "Count of Non-Normal Metabolites",
       fill = "Imputation Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#extract the Non_Normal_Count for the "Original" method
original_non_normal2 <- shapiro_summary2 %>%
  filter(Method == "Original") %>%
  pull(Non_Normal_Count)

ggplot(shapiro_summary2, aes(x = Missingness, y = Non_Normal_Count, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = original_non_normal2, linetype = "dashed", color = "lightgreen", size = 1) + # Dashed reference line
  theme_minimal(base_size = 16) +  # Increase overall font size
  labs(
    title = "Effect of Imputation on Normality of Metabolites (Visit 2)",
    x = "Missingness Percentage (%)",
    y = "Count of Non-Normal Metabolites",
    color = "Imputation Method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),  # Improve x-axis readability
    axis.text.y = element_text(size = 14, face = "bold"),  # Improve y-axis readability
    axis.title.x = element_text(size = 16, face = "bold"),  # Make x-axis title bigger and bold
    axis.title.y = element_text(size = 16, face = "bold"),  # Make y-axis title bigger and bold
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Center and enlarge title
    legend.title = element_text(size = 16, face = "bold"),  # Increase legend title size
    legend.text = element_text(size = 14)  # Increase legend text size
  )



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
t_test_half_min_10pct <- t_test_func(FAO_original, Halfmin_10pct_tot) #0/34
t_test_half_min_20pct <- t_test_func(FAO_original, Halfmin_20pct_tot) #0/34
t_test_half_min_30pct <- t_test_func(FAO_original, Halfmin_30pct_tot) #0/34
t_test_half_min_40pct <- t_test_func(FAO_original, Halfmin_40pct_tot) #0/34
#KNN
t_test_KNN_10pct <- t_test_func(FAO_original, KNN_10pct_tot) #0/34
t_test_KNN_20pct <- t_test_func(FAO_original, KNN_20pct_tot) #0/34
t_test_KNN_30pct <- t_test_func(FAO_original, KNN_30pct_tot) #0/34
t_test_KNN_40pct <- t_test_func(FAO_original, KNN_40pct_tot) #10/34
#RF
t_test_RF_10pct <- t_test_func(FAO_original, RF_10pct_tot) #0/34
t_test_RF_20pct <- t_test_func(FAO_original, RF_20pct_tot) #0/34
t_test_RF_30pct <- t_test_func(FAO_original, RF_30pct_tot) #0/34
t_test_RF_40pct <- t_test_func(FAO_original, RF_40pct_tot) #0/34
#QRILC
t_test_QRILC_10pct <- t_test_func(FAO_original, QRILC_10pct_tot) #0/34
t_test_QRILC_20pct <- t_test_func(FAO_original, QRILC_20pct_tot) #0/34
t_test_QRILC_30pct <- t_test_func(FAO_original, QRILC_30pct_tot) #0/34
t_test_QRILC_40pct <- t_test_func(FAO_original, QRILC_40pct_tot) #0/34


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
wilcox_half_min_10pct <- wilcoxon_func(FAO_original, Halfmin_10pct_tot) #0/34
wilcox_half_min_20pct <- wilcoxon_func(FAO_original, Halfmin_20pct_tot) #0/34
wilcox_half_min_30pct <- wilcoxon_func(FAO_original, Halfmin_30pct_tot) #0/34
wilcox_half_min_40pct <- wilcoxon_func(FAO_original, Halfmin_40pct_tot) #0/34
#KNN
wilcox_KNN_10pct <- wilcoxon_func(FAO_original, KNN_10pct_tot) #0/34
wilcox_KNN_20pct <- wilcoxon_func(FAO_original, KNN_20pct_tot) #0/34
wilcox_KNN_30pct <- wilcoxon_func(FAO_original, KNN_30pct_tot) #0/34
wilcox_KNN_40pct <- wilcoxon_func(FAO_original, KNN_40pct_tot) #0/34
#RF
wilcox_RF_10pct <- wilcoxon_func(FAO_original, RF_10pct_tot) #0/34
wilcox_RF_20pct <- wilcoxon_func(FAO_original, RF_20pct_tot) #0/34
wilcox_RF_30pct <- wilcoxon_func(FAO_original, RF_30pct_tot) #0/34
wilcox_RF_40pct <- wilcoxon_func(FAO_original, RF_40pct_tot) #0/34
#QRILC
wilcox_QRILC_10pct <- wilcoxon_func(FAO_original, QRILC_10pct_tot) #0/34
wilcox_QRILC_20pct <- wilcoxon_func(FAO_original, QRILC_20pct_tot) #0/34
wilcox_QRILC_30pct <- wilcoxon_func(FAO_original, QRILC_30pct_tot) #0/34
wilcox_QRILC_40pct <- wilcoxon_func(FAO_original, QRILC_40pct_tot) #0/34

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
visit_Halfmin10pct_res <- visit_statistical_tests(Halfmin_10pct_tot) #0/34 W and 0/34 T
visit_Halfmin20pct_res <- visit_statistical_tests(Halfmin_20pct_tot) #0/34 W and 3/34 T
visit_Halfmin30pct_res <- visit_statistical_tests(Halfmin_30pct_tot) #0/34 W and 4/34 T
visit_Halfmin40pct_res <- visit_statistical_tests(Halfmin_40pct_tot) #0/34 W and 1/34 T
#KNN
visit_KNN10pct_res <- visit_statistical_tests(KNN_10pct_tot) #13/34 W 13/34 T
visit_KNN20pct_res <- visit_statistical_tests(KNN_20pct_tot) #12/34 W 11/34 T
visit_KNN30pct_res <- visit_statistical_tests(KNN_30pct_tot) #10/11 W 10/34 T
visit_KNN40pct_res <- visit_statistical_tests(KNN_40pct_tot) #6/34 W 7/34 T
#RF
visit_RF10pct_res <- visit_statistical_tests(RF_10pct_tot) #17/34 W 16/34 T
visit_RF20pct_res <- visit_statistical_tests(RF_20pct_tot) #21/34 W 21/34 T
visit_RF30pct_res <- visit_statistical_tests(RF_30pct_tot) #18/34 W 19/34 T
visit_RF40pct_res <- visit_statistical_tests(RF_40pct_tot) #21/34 W 20/34 T
#QRILC
visit_QRILC10pct_res <- visit_statistical_tests(QRILC_10pct_tot) #13/34 W 15/34 T
visit_QRILC20pct_res <- visit_statistical_tests(QRILC_20pct_tot) #13/34 W 12/34 T
visit_QRILC30pct_res <- visit_statistical_tests(QRILC_30pct_tot) #10/34 W 10/34 T
visit_QRILC40pct_res <- visit_statistical_tests(QRILC_40pct_tot) #0/34 W 4/34 T

#Create data frame for number of significant metabolites
data <- data.frame(
  Method = rep(c("Halfmin", "Halfmin", "Halfmin", "Halfmin", "Halfmin",
                 "KNN", "KNN", "KNN", "KNN", "KNN",
                 "RF", "RF", "RF", "RF", "RF",
                 "QRILC", "QRILC", "QRILC", "QRILC", "QRILC"), each = 1),
  Percentage = c("Original", "10%", "20%", "30%", "40%",
                 "Original", "10%", "20%", "30%", "40%",
                 "Original", "10%", "20%", "30%", "40%",
                 "Original", "10%", "20%", "30%", "40%"),
  Wilcoxon = c(18, 0, 0, 0, 0,   
               18, 13, 12, 10, 6, 
               18, 17, 21, 18, 21,  
               18, 13, 13, 10, 0),  
  TTest = c(16, 0, 3, 4, 1,  
            16, 13, 11, 10, 7,  
            16, 16, 21, 19, 20,  
            16, 15, 12, 10, 4)  
)

# Convert Percentage to factor to maintain ordering
data$Percentage <- factor(data$Percentage, levels = c("Original", "10%", "20%", "30%", "40%"))

# Reshape data for plotting
data_long <- data %>% 
  pivot_longer(cols = c("Wilcoxon", "TTest"), names_to = "Test", values_to = "Significant_Metabolites")

#ceate grouped bar plot ensuring "Original" appears first in each facet
ggplot(data_long, aes(x = Percentage, y = Significant_Metabolites, fill = Test)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Significant_Metabolites), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 6, fontface = "bold") +  
  facet_wrap(~ Method, scales = "free_x", nrow = 2) +  
  labs(title = "Number of Significant Metabolites (Visit 1 vs Visit 2)",
       x = "Missingness Percentage",
       y = "Number of Significant Metabolites",
       fill = "Statistical Test") +
  theme_minimal(base_size = 16) +  # Increase overall font size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),  # Larger x-axis labels
    axis.text.y = element_text(size = 14, face = "bold"),  # Larger y-axis labels
    axis.title.x = element_text(size = 16, face = "bold"),  # Larger and bold x-axis title
    axis.title.y = element_text(size = 16, face = "bold"),  # Larger and bold y-axis title
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Centered and bold title
    legend.title = element_text(size = 16, face = "bold"),  # Larger and bold legend title
    legend.text = element_text(size = 14),  # Larger legend text
    strip.text = element_text(size = 16, face = "bold"),  # Larger facet label size
    panel.spacing = unit(2, "lines")  # More spacing between facets
  )


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
#visit 1
nrmse_res1_halfmin_10pct <- calculate_weighted_nrmse(visit1_data, Halfmin_10pct_1, "Halfmin", 0.1)
nrmse_res1_halfmin_20pct <- calculate_weighted_nrmse(visit1_data, Halfmin_20pct_1,"Halfmin", 0.2)
nrmse_res1_halfmin_30pct <- calculate_weighted_nrmse(visit1_data, Halfmin_30pct_1, "Halfmin", 0.3)
nrmse_res1_halfmin_40pct <- calculate_weighted_nrmse(visit1_data, Halfmin_40pct_1, "Halfmin", 0.4)
#visit 2
nrmse_res2_halfmin_10pct <- calculate_weighted_nrmse(visit2_data, Halfmin_10pct_2, "Halfmin", 0.1)
nrmse_res2_halfmin_20pct <- calculate_weighted_nrmse(visit2_data, Halfmin_20pct_2,"Halfmin", 0.2)
nrmse_res2_halfmin_30pct <- calculate_weighted_nrmse(visit2_data, Halfmin_30pct_2, "Halfmin", 0.3)
nrmse_res2_halfmin_40pct <- calculate_weighted_nrmse(visit2_data, Halfmin_40pct_2, "Halfmin", 0.4)
#KNN
#visit 1
nrmse_res1_KNN_10pct <- calculate_weighted_nrmse(visit1_data, KNN_10pct_1, "KNN", 0.1)
nrmse_res1_KNN_20pct <- calculate_weighted_nrmse(visit1_data, KNN_20pct_1, "KNN", 0.2)
nrmse_res1_KNN_30pct <- calculate_weighted_nrmse(visit1_data, KNN_30pct_1, "KNN", 0.3)
nrmse_res1_KNN_40pct <- calculate_weighted_nrmse(visit1_data, KNN_40pct_1, "KNN", 0.4)
#visit 2
nrmse_res2_KNN_10pct <- calculate_weighted_nrmse(visit2_data, KNN_10pct_2, "KNN", 0.1)
nrmse_res2_KNN_20pct <- calculate_weighted_nrmse(visit2_data, KNN_20pct_2, "KNN", 0.2)
nrmse_res2_KNN_30pct <- calculate_weighted_nrmse(visit2_data, KNN_30pct_2, "KNN", 0.3)
nrmse_res2_KNN_40pct <- calculate_weighted_nrmse(visit2_data, KNN_40pct_2, "KNN", 0.4)
#RF
#visit 1
nrmse_res1_RF_10pct <- calculate_weighted_nrmse(visit1_data, RF_10pct_1, "RF", 0.1)
nrmse_res1_RF_20pct <- calculate_weighted_nrmse(visit1_data, RF_20pct_1, "RF", 0.2)
nrmse_res1_RF_30pct <- calculate_weighted_nrmse(visit1_data, RF_30pct_1, "RF", 0.3)
nrmse_res1_RF_40pct <- calculate_weighted_nrmse(visit1_data, RF_40pct_1, "RF", 0.4)
#visit 2
nrmse_res2_RF_10pct <- calculate_weighted_nrmse(visit2_data, RF_10pct_2, "RF", 0.1)
nrmse_res2_RF_20pct <- calculate_weighted_nrmse(visit2_data, RF_20pct_2, "RF", 0.2)
nrmse_res2_RF_30pct <- calculate_weighted_nrmse(visit2_data, RF_30pct_2, "RF", 0.3)
nrmse_res2_RF_40pct <- calculate_weighted_nrmse(visit2_data, RF_40pct_2, "RF", 0.4)
#QRILC
#visit 1
nrmse_res1_QRILC_10pct <- calculate_weighted_nrmse(visit1_data, QRILC_10pct_1, "QRILC", 0.1)
nrmse_res1_QRILC_20pct <- calculate_weighted_nrmse(visit1_data, QRILC_20pct_1, "QRILC", 0.2)
nrmse_res1_QRILC_30pct <- calculate_weighted_nrmse(visit1_data, QRILC_30pct_1, "QRILC", 0.3)
nrmse_res1_QRILC_40pct <- calculate_weighted_nrmse(visit1_data, QRILC_40pct_1, "QRILC", 0.4)
#visit 2
nrmse_res2_QRILC_10pct <- calculate_weighted_nrmse(visit2_data, QRILC_10pct_2, "QRILC", 0.1)
nrmse_res2_QRILC_20pct <- calculate_weighted_nrmse(visit2_data, QRILC_20pct_2, "QRILC", 0.2)
nrmse_res2_QRILC_30pct <- calculate_weighted_nrmse(visit2_data, QRILC_30pct_2, "QRILC", 0.3)
nrmse_res2_QRILC_40pct <- calculate_weighted_nrmse(visit2_data, QRILC_40pct_2, "QRILC", 0.4)

# ------------------------------------
# Part 3.1: NRMSE Plot 
# ------------------------------------

#for plotting combine all nrmse results into one dataframe
#visit 1
nrmse_data1 <- bind_rows(
  nrmse_res1_halfmin_10pct, nrmse_res1_halfmin_20pct, 
  nrmse_res1_halfmin_30pct, nrmse_res1_halfmin_40pct,
  
  nrmse_res1_KNN_10pct, nrmse_res1_KNN_20pct, 
  nrmse_res1_KNN_30pct, nrmse_res1_KNN_40pct,
  
  nrmse_res1_RF_10pct, nrmse_res1_RF_20pct, 
  nrmse_res1_RF_30pct, nrmse_res1_RF_40pct,
  
  nrmse_res1_QRILC_10pct, nrmse_res1_QRILC_20pct, 
  nrmse_res1_QRILC_30pct, nrmse_res1_QRILC_40pct
)

#convert MCAR_Proportion to a factor for correct ordering in the plot
nrmse_data1$MCAR_Proportion <- as.factor(nrmse_data1$MCAR_Proportion)

#Generate Boxplot 
ggplot(nrmse_data1, aes(x = MCAR_Proportion, y = Weighted_NRMSE, fill = Imputation_Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + #boxplot wiht transparency and w/o outliers
  scale_fill_manual(values = c("lightblue", "orange", "blue", "magenta")) +
  labs(
    title = "Weighted NRMSE across Imputation Method and MCAR Proportions (Visit 1)",
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


#visit 2
nrmse_data2 <- bind_rows(
  nrmse_res2_halfmin_10pct, nrmse_res2_halfmin_20pct, 
  nrmse_res2_halfmin_30pct, nrmse_res2_halfmin_40pct,
  
  nrmse_res2_KNN_10pct, nrmse_res2_KNN_20pct, 
  nrmse_res2_KNN_30pct, nrmse_res2_KNN_40pct,
  
  nrmse_res2_RF_10pct, nrmse_res2_RF_20pct, 
  nrmse_res2_RF_30pct, nrmse_res2_RF_40pct,
  
  nrmse_res2_QRILC_10pct, nrmse_res2_QRILC_20pct, 
  nrmse_res2_QRILC_30pct, nrmse_res2_QRILC_40pct
)

#convert MCAR_Proportion to a factor for correct ordering in the plot
nrmse_data2$MCAR_Proportion <- as.factor(nrmse_data2$MCAR_Proportion)

#Generate Boxplot 
ggplot(nrmse_data2, aes(x = MCAR_Proportion, y = Weighted_NRMSE, fill = Imputation_Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + #boxplot wiht transparency and w/o outliers
  scale_fill_manual(values = c("lightblue", "orange", "blue", "magenta")) +
  labs(
    title = "Weighted NRMSE across Imputation Method and MCAR Proportions (Visit 2)",
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
norm_mean_diff <- function(original, imputed, method, percentage, visit) {
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
  plot_title <- paste0("Normalized Difference with ", percentage, "% Missing Values using ", method, " Imputation ", "(Visit ", visit, ")")
  
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
#visit 1
norm_diff1_Halfmin_10pct <- norm_mean_diff(visit1_data, Halfmin_10pct_1, "Half-min", 10, 1)
norm_diff1_Halfmin_20pct <- norm_mean_diff(visit1_data, Halfmin_20pct_1, "Half-min", 20, 1)
norm_diff1_Halfmin_30pct <- norm_mean_diff(visit1_data, Halfmin_30pct_1, "Half-min", 30, 1)
norm_diff1_Halfmin_40pct <- norm_mean_diff(visit1_data, Halfmin_40pct_1, "Half-min", 40, 1)
#visit 2
norm_diff2_Halfmin_10pct <- norm_mean_diff(visit2_data, Halfmin_10pct_2, "Half-min", 10, 2)
norm_diff2_Halfmin_20pct <- norm_mean_diff(visit2_data, Halfmin_20pct_2, "Half-min", 20, 2)
norm_diff2_Halfmin_30pct <- norm_mean_diff(visit2_data, Halfmin_30pct_2, "Half-min", 30, 2)
norm_diff2_Halfmin_40pct <- norm_mean_diff(visit2_data, Halfmin_40pct_2, "Half-min", 40, 2)

#KNN
#Visit 1
norm_diff1_KNN_10pct <- norm_mean_diff(visit1_data, KNN_10pct_1, "KNN", 10, 1)
norm_diff1_KNN_20pct <- norm_mean_diff(visit1_data, KNN_20pct_1, "KNN", 20, 1)
norm_diff1_KNN_30pct <- norm_mean_diff(visit1_data, KNN_30pct_1, "KNN", 30, 1)
norm_diff1_KNN_40pct <- norm_mean_diff(visit1_data, KNN_40pct_1, "KNN", 40, 1)
#Visit 2
norm_diff2_KNN_10pct <- norm_mean_diff(visit2_data, KNN_10pct_2, "KNN", 10, 2)
norm_diff2_KNN_20pct <- norm_mean_diff(visit2_data, KNN_20pct_2, "KNN", 20, 2)
norm_diff2_KNN_30pct <- norm_mean_diff(visit2_data, KNN_30pct_2, "KNN", 30, 2)
norm_diff2_KNN_40pct <- norm_mean_diff(visit2_data, KNN_40pct_2, "KNN", 40, 2)

#RF
#Visit 1
norm_diff1_RF_10pct <- norm_mean_diff(visit1_data, RF_10pct_1, "RF", 10, 1)
norm_diff1_RF_20pct <- norm_mean_diff(visit1_data, RF_20pct_1, "RF", 20, 1)
norm_diff1_RF_30pct <- norm_mean_diff(visit1_data, RF_30pct_1, "RF", 30, 1)
norm_diff1_RF_40pct <- norm_mean_diff(visit1_data, RF_40pct_1, "RF", 40, 1)
#Visit 2
norm_diff2_RF_10pct <- norm_mean_diff(visit2_data, RF_10pct_2, "RF", 10, 2)
norm_diff2_RF_20pct <- norm_mean_diff(visit2_data, RF_20pct_2, "RF", 20, 2)
norm_diff2_RF_30pct <- norm_mean_diff(visit2_data, RF_30pct_2, "RF", 30, 2)
norm_diff2_RF_40pct <- norm_mean_diff(visit2_data, RF_40pct_2, "RF", 40, 2)

#QRILC
#Visit 1
norm_diff1_QRILC_10pct <- norm_mean_diff(visit1_data, QRILC_10pct_1, "QRILC", 10, 1)
norm_diff1_QRILC_20pct <- norm_mean_diff(visit1_data, QRILC_20pct_1, "QRILC", 20, 1)
norm_diff1_QRILC_30pct <- norm_mean_diff(visit1_data, QRILC_30pct_1, "QRILC", 30, 1)
norm_diff1_QRILC_40pct <- norm_mean_diff(visit1_data, QRILC_40pct_1, "QRILC", 40, 1)
#Visit 2
norm_diff2_QRILC_10pct <- norm_mean_diff(visit2_data, QRILC_10pct_2, "QRILC", 10, 2)
norm_diff2_QRILC_20pct <- norm_mean_diff(visit2_data, QRILC_20pct_2, "QRILC", 20, 2)
norm_diff2_QRILC_30pct <- norm_mean_diff(visit2_data, QRILC_30pct_2, "QRILC", 30, 2)
norm_diff2_QRILC_40pct <- norm_mean_diff(visit2_data, QRILC_40pct_2, "QRILC", 40, 2)

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
#visit 1
halfmin_data1 <- bind_rows(
  norm_mean_diff_data(visit1_data, Halfmin_10pct_1, "Half-min", 10),
  norm_mean_diff_data(visit1_data, Halfmin_20pct_1, "Half-min", 20),
  norm_mean_diff_data(visit1_data, Halfmin_30pct_1, "Half-min", 30),
  norm_mean_diff_data(visit1_data, Halfmin_40pct_1, "Half-min", 40)
)

knn_data1 <- bind_rows(
  norm_mean_diff_data(visit1_data, KNN_10pct_1, "KNN", 10),
  norm_mean_diff_data(visit1_data, KNN_20pct_1, "KNN", 20),
  norm_mean_diff_data(visit1_data, KNN_30pct_1, "KNN", 30),
  norm_mean_diff_data(visit1_data, KNN_40pct_1, "KNN", 40)
)

rf_data1 <- bind_rows(
  norm_mean_diff_data(visit1_data, RF_10pct_1, "RF", 10),
  norm_mean_diff_data(visit1_data, RF_20pct_1, "RF", 20),
  norm_mean_diff_data(visit1_data, RF_30pct_1, "RF", 30),
  norm_mean_diff_data(visit1_data, RF_40pct_1, "RF", 40)
)

qrilc_data1 <- bind_rows(
  norm_mean_diff_data(visit1_data, QRILC_10pct_1, "QRILC", 10),
  norm_mean_diff_data(visit1_data, QRILC_20pct_1, "QRILC", 20),
  norm_mean_diff_data(visit1_data, QRILC_30pct_1, "QRILC", 30),
  norm_mean_diff_data(visit1_data, QRILC_40pct_1, "QRILC", 40)
)

#function to plot density for each method
plot_density <- function(data, method, visit) {
  ggplot(data, aes(x = Normalized_Difference, fill = Percentage, color = Percentage)) +
    geom_density(alpha = 0.3) +
    theme_minimal() +
    labs(title = paste("Normalized Difference for", method, "Imputation (Visit ", visit, ")"),
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.2, 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    theme(legend.title = element_blank(), legend.position = "right")
}

#generate plots
plot_halfmin1 <- plot_density(halfmin_data1, "Half-min", 1)
plot_knn1 <- plot_density(knn_data1, "KNN", 1)
plot_rf1 <- plot_density(rf_data1, "RF", 1)
plot_qrilc1 <- plot_density(qrilc_data1, "QRILC", 1)

#display plots
print(plot_halfmin1)
print(plot_knn1)
print(plot_rf1)
print(plot_qrilc1)

#visit 2
halfmin_data2 <- bind_rows(
  norm_mean_diff_data(visit2_data, Halfmin_10pct_2, "Half-min", 10),
  norm_mean_diff_data(visit2_data, Halfmin_20pct_2, "Half-min", 20),
  norm_mean_diff_data(visit2_data, Halfmin_30pct_2, "Half-min", 30),
  norm_mean_diff_data(visit2_data, Halfmin_40pct_2, "Half-min", 40)
)

knn_data2 <- bind_rows(
  norm_mean_diff_data(visit2_data, KNN_10pct_2, "KNN", 10),
  norm_mean_diff_data(visit2_data, KNN_20pct_2, "KNN", 20),
  norm_mean_diff_data(visit2_data, KNN_30pct_2, "KNN", 30),
  norm_mean_diff_data(visit2_data, KNN_40pct_2, "KNN", 40)
)

rf_data2 <- bind_rows(
  norm_mean_diff_data(visit2_data, RF_10pct_2, "RF", 10),
  norm_mean_diff_data(visit2_data, RF_20pct_2, "RF", 20),
  norm_mean_diff_data(visit2_data, RF_30pct_2, "RF", 30),
  norm_mean_diff_data(visit2_data, RF_40pct_2, "RF", 40)
)

qrilc_data2 <- bind_rows(
  norm_mean_diff_data(visit2_data, QRILC_10pct_2, "QRILC", 10),
  norm_mean_diff_data(visit2_data, QRILC_20pct_2, "QRILC", 20),
  norm_mean_diff_data(visit2_data, QRILC_30pct_2, "QRILC", 30),
  norm_mean_diff_data(visit2_data, QRILC_40pct_2, "QRILC", 40)
)

#function to plot density for each method
plot_density <- function(data, method, visit) {
  ggplot(data, aes(x = Normalized_Difference, fill = Percentage, color = Percentage)) +
    geom_density(alpha = 0.3) +
    theme_minimal() +
    labs(title = paste("Normalized Difference for", method, "Imputation (Visit ", visit, ")"),
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.2, 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    theme(legend.title = element_blank(), legend.position = "right")
}

#generate plots
plot_halfmin2 <- plot_density(halfmin_data2, "Half-min", 2)
plot_knn2 <- plot_density(knn_data2, "KNN", 2)
plot_rf2 <- plot_density(rf_data2, "RF", 2)
plot_qrilc2 <- plot_density(qrilc_data2, "QRILC", 2)

#display plots
print(plot_halfmin2)
print(plot_knn2)
print(plot_rf2)
print(plot_qrilc2)

#combine and make one plot for visit 1 and one for visit 2
all_data1 <- bind_rows(halfmin_data1, knn_data1, rf_data1, qrilc_data1)
all_data2 <- bind_rows(halfmin_data2, knn_data2, rf_data2, qrilc_data2)

ggplot(all_data1, aes(x = Normalized_Difference, fill = Percentage, color = Percentage)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Normalized Difference Across Imputation Methods (Visit 1)",
       x = "Normalized Difference",
       y = "Density") +
  xlim(-0.2, 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~Method, scales = "free") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),  # Makes facet labels bold
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold")
  )


ggplot(all_data2, aes(x = Normalized_Difference, fill = Percentage, color = Percentage)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Normalized Difference Across Imputation Methods (Visit 2)",
       x = "Normalized Difference",
       y = "Density") +
  xlim(-0.2, 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~Method, scales = "free") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),  # Makes facet labels bold
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold")
  )


# ------------------------------------
# Part 5: Plot Distribution
# ------------------------------------

#function to plot the distribution before and after imputation and only imputed values 
# Function to plot the distribution before and after imputation, filtering only significant metabolites
plot_distribution <- function(original, imputed, method, percentage, visit, significant_metabolites) {
  # Select numeric metabolite columns
  numeric_original <- original[, 6:ncol(original)]
  numeric_imputed <- imputed[, 6:ncol(imputed)]
  
  # Convert to long format for plotting
  original_long <- numeric_original %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    mutate(Data = "Original")
  
  imputed_long <- numeric_imputed %>%
    pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    mutate(Data = "Imputed")
  
  # Combine data
  combined_data <- bind_rows(original_long, imputed_long)
  
  # **Filter only the significant metabolites**
  combined_data <- combined_data %>% filter(Metabolite %in% significant_metabolites)
  
  # Calculate mean values for plotting mean lines
  mean_data <- combined_data %>%
    group_by(Metabolite, Data) %>%
    summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  # Plot
  # Modify the plot to include a custom legend for the mean lines
  plot <- ggplot(combined_data, aes(x = Value, fill = Data)) +
    geom_density(alpha = 0.5, size = 1) +  # Transparency for overlapping distributions
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal(base_size = 16) +  # Increase overall font size
    labs(title = paste0("Density Distribution Before and After ", method, 
                        " Imputation (", percentage, "% Missing, Visit ", visit, ")"),
         x = "Value",
         y = "Density",
         fill = "Data") +  # Legend title for data type
    # Add mean lines for Original and Imputed data with a custom color legend
    geom_vline(data = mean_data %>% filter(Data == "Original"),
               aes(xintercept = mean_value, color = "Original Mean"), 
               linewidth = 0.8, linetype = "dashed") +
    geom_vline(data = mean_data %>% filter(Data == "Imputed"),
               aes(xintercept = mean_value, color = "Imputed Mean"), 
               linewidth = 0.8, linetype = "dashed") +
    scale_color_manual(name = "Mean Values",  # Custom legend for mean lines
                       values = c("Original Mean" = "blue", "Imputed Mean" = "red")) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 14, face = "bold"),  # Bigger and bold x-axis labels
      axis.text.y = element_text(size = 14, face = "bold"),  # Bigger and bold y-axis labels
      axis.title.x = element_text(size = 16, face = "bold"),  # Bigger and bold x-axis title
      axis.title.y = element_text(size = 16, face = "bold"),  # Bigger and bold y-axis title
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Bigger, centered, bold title
      legend.title = element_text(size = 16, face = "bold"),  # Bigger and bold legend title
      legend.text = element_text(size = 14),  # Bigger legend text
      strip.text = element_text(size = 16, face = "bold")  # Bigger facet labels
    )
  
  
  
  print(plot)
  return(combined_data)
}

# Extract significant metabolites from Wilcoxon test (FDR-adjusted p-value < 0.05)
significant_metabolites_wilcox <- visit_original_res %>%
  filter(Wilcoxon_Adjusted_P < 0.05) %>%
  pull(Metabolite)  # Extract the names of significant metabolites



#call function to plot distirbution
#Halfmin
# Visit 1
dist1_Halfmin_10pct <- plot_distribution(visit1_data, Halfmin_10pct_1, "Half-min", 10, 1, significant_metabolites_wilcox)
dist1_Halfmin_20pct <- plot_distribution(visit1_data, Halfmin_20pct_1, "Half-min", 20, 1, significant_metabolites_wilcox)
dist1_Halfmin_30pct <- plot_distribution(visit1_data, Halfmin_30pct_1, "Half-min", 30, 1, significant_metabolites_wilcox)
dist1_Halfmin_40pct <- plot_distribution(visit1_data, Halfmin_40pct_1, "Half-min", 40, 1, significant_metabolites_wilcox)

# Visit 2
dist2_Halfmin_10pct <- plot_distribution(visit2_data, Halfmin_10pct_2, "Half-min", 10, 2, significant_metabolites_wilcox)
dist2_Halfmin_20pct <- plot_distribution(visit2_data, Halfmin_20pct_2, "Half-min", 20, 2, significant_metabolites_wilcox)
dist2_Halfmin_30pct <- plot_distribution(visit2_data, Halfmin_30pct_2, "Half-min", 30, 2, significant_metabolites_wilcox)
dist2_Halfmin_40pct <- plot_distribution(visit2_data, Halfmin_40pct_2, "Half-min", 40, 2, significant_metabolites_wilcox)

#KNN
# Visit 1
dist1_KNN_10pct <- plot_distribution(visit1_data, KNN_10pct_1, "KNN", 10, 1, significant_metabolites_wilcox)
dist1_KNN_20pct <- plot_distribution(visit1_data, KNN_20pct_1, "KNN", 20, 1, significant_metabolites_wilcox)
dist1_KNN_30pct <- plot_distribution(visit1_data, KNN_30pct_1, "KNN", 30, 1, significant_metabolites_wilcox)
dist1_KNN_40pct <- plot_distribution(visit1_data, KNN_40pct_1, "KNN", 40, 1, significant_metabolites_wilcox)

# Visit 2
dist2_KNN_10pct <- plot_distribution(visit2_data, KNN_10pct_2, "KNN", 10, 2, significant_metabolites_wilcox)
dist2_KNN_20pct <- plot_distribution(visit2_data, KNN_20pct_2, "KNN", 20, 2, significant_metabolites_wilcox)
dist2_KNN_30pct <- plot_distribution(visit2_data, KNN_30pct_2, "KNN", 30, 2, significant_metabolites_wilcox)
dist2_KNN_40pct <- plot_distribution(visit2_data, KNN_40pct_2, "KNN", 40, 2, significant_metabolites_wilcox)


#RF
# Visit 1
dist1_RF_10pct <- plot_distribution(visit1_data, RF_10pct_1, "RF", 10, 1, significant_metabolites_wilcox)
dist1_RF_20pct <- plot_distribution(visit1_data, RF_20pct_1, "RF", 20, 1, significant_metabolites_wilcox)
dist1_RF_30pct <- plot_distribution(visit1_data, RF_30pct_1, "RF", 30, 1, significant_metabolites_wilcox)
dist1_RF_40pct <- plot_distribution(visit1_data, RF_40pct_1, "RF", 40, 1, significant_metabolites_wilcox)

# Visit 2
dist2_RF_10pct <- plot_distribution(visit2_data, RF_10pct_2, "RF", 10, 2, significant_metabolites_wilcox)
dist2_RF_20pct <- plot_distribution(visit2_data, RF_20pct_2, "RF", 20, 2, significant_metabolites_wilcox)
dist2_RF_30pct <- plot_distribution(visit2_data, RF_30pct_2, "RF", 30, 2, significant_metabolites_wilcox)
dist2_RF_40pct <- plot_distribution(visit2_data, RF_40pct_2, "RF", 40, 2, significant_metabolites_wilcox)


#QRILC
#Visit 1
dist1_QRILC_10pct <- plot_distribution(visit1_data, QRILC_10pct_1, "QRILC", 10,1,significant_metabolites_wilcox)
dist1_QRILC_20pct <- plot_distribution(visit1_data, QRILC_20pct_1, "QRILC", 20,1,significant_metabolites_wilcox)
dist1_QRILC_30pct <- plot_distribution(visit1_data, QRILC_30pct_1, "QRILC", 30,1,significant_metabolites_wilcox)
dist1_QRILC_40pct <- plot_distribution(visit1_data, QRILC_40pct_1, "QRILC", 40,1,significant_metabolites_wilcox)
#Visit 2
dist2_QRILC_10pct <- plot_distribution(visit2_data, QRILC_10pct_2, "QRILC", 10,2,significant_metabolites_wilcox)
dist2_QRILC_20pct <- plot_distribution(visit2_data, QRILC_20pct_2, "QRILC", 20,2, significant_metabolites_wilcox)
dist2_QRILC_30pct <- plot_distribution(visit2_data, QRILC_30pct_2, "QRILC", 30,2,significant_metabolites_wilcox)
dist2_QRILC_40pct <- plot_distribution(visit2_data, QRILC_40pct_2, "QRILC", 40,2, significant_metabolites_wilcox)


# ------------------------------------
# Part 5.1: Distribution Plot (Whole Dataset)
# ------------------------------------

#function to plot distribution before and after imputation (entire dataset)
plot_whole_distribution <- function(original, imputed, method, percentage, visit) {
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
    labs(title = paste0("Overall Density Distribution Before and After ", method, " Imputation (", percentage, "% Missing, Visit ", ")"),
         x = "Value",
         y = "Density") +
    geom_vline(data = mean_data %>% filter(Data == "Original Data"),
               aes(xintercept = mean_value, color = "Original Mean"), linewidth = 0.5, linetype = "dashed") + 
    geom_vline(data = mean_data %>% filter(Data == "Imputed Data"),
               aes(xintercept = mean_value, color = "Imputed Mean"), linewidth = 0.5, linetype = "dashed") +
    scale_color_manual(name = "Mean Values", 
                       values = c("Original Mean" = "blue", "Imputed Mean" = "red")) + 
    xlim(-100, 400) + 
    theme(legend.position = "bottom")
  
  
  print(plot)
  return(combined_data)
}

#call function to plot whole distirbution
#Halfmin
#visit 1
whole_dist1_Halfmin_10pct <- plot_whole_distribution(visit1_data, Halfmin_10pct_1, "Half-min", 10, 1)
whole_dist1_Halfmin_40pct <- plot_whole_distribution(visit1_data, Halfmin_40pct_1, "Half-min", 40, 1)
#visit 2
whole_dist2_Halfmin_10pct <- plot_whole_distribution(visit2_data, Halfmin_10pct_2, "Half-min", 10, 2)
whole_dist2_Halfmin_40pct <- plot_whole_distribution(visit2_data, Halfmin_40pct_2, "Half-min", 40, 2)

#KNN
#visit 1
whole_dist1_KNN_10pct <- plot_whole_distribution(visit1_data, KNN_10pct_1, "KNN", 10, 1)
whole_dist1_KNN_40pct <- plot_whole_distribution(visit1_data, KNN_40pct_1, "KNN", 40, 1)
#visit 2
whole_dist2_KNN_10pct <- plot_whole_distribution(visit2_data, KNN_10pct_2, "KNN", 10, 2)
whole_dist2_KNN_40pct <- plot_whole_distribution(visit2_data, KNN_40pct_2, "KNN", 40, 2)

#RF
#visit 1
whole_dist1_RF_10pct <- plot_whole_distribution(visit1_data, RF_10pct_1, "RF", 10, 1)
whole_dist1_RF_40pct <- plot_whole_distribution(visit1_data, RF_40pct_1, "RF", 40, 1)
#visit 2
whole_dist2_RF_10pct <- plot_whole_distribution(visit2_data, RF_10pct_2, "RF", 10, 2)
whole_dist2_RF_40pct <- plot_whole_distribution(visit2_data, RF_40pct_2, "RF", 40, 2)

#QRILC
#visit 1
whole_dist1_QRILC_10pct <- plot_whole_distribution(visit1_data, QRILC_10pct_1, "QRILC", 10, 1)
whole_dist1_QRILC_40pct <- plot_whole_distribution(visit1_data, QRILC_40pct_1, "QRILC", 40, 1)
#visit 2
whole_dist2_QRILC_10pct <- plot_whole_distribution(visit2_data, QRILC_10pct_2, "QRILC", 10, 2)
whole_dist2_QRILC_40pct <- plot_whole_distribution(visit2_data, QRILC_40pct_2, "QRILC", 40, 2)

# ------------------------------------
# Part 6: ANOVA
# ------------------------------------

#fit an ANOVA model with an interaction term between imputation method and missingness level 
#apply log to normalize data


#ANOVA for Visit 1 and VIsit 2 (effect of imputation and mcar proportion)
#Visit 1
anova_visit1 <- aov(log(Weighted_NRMSE) ~ Imputation_Method * MCAR_Proportion, data = nrmse_data1)
summary(anova_visit1)  

#check assumptions (Homogeneity of Variance)
#leveneTest(Weighted_NRMSE ~ Imputation_Method * MCAR_Proportion, data = nrmse_data1)

#ANOVA for Visit 2
anova_visit2 <- aov(log(Weighted_NRMSE) ~ Imputation_Method * MCAR_Proportion, data = nrmse_data2)
summary(anova_visit2)  



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

#visit 1
ggplot(data.frame(residuals = residuals(anova_visit1)), aes(x = residuals)) +
  geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency") +
  theme_minimal()
#visit 2
ggplot(data.frame(residuals = residuals(anova_visit2)), aes(x = residuals)) +
  geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency") +
  theme_minimal()

#Q-Q plot of residuals
#plot residula against theoretical normal distirbutions
#red line shows perfect normal distirbution
#S-shaped pattern = possibles skewness

#visit 1
qqnorm(residuals(anova_visit1), col = "blue")
qqline(residuals(anova_visit1), col = "red")

#visit 2
qqnorm(residuals(anova_visit2), col = "blue")
qqline(residuals(anova_visit2), col = "red")

#Tukey-Anscombe plot to check residuals vs fitted values
#x-axis = fitted values (predicted) and y-axis = residueals (error)
#residuals should be evenly spread, if the points fan out or forma pattern the assumption of homoscedascity is violated 

#visit 1
plot(fitted(anova_visit1), resid(anova_visit1), 
     main = "Tukey-Anscombe Plot", 
     col = "blue", 
     xlab = "Fitted Values (Predicted by ANOVA Model)", 
     ylab = "Residuals (Errors)")

#run shapiro test on anova model
shapiro.test(resid(anova_visit1))

#visit 2
plot(fitted(anova_visit2), resid(anova_visit2), 
     main = "Tukey-Anscombe Plot", 
     col = "blue", 
     xlab = "Fitted Values (Predicted by ANOVA Model)", 
     ylab = "Residuals (Errors)")

#run shapiro test on anova model
shapiro.test(resid(anova_visit2))

# ------------------------------------
# Part 7: Kruskal-Wallis
# ------------------------------------

#perform kruskal-wallis test

#visit 1
v1_kruskal_test <- kruskal.test(Weighted_NRMSE ~ Imputation_Method, data = nrmse_data1)
v1_kruskal_test2 <- kruskal.test(Weighted_NRMSE ~ MCAR_Proportion, data = nrmse_data1)
#show results
print(v1_kruskal_test)
print(v1_kruskal_test2) 

#visit 2
v2_kruskal_test <- kruskal.test(Weighted_NRMSE ~ Imputation_Method, data = nrmse_data2)
v2_kruskal_test2 <- kruskal.test(Weighted_NRMSE ~ MCAR_Proportion, data = nrmse_data2)
#show results
print(v2_kruskal_test)
print(v2_kruskal_test2) 


#p-value storngly signficant 
#chi-square = 98.713 (higher value means larger difference between groups)
#at least one imputation method significantly differs from the others in terms of NRMSE

# ------------------------------------
# Part 8: Dunn's Test for each imputation
# ------------------------------------

#perform Dunn's Test for pairwise comparison (BH correction for multiple testing)
#visit 1
v1_dunn_test <- dunnTest(Weighted_NRMSE ~ Imputation_Method, data = nrmse_data1, method = "bh")
v1_dunn_test2 <- dunnTest(Weighted_NRMSE ~ MCAR_Proportion, data = nrmse_data1, method = "bh") #for percentage of mcar

#print the results
print(v1_dunn_test)
print(v2_dunn_test)

#visit 2
v2_dunn_test <- dunnTest(Weighted_NRMSE ~ Imputation_Method, data = nrmse_data2, method = "bh")
v2_dunn_test2 <- dunnTest(Weighted_NRMSE ~ MCAR_Proportion, data = nrmse_data2, method = "bh") #for percentage of mcar

#print the results
print(v2_dunn_test)
print(v2_dunn_test)

#plot Dunns test results
ggplot(nrmse_data1, aes(x = Imputation_Method, y = Weighted_NRMSE, fill = Imputation_Method)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Pairwise Comparisons of Imputation Methods",
       x = "Imputation Method",
       y = "Weighted NRMSE") +
  ylim(0,1)

#plot Dunns test results
ggplot(nrmse_data2, aes(x = Imputation_Method, y = Weighted_NRMSE, fill = Imputation_Method)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Pairwise Comparisons of Imputation Methods",
       x = "Imputation Method",
       y = "Weighted NRMSE") +
  ylim(0,1)

# ------------------------------------
# Part 9: Compare Metabolite 
# ------------------------------------

#function to compute summary statistics
compute_stats <- function(data, visit_label, imputation_method, mcar_percent) {
  data_long <- data %>%
    pivot_longer(cols = 6:ncol(data), names_to = "Metabolite", values_to = "Value") %>%
    mutate(Visit = visit_label, Imputation_Method = imputation_method, MCAR_Percent = mcar_percent)
  
  stats <- data_long %>%
    group_by(Metabolite, Visit, Imputation_Method, MCAR_Percent) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      .groups = 'drop'
    )
  
  return(stats)
}

#compute original statistics for Visit 1 and Visit 2
original_visit1_stats <- compute_stats(visit1_data, "Visit 1", "Original", "0%")
original_visit2_stats <- compute_stats(visit2_data, "Visit 2", "Original", "0%")

#halfmin
#visit 1
halfmin10_visit1_stats <- compute_stats(Halfmin_10pct_1, "Visit 1", "Imputed", "10%")
halfmin20_visit1_stats <- compute_stats(Halfmin_20pct_1, "Visit 1", "Imputed", "20%")
halfmin30_visit1_stats <- compute_stats(Halfmin_30pct_1, "Visit 1", "Imputed", "30%")
halfmin40_visit1_stats <- compute_stats(Halfmin_40pct_1, "Visit 1", "Imputed", "40%")
#visit 2
halfmin10_visit2_stats <- compute_stats(Halfmin_10pct_2, "Visit 2", "Imputed", "10%")
halfmin20_visit2_stats <- compute_stats(Halfmin_20pct_2, "Visit 2", "Imputed", "20%")
halfmin30_visit2_stats <- compute_stats(Halfmin_30pct_2, "Visit 2", "Imputed", "30%")
halfmin40_visit2_stats <- compute_stats(Halfmin_40pct_2, "Visit 2", "Imputed", "40%")

#knn
#visit 1
knn10_visit1_stats <- compute_stats(KNN_10pct_1, "Visit 1", "Imputed", "10%")
knn20_visit1_stats <- compute_stats(KNN_20pct_1, "Visit 1", "Imputed", "20%")
knn30_visit1_stats <- compute_stats(KNN_30pct_1, "Visit 1", "Imputed", "30%")
knn40_visit1_stats <- compute_stats(KNN_40pct_1, "Visit 1", "Imputed", "40%")
#visit 2
knn10_visit2_stats <- compute_stats(KNN_10pct_2, "Visit 2", "Imputed", "10%")
knn20_visit2_stats <- compute_stats(KNN_20pct_2, "Visit 2", "Imputed", "20%")
knn30_visit2_stats <- compute_stats(KNN_30pct_2, "Visit 2", "Imputed", "30%")
knn40_visit2_stats <- compute_stats(KNN_40pct_2, "Visit 2", "Imputed", "40%")

#RF
#visit 1
rf10_visit1_stats <- compute_stats(RF_10pct_1, "Visit 1", "Imputed", "10%")
rf20_visit1_stats <- compute_stats(RF_20pct_1, "Visit 1", "Imputed", "20%")
rf30_visit1_stats <- compute_stats(RF_30pct_1, "Visit 1", "Imputed", "30%")
rf40_visit1_stats <- compute_stats(RF_40pct_1, "Visit 1", "Imputed", "40%")
#visit 2
rf10_visit2_stats <- compute_stats(RF_10pct_2, "Visit 2", "Imputed", "10%")
rf20_visit2_stats <- compute_stats(RF_20pct_2, "Visit 2", "Imputed", "20%")
rf30_visit2_stats <- compute_stats(RF_30pct_2, "Visit 2", "Imputed", "30%")
rf40_visit2_stats <- compute_stats(RF_40pct_2, "Visit 2", "Imputed", "40%")

#QRILC
#visit 1
qrilc10_visit1_stats <- compute_stats(QRILC_10pct_1, "Visit 1", "Imputed", "10%")
qrilc20_visit1_stats <- compute_stats(QRILC_20pct_1, "Visit 1", "Imputed", "20%")
qrilc30_visit1_stats <- compute_stats(QRILC_30pct_1, "Visit 1", "Imputed", "30%")
qrilc40_visit1_stats <- compute_stats(QRILC_40pct_1, "Visit 1", "Imputed", "40%")
#visit 2
qrilc10_visit2_stats <- compute_stats(QRILC_10pct_2, "Visit 2", "Imputed", "10%")
qrilc20_visit2_stats <- compute_stats(QRILC_20pct_2, "Visit 2", "Imputed", "20%")
qrilc30_visit2_stats <- compute_stats(QRILC_30pct_2, "Visit 2", "Imputed", "30%")
qrilc40_visit2_stats <- compute_stats(QRILC_40pct_2, "Visit 2", "Imputed", "40%")

# Combine statistics for each imputation method separately, including the original (0%) data
halfmin_stats <- bind_rows(
  original_visit1_stats, original_visit2_stats,
  halfmin10_visit1_stats, halfmin20_visit1_stats, halfmin30_visit1_stats, halfmin40_visit1_stats,
  halfmin10_visit2_stats, halfmin20_visit2_stats, halfmin30_visit2_stats, halfmin40_visit2_stats
)

knn_stats <- bind_rows(
  original_visit1_stats, original_visit2_stats,
  knn10_visit1_stats, knn20_visit1_stats, knn30_visit1_stats, knn40_visit1_stats,
  knn10_visit2_stats, knn20_visit2_stats, knn30_visit2_stats, knn40_visit2_stats
)

rf_stats <- bind_rows(
  original_visit1_stats, original_visit2_stats,
  rf10_visit1_stats, rf20_visit1_stats, rf30_visit1_stats, rf40_visit1_stats,
  rf10_visit2_stats, rf20_visit2_stats, rf30_visit2_stats, rf40_visit2_stats
)

qrilc_stats <- bind_rows(
  original_visit1_stats, original_visit2_stats,
  qrilc10_visit1_stats, qrilc20_visit1_stats, qrilc30_visit1_stats, qrilc40_visit1_stats,
  qrilc10_visit2_stats, qrilc20_visit2_stats, qrilc30_visit2_stats, qrilc40_visit2_stats
)

# Convert MCAR_Percent to numeric (keeping 0% for original data)
halfmin_stats$MCAR_Percent <- as.numeric(gsub("%", "", halfmin_stats$MCAR_Percent))
knn_stats$MCAR_Percent <- as.numeric(gsub("%", "", knn_stats$MCAR_Percent))
rf_stats$MCAR_Percent <- as.numeric(gsub("%", "", rf_stats$MCAR_Percent))
qrilc_stats$MCAR_Percent <- as.numeric(gsub("%", "", qrilc_stats$MCAR_Percent))

# Define function to generate plots for Mean, Median, and SD for a given imputation method
plot_stats <- function(data, method) {
  plot_mean <- ggplot(data, aes(x = MCAR_Percent, y = Mean, color = Visit, group = interaction(Metabolite, Visit))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_wrap(~Metabolite, scales = "free_y") +
    theme_minimal() +
    labs(title = paste("Mean Changes -", method, "Imputation (Including Original)"),
         x = "MCAR Percentage",
         y = "Mean") +
    scale_color_manual(values = c("Visit 1" = "blue", "Visit 2" = "red")) +
    theme(legend.position = "bottom")
  
  plot_median <- ggplot(data, aes(x = MCAR_Percent, y = Median, color = Visit, group = interaction(Metabolite, Visit))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_wrap(~Metabolite, scales = "free_y") +
    theme_minimal() +
    labs(title = paste("Median Changes -", method, "Imputation (Including Original)"),
         x = "MCAR Percentage",
         y = "Median") +
    scale_color_manual(values = c("Visit 1" = "blue", "Visit 2" = "red")) +
    theme(legend.position = "bottom")
  
  plot_sd <- ggplot(data, aes(x = MCAR_Percent, y = SD, color = Visit, group = interaction(Metabolite, Visit))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_wrap(~Metabolite, scales = "free_y") +
    theme_minimal() +
    labs(title = paste("SD Changes -", method, "Imputation (Including Original)"),
         x = "MCAR Percentage",
         y = "Standard Deviation") +
    scale_color_manual(values = c("Visit 1" = "blue", "Visit 2" = "red")) +
    theme(legend.position = "bottom")
  
  print(plot_mean)
  print(plot_median)
  print(plot_sd)
}

# Generate plots for each imputation method separately, now including original (0%) values
plot_stats(halfmin_stats, "Halfmin")
plot_stats(knn_stats, "KNN")
plot_stats(rf_stats, "RF")
plot_stats(qrilc_stats, "QRILC")


#MEAN DIFFERENCE PLOT 

#extract significant metabolites based on Wilcoxon test (FDR-adjusted p-value < 0.05)
significant_metabolites_wilcox <- visit_original_res %>%
  filter(Wilcoxon_Adjusted_P < 0.05) %>%
  pull(Metabolite)  # Extract the names of significant metabolites


# Function to compute mean differences
compute_mean_difference <- function(data, method) {
  mean_values <- data %>%
    pivot_longer(cols = 6:ncol(.), names_to = "Metabolite", values_to = "Value") %>%
    group_by(Metabolite, Visit) %>%
    summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Visit, values_from = Mean_Value) %>%
    mutate(Method = method,
           Mean_Difference = `Visit 2` - `Visit 1`)
  
  return(mean_values)
}

# Compute mean differences for all methods
mean_differences <- bind_rows(
  compute_mean_difference(FAO_original, "Original"),
  compute_mean_difference(Halfmin_10pct_tot, "Halfmin (10%)"),
  compute_mean_difference(Halfmin_20pct_tot, "Halfmin (20%)"),
  compute_mean_difference(Halfmin_30pct_tot, "Halfmin (30%)"),
  compute_mean_difference(Halfmin_40pct_tot, "Halfmin (40%)"),
  compute_mean_difference(KNN_10pct_tot, "KNN (10%)"),
  compute_mean_difference(KNN_20pct_tot, "KNN (20%)"),
  compute_mean_difference(KNN_30pct_tot, "KNN (30%)"),
  compute_mean_difference(KNN_40pct_tot, "KNN (40%)"),
  compute_mean_difference(RF_10pct_tot, "RF (10%)"),
  compute_mean_difference(RF_20pct_tot, "RF (20%)"),
  compute_mean_difference(RF_30pct_tot, "RF (30%)"),
  compute_mean_difference(RF_40pct_tot, "RF (40%)"),
  compute_mean_difference(QRILC_10pct_tot, "QRILC (10%)"),
  compute_mean_difference(QRILC_20pct_tot, "QRILC (20%)"),
  compute_mean_difference(QRILC_30pct_tot, "QRILC (30%)"),
  compute_mean_difference(QRILC_40pct_tot, "QRILC (40%)")
)

# Convert Method to factor to maintain ordering
mean_differences$Method <- factor(mean_differences$Method, 
                                  levels = c("Original", 
                                             "Halfmin (10%)", "Halfmin (20%)", "Halfmin (30%)", "Halfmin (40%)",
                                             "KNN (10%)", "KNN (20%)", "KNN (30%)", "KNN (40%)",
                                             "RF (10%)", "RF (20%)", "RF (30%)", "RF (40%)",
                                             "QRILC (10%)", "QRILC (20%)", "QRILC (30%)", "QRILC (40%)"))

# Filter mean differences to include only significant metabolites from Wilcoxon test
filtered_mean_differences <- mean_differences %>%
  filter(Metabolite %in% significant_metabolites_wilcox)

# Function to plot mean differences for a specific imputation method (only significant metabolites)
plot_mean_differences_per_method <- function(data, method_name) {
  ggplot(data %>% filter(Imputation_Method %in% c("Original", method_name)), 
         aes(x = Metabolite, y = Mean_Difference, fill = Percentage)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste("Mean Differences Between Visit 1 and Visit 2 (", method_name, ")", sep=""),
         x = "Metabolite",
         y = "Mean Difference (Visit 2 - Visit 1)",
         fill = "Missingness %") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    )
}

# Generate separate plots for each imputation method using the filtered data
plot_halfmin <- plot_mean_differences_per_method(filtered_mean_differences, "Halfmin")
plot_knn <- plot_mean_differences_per_method(filtered_mean_differences, "KNN")
plot_rf <- plot_mean_differences_per_method(filtered_mean_differences, "RF")
plot_qrilc <- plot_mean_differences_per_method(filtered_mean_differences, "QRILC")

# Display the filtered plots
print(plot_halfmin)
print(plot_knn)
print(plot_rf)
print(plot_qrilc)


# Function to modify plots for layout adjustments
adjust_plot_for_layout <- function(plot, remove_x_axis = TRUE, remove_y_axis = TRUE, remove_legend = TRUE) {
  plot + theme(
    axis.text.x = if (remove_x_axis) element_blank() else element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, face = "bold"),
    axis.title.x = if (remove_x_axis) element_blank() else element_text(size = 14, face = "bold"),
    axis.title.y = if (remove_y_axis) element_blank() else element_text(size = 14, face = "bold"),
    legend.position = if (remove_legend) "none" else "right"  # Keep legend only for the last plot
  )
}

# Modify plots: Remove x-axis from all except the bottom row
plot_halfmin_adj <- adjust_plot_for_layout(plot_halfmin, remove_x_axis = TRUE, remove_y_axis = TRUE, remove_legend = TRUE)  # Keep y-axis here
plot_knn_adj <- adjust_plot_for_layout(plot_knn, remove_x_axis = TRUE, remove_y_axis = TRUE, remove_legend = TRUE)
plot_rf_adj <- adjust_plot_for_layout(plot_rf, remove_x_axis = TRUE, remove_y_axis = FALSE, remove_legend = TRUE)
plot_qrilc_adj <- adjust_plot_for_layout(plot_qrilc, remove_x_axis = FALSE, remove_y_axis = TRUE, remove_legend = FALSE)  # Keep legend here

# Combine plots in a vertical layout with shared axes
combined_plot <- plot_halfmin_adj / plot_knn_adj / plot_rf_adj / plot_qrilc_adj + 
  plot_layout(ncol = 1, guides = "collect") &  # Ensures only one legend
  theme(legend.position = "right")  # Moves legend to the side

# Display the final combined plot
print(combined_plot)



