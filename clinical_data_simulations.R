#load necessary libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(impute) #for knn

#load FAO data
FAO_data <- read.csv("/Users/marcinebessire/Desktop/project/FAO_data.csv", check.names = FALSE)


# ---------------------------------
# Part 1: Plot Distribution of Metabolites 
# ---------------------------------

#plot original distribution as density plots
plot_original_distribution <- function(data) {
  #select numeric columns
  numeric_data <- data[, 6:ncol(data)]
  
  #convert data to long format for ggplot
  long_data <- pivot_longer(numeric_data, cols = everything(), names_to = "Variable", values_to = "Metabolites")
  
  #plot density plots
  ggplot(long_data, aes(x = Metabolites)) +
    geom_density(fill = "lightblue", alpha = 0.7) +  
    facet_wrap(~Variable, scales = "free") +  #plot for each variable
    theme_minimal() +
    labs(title = "Original Data Distribution (Density Plots)", x = "Metabolites", y = "Density") +
    theme(strip.text = element_text(size = 10, face = "bold"))
}


#plot the original distribution
plot_original_distribution(FAO_data)


# ---------------------------------
# Part 2: Generate MCAR in dataset
# MCAR because randomly assign MV to each column without any dependency on the observed or unobserved data
# ---------------------------------
MCAR_manipulation <- function(data, missing_percentage){
  #copy dataset (avoid overwriting)
  data_copy <- data
  
  #iterate through each column
  for (col in 6:ncol(data_copy)) {
    #calculate number of missing value per column
    num_mv <- round(nrow(data_copy) * missing_percentage)
    
    if (num_mv > 0) {
      #select random rows
      missing_indices <- sample(1:nrow(data_copy), num_mv, replace = FALSE)
      data_copy[missing_indices, col] <- NA #set to NA 
    }
  }
  
  return(data_copy)
}

#call function to generate MCAR
FAO_data_5pct <- MCAR_manipulation(FAO_data, 0.05) #5% missing values
FAO_data_10pct <- MCAR_manipulation(FAO_data, 0.1) #10% missing values
FAO_data_20pct <- MCAR_manipulation(FAO_data, 0.2) #20% missing values
FAO_data_25pct <- MCAR_manipulation(FAO_data, 0.25) #25% missing values
FAO_data_30pct <- MCAR_manipulation(FAO_data, 0.3) #30% missing values
FAO_data_40pct <- MCAR_manipulation(FAO_data, 0.4) #40% missing values

#save them to csv file 
write_csv(FAO_data_5pct, "/Users/marcinebessire/Desktop/project/FAO_5pct.csv")
write_csv(FAO_data_10pct, "/Users/marcinebessire/Desktop/project/FAO_10pct.csv")
write_csv(FAO_data_20pct, "/Users/marcinebessire/Desktop/project/FAO_20pct.csv")
write_csv(FAO_data_25pct, "/Users/marcinebessire/Desktop/project/FAO_25pct.csv")
write_csv(FAO_data_30pct, "/Users/marcinebessire/Desktop/project/FAO_30pct.csv")
write_csv(FAO_data_40pct, "/Users/marcinebessire/Desktop/project/FAO_40pct.csv")

