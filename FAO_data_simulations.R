#load necessary libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(impute) #for knn
library(reshape2) #for melt

#load FAO data
FAO_data <- read.csv("/Users/marcinebessire/Desktop/project/FAO_data.csv", check.names = FALSE)

# ---------------------------------
# Part 0: Plot Whole Distribution
# ---------------------------------

# Function to plot the overall distribution of all numeric values
plot_overall_distribution <- function(original) {
  # Select numeric columns
  numeric_original <- original[, 6:ncol(original)]
  
  # Convert the entire numeric data into a single vector
  all_values <- unlist(numeric_original, use.names = FALSE)
  
  # Create a data frame for plotting
  df_all_values <- data.frame(Value = all_values)
  
  # Plot the overall density distribution
  plot <- ggplot(df_all_values, aes(x = Value)) +
    geom_density(fill = "gray", alpha = 0.5, color = "black") +
    theme_minimal() +
    labs(title = "Overall Distribution of FAO Data",
         x = "Measurement",
         y = "Density") +
    xlim(-100,500)
  
  print(plot)
}

# Example usage
plot_overall_distribution(FAO_data)


# ---------------------------------
# Part 1: Plot Distribution of Metabolites 
# ---------------------------------

#plot original distribution as density plots
plot_original_distribution <- function(data) {
  #select numeric columns
  numeric_data <- data[, 6:ncol(data)]
  
  #convert data to long format for ggplot
  long_data <- pivot_longer(numeric_data, cols = everything(), names_to = "Variable", values_to = "Metabolites")
  
  #compute mean and median for each variable 
  stats <- long_data %>%
    group_by(Variable) %>%
    summarize(
      Mean = mean(Metabolites, na.rm = TRUE), #mean
      Median = median(Metabolites, na.rm = TRUE) #median
    )
  
  #plot density plots
  ggplot(long_data, aes(x = Metabolites)) +
    geom_density(fill = "lightblue", alpha = 0.7, color = "blue") +  
    geom_vline(data = stats, aes(xintercept = Mean, color = "Mean"), linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = stats, aes(xintercept = Median, color = "Median"), linetype = "solid", linewidth = 0.5) +
    facet_wrap(~Variable, scales = "free") +  #plot for each variable
    theme_minimal() +
    labs(title = "Original Data Distribution (Density Plots)", x = "Metabolites", y = "Density") +
    scale_color_manual(values = c("Mean"= "red", "Median" = "darkgreen"), name = "Statistics") +
    theme(
      strip.text = element_text(size = 10, face = "bold"), 
      legend.position = "bottom"
    )
}


#plot the original distribution
plot_original_distribution(FAO_data)


# ---------------------------------
# Part 2: Generate MCAR in data set
# MCAR because randomly assign MV to each column without any dependency on the observed or unobserved data
# ---------------------------------

#function to introduce MCAR missing values with balance across Visit 1 and Visit 2
MCAR_manipulation_balanced <- function(data, missing_percentage){
  #copy dataset to avoid modifying the original
  data_copy <- data
  
  #get the Visit 1 and Visit 2 indices
  visit1_indices <- which(data_copy$Visit == "Visit 1")
  visit2_indices <- which(data_copy$Visit == "Visit 2")
  
  #iterate through each numeric column
  for (col in 6:ncol(data_copy)) {
    #get number of missing values to introduce (rounded down to ensure pairing)
    num_mv_total <- round(nrow(data_copy) * missing_percentage)
    num_mv <- floor(num_mv_total / 2)  #even split between Visit 1 and Visit 2
    
    if (num_mv > 0 && length(visit1_indices) > 0 && length(visit2_indices) > 0) {
      #randomly choose indices for missing values in each visit group
      missing_indices_v1 <- sample(visit1_indices, num_mv, replace = FALSE)
      missing_indices_v2 <- sample(visit2_indices, num_mv, replace = FALSE)
      
      #set selected values to NA
      data_copy[missing_indices_v1, col] <- NA
      data_copy[missing_indices_v2, col] <- NA
    }
  }
  
  return(data_copy)
}

#call function to generate MCAR (evenly between visit 1 and visit 2)
FAO_data_10pct <- MCAR_manipulation_balanced(FAO_data, 0.1) #10% missing values
FAO_data_20pct <- MCAR_manipulation_balanced(FAO_data, 0.2) #20% missing values
FAO_data_30pct <- MCAR_manipulation_balanced(FAO_data, 0.3) #30% missing values
FAO_data_40pct <- MCAR_manipulation_balanced(FAO_data, 0.4) #40% missing values


#save them to csv file 
write_csv(FAO_data_10pct, "/Users/marcinebessire/Desktop/project/FAO_10pct.csv")
write_csv(FAO_data_20pct, "/Users/marcinebessire/Desktop/project/FAO_20pct.csv")
write_csv(FAO_data_30pct, "/Users/marcinebessire/Desktop/project/FAO_30pct.csv")
write_csv(FAO_data_40pct, "/Users/marcinebessire/Desktop/project/FAO_40pct.csv")


#plot missing values
plot_missing_values <- function(data, title) {
  #convert data to long format
  data_long <- melt(data, id.vars = c("ID", "Participant", "MonthDay", "Year", "Visit"))
  
  #create a column indicating missing values
  data_long$Missing <- ifelse(is.na(data_long$value), "Missing", "Present")
  
  #correctly order Participant IDs 
  sorted_participants <- mixedsort(unique(data_long$Participant))  
  data_long$Participant <- factor(data_long$Participant, levels = sorted_participants)  #apply order
  
  #plot missing values using a heatmap
  ggplot(data_long, aes(x = variable, y = Participant, fill = Missing)) +
    geom_tile(color = "grey") +  #add borders for clarity
    facet_wrap(~Visit, ncol = 2) +  #separate by Visit 1 and Visit 2
    scale_fill_manual(values = c("Present" = "white", "Missing" = "red")) +
    theme_minimal(base_size = 16) +  #set global font size
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),  
      axis.text.y = element_text(size = 12, face = "bold"),  
      axis.title.x = element_text(size = 16, face = "bold"),  
      axis.title.y = element_text(size = 16, face = "bold"),  
      strip.text = element_text(size = 18, face = "bold"),  
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
      legend.title = element_text(size = 16, face = "bold"),  
      legend.text = element_text(size = 14)  
    ) +
    labs(title = title, x = "Metabolite", y = "Patient ID", fill = "Data Status")
}

#call function for plotting
plot_missing_values(FAO_data_10pct, "Missing Data Pattern - 10% MCAR")
plot_missing_values(FAO_data_20pct, "Missing Data Pattern - 20% MCAR")
plot_missing_values(FAO_data_30pct, "Missing Data Pattern - 30% MCAR")
plot_missing_values(FAO_data_40pct, "Missing Data Pattern - 40% MCAR")

# # ---------- do not do it like this --------
# # Function to introduce MCAR missing values with balance across Visit 1 and Visit 2
# MCAR_manipulation_balanced2 <- function(data, missing_percentage){
#   # Copy dataset to avoid modifying the original
#   data_copy <- data
#   
#   # Get Visit 1 and Visit 2 row indices
#   visit1_indices <- which(data_copy$Visit == "Visit 1")
#   visit2_indices <- which(data_copy$Visit == "Visit 2")
#   
#   # Iterate through each numeric column (assuming from 6th column onwards)
#   for (col in 6:ncol(data_copy)) {
#     # Compute total number of missing values to introduce
#     num_mv_total <- max(1, round(nrow(data_copy) * missing_percentage))  # Ensure at least 1 missing value
#     
#     # Split evenly
#     num_mv <- num_mv_total %/% 2  # Integer division to get pairs
#     
#     # Determine if there's an extra missing value (for odd cases)
#     extra_mv <- num_mv_total %% 2  # 1 if odd, 0 if even
#     
#     missing_indices_v1 <- c()
#     missing_indices_v2 <- c()
#     
#     if (num_mv > 0 && length(visit1_indices) > 0 && length(visit2_indices) > 0) {
#       # Assign missing values evenly
#       missing_indices_v1 <- sample(visit1_indices, min(num_mv, length(visit1_indices)), replace = FALSE)
#       missing_indices_v2 <- sample(visit2_indices, min(num_mv, length(visit2_indices)), replace = FALSE)
#     }
#     
#     # If there's an extra missing value, randomly assign it
#     if (extra_mv == 1) {
#       extra_visit <- sample(c("Visit 1", "Visit 2"), 1)  # Choose where to place extra NA
#       if (extra_visit == "Visit 1" && length(visit1_indices) > length(missing_indices_v1)) {
#         extra_index <- sample(setdiff(visit1_indices, missing_indices_v1), 1)
#         missing_indices_v1 <- c(missing_indices_v1, extra_index)
#       } else if (length(visit2_indices) > length(missing_indices_v2)) {
#         extra_index <- sample(setdiff(visit2_indices, missing_indices_v2), 1)
#         missing_indices_v2 <- c(missing_indices_v2, extra_index)
#       }
#     }
#     
#     # Set selected values to NA
#     data_copy[missing_indices_v1, col] <- NA
#     data_copy[missing_indices_v2, col] <- NA
#   }
#   
#   return(data_copy)
# }
# 
# 
# 
# #call function to generate MCAR (evenly between visit 1 and visit 2)
# FAO_data_5pct2 <- MCAR_manipulation_balanced2(FAO_data, 0.05) #5% missing values
# FAO_data_10pct2 <- MCAR_manipulation_balanced2(FAO_data, 0.1) #10% missing values
# FAO_data_20pct2 <- MCAR_manipulation_balanced2(FAO_data, 0.2) #20% missing values
# FAO_data_25pct2 <- MCAR_manipulation_balanced2(FAO_data, 0.25) #25% missing values
# FAO_data_30pct2 <- MCAR_manipulation_balanced2(FAO_data, 0.3) #30% missing values
# FAO_data_40pct2 <- MCAR_manipulation_balanced2(FAO_data, 0.4) #40% missing values

