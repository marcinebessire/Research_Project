#Load required library
library(tidyverse)
library(ggplot2)
library(corrplot)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)

#load imputed data
halfmin23 <- read.csv("/Users/marcinebessire/Desktop/project/HalfMin_Imputation23.csv", check.names = FALSE)
halfmin24 <- read.csv("/Users/marcinebessire/Desktop/project/HalfMin_Imputation24.csv", check.names = FALSE)

knn23 <- read.csv("/Users/marcinebessire/Desktop/project/KNN_Imputation23.csv", check.names = FALSE)
knn24 <- read.csv("/Users/marcinebessire/Desktop/project/KNN_Imputation24.csv", check.names = FALSE)

QRILC23 <- read.csv("/Users/marcinebessire/Desktop/project/QRILC_Imputation23.csv", check.names = FALSE)
QRILC24 <- read.csv("/Users/marcinebessire/Desktop/project/QRILC_Imputation24.csv", check.names = FALSE)

RF23 <- read.csv("/Users/marcinebessire/Desktop/project/RF_Imputation23.csv", check.names = FALSE)
RF24 <- read.csv("/Users/marcinebessire/Desktop/project/RF_Imputation24.csv", check.names = FALSE)

# Part 1 -----
# Kolomogorov-Smirnov test: nonparametric test to test whether two samples came from same distribution

#HALFMIN
#vs KNN
#2023
ks.test(halfmin23$Imputed_Data, knn23$Imputed_Data) #p-value = 0.6209
#2024
ks.test(halfmin24$Imputed_Data, knn24$Imputed_Data) #p-value = 7.398e-12
#vs QRILC
#2023
ks.test(halfmin23$Imputed_Data, QRILC23$Imputed_Data) #p-value = 0.9994
#2024
ks.test(halfmin24$Imputed_Data, QRILC24$Imputed_Data) #p-value = 0.03766
#vs RF
#2023
ks.test(halfmin23$Imputed_Data, RF23$Imputed_Data) #p-value = 0.762
#2024
ks.test(halfmin24$Imputed_Data, RF24$Imputed_Data) #p-value = 0.0001061

#KNN 
#vs RF 
#2023
ks.test(knn23$Imputed_Data, RF23$Imputed_Data) #p-value = 1
#2024
ks.test(knn24$Imputed_Data, RF24$Imputed_Data) #p-value = 6.902e-10

#vs QRILC
#2023
ks.test(knn23$Imputed_Data, QRILC23$Imputed_Data) #p-value = 0.9585
#2024
ks.test(knn24$Imputed_Data, QRILC24$Imputed_Data) #p-value = 1.425e-11

#QRILC 
#vs RF
#2023
ks.test(QRILC23$Imputed_Data, RF23$Imputed_Data) #p-value = 0.985
#2024
ks.test(QRILC24$Imputed_Data, RF24$Imputed_Data) #p-value = 0.2602



