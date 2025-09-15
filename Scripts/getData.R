# script to get data from airway package

setwd("C:/Users/dasan/Downloads/Bioinformatics_Project/Differencial_expression_analysis/Scripts")

rm(list = ls())

# Installing the library
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomeInfoDbData")
BiocManager::install("airway")

# Loading the library
library(airway)

# Loading the data
data(airway)
airway

# Cleaning the saving sample data
sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info,
            file = "C:/Users/dasan/Downloads/Bioinformatics_Project/Differencial_expression_analysis/Data/sample_info.csv",
            sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE)

# Cleaning the saving the counts data
countsData <- assay(airway)
write.table(countsData, file = "C:/Users/dasan/Downloads/Bioinformatics_Project/Differencial_expression_analysis/Data/counts_data.csv", sep = ',', col.names = TRUE, row.names = TRUE, quote = FALSE)
