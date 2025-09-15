# script to perform differential gene expression analysis using DESeq2 package

setwd("C:/Users/dasan/Downloads/Bioinformatics_Project/Differencial_expression_analysis/Scripts")

rm(list = ls())

# Installing the package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# load libraries
library(DESeq2)
library(tidyverse)
library(airway)
library(ggplot2)

# Step 1: preparing count data ----------------

# read in counts data
counts_data <- read.csv("C:/Users/dasan/Downloads/Bioinformatics_Project/Differencial_expression_analysis/Data/counts_data.csv")
head(counts_data)


# read in sample info
colData <- read.csv("C:/Users/dasan/Downloads/Bioinformatics_Project/Differencial_expression_analysis/Data/sample_info.csv")


# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# checking if the columns are in the same order
all(colnames(counts_data) == rownames(colData))


# Step 2: construction of DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res



# Explore Results ----------------

summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)

# e.g.: treated_4hrs, treated_8hrs, untreated

results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

# output directory for plots
out_dir <- "C:/Users/dasan/Downloads/Bioinformatics_Project/Differencial_expression_analysis/Plots/"
dir.create(out_dir, showWarnings = FALSE)  # create folder if not exists

# MA plot
png(filename = paste0(out_dir, "MA_plot.png"), width = 2000, height = 1600, res = 300)
plotMA(res)
dev.off()


