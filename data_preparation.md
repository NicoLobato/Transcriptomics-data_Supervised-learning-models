# Data Preparation

## Loading libraries

```
library("data.table")
library("edgeR")
library("tidyverse")
```

## Reading the data

```
data <- fread("data/allcountsSalmon.csv", sep = ",")
dim(data) # 206574     26
head(data[,0:10])
```

The data contains 206,574 gene expression levels in 24 individuals, named according to their drought resistance as follows:
- Root type: PoT (Tolerant type) / PoS (Sensitive type)
- Stem type: PuT (Tolerant type) / PuS (Sensitive type)
- Treatment type: C (Control type) / S (Drought type)

This data will be prepared and reformatted for further analysis as follows:

```
genes_names <- data$target_id # save the gene IDs
samples_names <- names(data)[3:ncol(data)] # save the sample IDs

data <- as.matrix(data[,-c("V1","target_id")])
rownames(data) <- genes_names
rm(genes_names)
```

Before starting the analysis, the matrix of read counts is filtered so that I keep genes where CPM >= 1 in at least 3 samples (CPM = counts per million).

```
keep=rowSums(cpm(data)>=1)>=3
length(keep)
data <- data[keep,]
rm(keep)
dim(data) # 44881    24
# The number of genes is now significantly reduced
```

We now create the data structure:
- 1st col: Sample IDs
- 2nd col: Treatment type (Control/Drought)
- 3rd col: Stem type (Tolerant/Sensitive)
- 4th col: Root type (Tolerant/Sensitive)
- Remaining cols: Gene expression levels

```
data <- t(data) %>%
  as.data.table(.)

treatment_type <- ifelse(endsWith(samples_names, "C"), "Control", "Drought")
stem_type <- ifelse(grepl("PuS", samples_names), "Sensitive", "Tolerant")
root_type <- ifelse(grepl("PoS", samples_names), "Sensitive", "Tolerant")

prepared_data <- cbind(Sample_ID = samples_names, treatment = treatment_type,
              stem = stem_type, root = root_type, data)
              
head(prepared_data[,0:10])
```

We need to save the prepared data for the analysis:

```
save(prepared_data, file = "data/prepared_data.RData")
```

Our objective is to predict the sample classes by implementing classification algorithms, including Lasso-penalised Logistic Regression, Random Forest and Support Vector Machines (SVM).

First, we will perform a binary analysis to determine whether our model can predict treatment type based on gene expression levels, identifying the most important genes.

Next, we will assess the model's predictive power, focusing on stem type independently of treatment and root type. If successful, this will help identify the key genes involved in the drought response. Lastly, we will conduct a similar analysis considering the root type.
