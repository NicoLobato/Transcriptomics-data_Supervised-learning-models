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

The data contains 206,574 gene expression levels in 24 samples, named according to their drought-related characteristics as follows:
- Rootstock type: PoT (Tolerant) / PoS (Sensitive)
- Scion type: PuT (Tolerant) / PuS (Sensitive)
- Treatment type: C (Control) / S (Drought)

This data will be prepared and reformatted for further analysis as follows:

```
genes_names <- data$target_id # save the transcript IDs
samples_names <- names(data)[3:ncol(data)] # save the sample IDs

data <- as.matrix(data[,-c("V1","target_id")])
rownames(data) <- genes_names
rm(genes_names)
```

Before starting the analysis, the matrix of read counts is filtered so that I keep transcripts where CPM >= 1 in at least 3 samples (CPM = counts per million).

```
keep=rowSums(cpm(data)>=1)>=3
length(keep)
data <- data[keep,]
rm(keep)
dim(data) # 44881    24
# The number of transcripts is now significantly reduced
```

We now create the data structure:
- 1st col: Sample IDs
- 2nd col: Treatment type (Control/Drought)
- 3rd col: Scion type (Tolerant/Sensitive)
- 4th col: Rootstock type (Tolerant/Sensitive)
- Remaining cols: Gene expression levels

```
data <- t(data) %>%
  as.data.table(.)

treatment_type <- ifelse(endsWith(samples_names, "C"), "Control", "Drought")
scion_type <- ifelse(grepl("PuS", samples_names), "Sensitive", "Tolerant")
rootstock_type <- ifelse(grepl("PoS", samples_names), "Sensitive", "Tolerant")

prepared_data <- cbind(Sample_ID = samples_names, treatment = treatment_type,
              scion = scion_type, rootstock = rootstock_type, data)
              
head(prepared_data[,0:10])
```

We need to save the prepared data for the analysis:

```
save(prepared_data, file = "data/prepared_data.RData")
```

Our objective is to predict the categories by implementing classification algorithms, including Lasso-penalized Logistic Regression, Random Forest, and Support Vector Machines.

First, we will perform a binary analysis to determine whether our models can predict treatment type based on gene expression levels, regardless of scion and rootstock type, while identifying the most important variables.

Next, we will assess the models' predictive power, focusing on scion type independently of treatment and rootstock type. If successful, this will help identify the key genes involved in drought response. Lastly, we will conduct a similar analysis considering rootstock type.
