# Lasso-penalised Logistic Regression for Root Type

## Loading libraries

```
library("tidyverse")
library("caret")
library("NOISeq")
library("limma")
library("glmnet")
library("ROCit")
```

## Reading the data

We load the prepared data from the "0.data_preparation.md" file:

```
load("data/prepared_data.RData")
```

## Exploratory data analysis

```
dim(prepared_data) # 24 44885
head(prepared_data[,0:10])
str(prepared_data)
```

This analysis will focus on the root type. As we see, the data consists of 206,574 gene expression levels from 24 individuals. This indicates that we have significantly more variables than samples, necessitating the use of a Logistic Regression model with Lasso penalisation.

## Study of the root type effect

Before starting, we need to select only the columns relevant to this binary analysis, discarding Treatment and Stem types:

```
data <- prepared_data %>%
  select(-c(treatment,stem))
rm(prepared_data)
```

We check if our data is balanced:

```
table(data$root) # 2 classes of 12 observations each
```

As we can see, our data is perfectly balanced, containing 2 classes with 12 observations each. This small sample size results in insufficient samples in the test set after splitting the data, leading to prediction errors. This issue was previously observed in the script "treatment_lasso.md", where the model's performance was affected. To overcome these limitations, we will perform k-fold cross-validation with k=4, resulting in 4 folds of 6 samples each and providing a total of 24 predictions.

```
n <- nrow(data)
k <- 4  # number of folds
folds <- createFolds(data$root, k = k, list = TRUE, returnTrain = FALSE) # create stratified folds indices
cv_results <- data.frame(Actual = character(n), Predicted = character(n), stringsAsFactors = FALSE)
all_probs <- numeric(n) # vector to store predicted probabilities
```

## Model creation

It is also important to mention that during the preprocessing steps (data normalisation is essential for Logistic Regression), we will use the TMM + VOOM normalisation method directly, as it is the standard approach in transcriptomics and has been shown to be the most effective for our data.

```
for (fold in 1:k) {

  test_indices <- folds[[fold]] # select the indices for the current test set
  
  # Splitting the data into training and test sets
  test_set <- data[test_indices, ]
  train_set <- data[-test_indices, ] # select all rows except those in the test set
  
  train_id <- train_set$Sample_ID
  test_id <- test_set$Sample_ID
  train_set$Sample_ID <- NULL
  test_set$Sample_ID <- NULL
  
  x_train <- as.matrix(train_set[,-1])
  y_train <- as.factor(train_set$root)
  x_test <- as.matrix(test_set[,-1])
  y_test <- as.factor(test_set$root)
  
  # Training set preprocessing
  x_train_prep <- t(x_train)
  colnames(x_train_prep) <- train_id
  x_train_prep_norm_tmm <- tmm(x_train_prep)
  x_train_prep_norm_tmm_voom <- voom(x_train_prep_norm_tmm)
  x_train_prep_norm_tmm_voom_df <- as.data.frame(x_train_prep_norm_tmm_voom)
  
  varianza <- apply(x_train_prep_norm_tmm_voom_df, 1, stats::var)
  varianza <- sort(varianza, decreasing=TRUE) # decreasing = TRUE to select the 1500 genes with more variance
  milquinientosgenes <- varianza[1:1500]
  genes_train <- names(milquinientosgenes)
  x_train_final <- x_train_prep_norm_tmm_voom_df[genes_train,]
  
  # Test set preprocessing
  x_test_prep <- t(x_test)
  colnames(x_test_prep) <- test_id
  x_test_prep_norm_tmm <- tmm(x_test_prep)
  x_test_prep_norm_tmm_voom <- voom(x_test_prep_norm_tmm)
  x_test_prep_norm_tmm_voom_df <- as.data.frame(x_test_prep_norm_tmm_voom)
  
  x_test_final <- x_test_prep_norm_tmm_voom_df[genes_train,]
  sum(is.na(x_test_final)) # ensure there are no missing values
  
  x_train_final <- t(x_train_final) # transpose to have the genes by columns
  x_test_final <- t(x_test_final)
  
  # Fit the model
  cv.lasso <- cv.glmnet(x_train_final, y_train, family = "binomial", type.measure = "class", alpha = 1)
  
  model <- glmnet(x_train_final, y_train, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.min)
  
  # Make predictions
  probs <- predict(model, newx = x_test_final, type = "response")
  
  all_probs[test_indices] <- probs # store predicted probabilities

  predictions <- ifelse(probs > 0.5, "Tolerant", "Sensitive")

  cv_results[test_indices, ] <- data.frame(Actual = as.character(y_test), Predicted = predictions, stringsAsFactors = FALSE) # store results
}
```

We obtain a data frame with all predictions, which allows us to create the confusion matrix:

```
table(cv_results$Predicted, cv_results$Actual)
```

From the confusion matrix, we achieve a perfect accuracy of 1.

Then, we can generate the ROC curve:

```
ROCit_lasso <- rocit(score = all_probs, class = as.factor(cv_results$Actual))
plot(ROCit_lasso, col = c(1,"gray50"), legend = FALSE, YIndex = FALSE)
legend("bottomright", col = 1, legend = paste("Lasso (AUC =",round(ROCit_lasso$AUC,2),")"), lwd = 2)
```

We observe that our model's AUC value is also perfect, being 1.

```
save(ROCit_lasso, file = "results/root/lasso_ROC.RData")
```

We save the ROC curve for comparison with other classifiers. It is crucial to ensure that the same number of variables is used in the models being compared to maintain the validity of the comparison.
