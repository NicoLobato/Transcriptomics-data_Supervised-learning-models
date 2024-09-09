# Lasso-penalised Logistic Regression

## Loading libraries

```
library("tidyverse")
library("limma")
library("NOISeq")
library("glmnet")
library("knitr")
library("ROCit")
library("pROC")
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

The initial analysis will focus on the treatment type. As mentioned in the previous script, the data consists of 206,574 gene expression levels from 24 individuals. This indicates that we have significantly more variables than samples, necessitating the use of a Logistic Regression model with Lasso penalisation.

## Study of the treatment type effect

Before we begin, we need to select only the columns relevant to this binary analysis, discarding the Stem and Root types:

```
data <- prepared_data %>%
  select(-c(stem,root))
rm(prepared_data)
```

We check if our data is balanced:

```
table(data$treatment) # 2 classes of 12 observations each
```

As we see, our data is perfectly balanced, containing 2 classes with 12 observations each. Due to the relatively small sample size, when splitting the data, we are left with very few samples in the test set. This was observed previously in the script "treatment_lasso.md", where the model's performance was affected and the confusion matrix results did not align with the AUC curve. To overcome these limitations, we will perform k-fold cross-validation with k=4, resulting in 4 folds of 6 samples each and providing a total of 24 predictions.

It is also important to mention that we will use the TMM + VOOM normalisation method directly, as it is the standard method in transcriptomics and has proven to be the best for our data.

```
n <- nrow(data)
k <- 4  # number of folds
folds <- sample(rep(1:k, length.out = n))  # assign each sample to a fold
cv_results <- data.frame(Actual = character(n), Predicted = character(n), stringsAsFactors = FALSE)
all_probs <- numeric(n) # vector to store predicted probabilities
```

## Model creation

```
for (fold in 1:k) {

  test_indices <- which(folds == fold) # select the indices for the current test set
  
  # Splitting the data into training and test sets
  test_set <- data[test_indices,] 
  training_set <- data[-test_indices,]  # select all rows except those in the test set
  
  training_id <- training_set$Sample_ID
  test_id <- test_set$Sample_ID
  training_set$Sample_ID <- NULL
  test_set$Sample_ID <- NULL
  
  # Preparing the data
  x_train <- as.matrix(training_set[,-1])
  y_train <- as.factor(training_set$treatment)
  x_test <- as.matrix(test_set[,-1])
  y_test <- as.factor(test_set$treatment)
  
  # Training set preprocessing
  x_train_prep <- t(x_train)
  colnames(x_train_prep) <- training_id
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
  
  predictions <- ifelse(probs > 0.5, "Drought", "Control")
  
  cv_results[test_indices, ] <- data.frame(Actual = as.character(y_test), Predicted = predictions,
                                stringsAsFactors = FALSE) # store results
}
```

Therefore, we obtain a data frame with all predictions, which allows us to create the confusion matrix:

```
table(cv_results$Predicted, cv_results$Actual) # TP = 10 ; FP = 1; FN= 2; TN = 11
```

We can generate the ROC curve:

```
ROCit_lasso <- rocit(score = all_probs, class = as.factor(cv_results$Actual))
plot(ROCit_lasso, col = c(1,"gray50"), legend = FALSE, YIndex = FALSE)
legend("bottomright", col = 1, legend = paste("Lasso (AUC =",round(ROCit_lasso$AUC,2),")"), lwd = 2)
```

We can now observe that the confusion matrix and ROC curve align well, with the accuracy of our Lasso-penalised Logistic Regression model being approximately 0.9.

```
save(ROCit_lasso, file = "data/treatment/lasso_ROC.RData")
```

We save the ROC curve to compare it with other classifiers. It is important to ensure that the same number of variables is used in the models being compared to make the comparison valid.
