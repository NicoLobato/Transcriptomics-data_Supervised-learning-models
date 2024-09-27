# Support Vector Machines (SVM)

## Loading libraries

```
library("tidyverse")
library(caret)
library("NOISeq")
library("limma")
library("e1071")
library("ROCit")
```

## Reading the data

We load the prepared data from the "0.data_preparation.md" file:

```
load("data/prepared_data.RData")
```

This model will focus on the treatment type.

## Study of the treatment type effect

Before starting, we need to select only the columns relevant to this binary analysis, discarding Stem and Root types:

```
data <- prepared_data %>%
  select(-c(stem,root))
rm(prepared_data)
```

We check if our data is balanced:

```
table(data$treatment) # 2 classes of 12 observations each
```

As we see, our data is perfectly balanced, containing 2 classes with 12 observations each. 

Once again, we will proceed with k-fold cross-validation with k=4 to overcome the limitation of the small sample size.

```
n <- nrow(data)
k <- 4  # number of folds
folds <- createFolds(data$treatment, k = k, list = TRUE, returnTrain = FALSE) # create stratified folds indices
cv_results <- data.frame(Actual = character(n), Predicted = character(n), stringsAsFactors = FALSE)
all_probs <- numeric(n) # vector to store predicted probabilities
```

## Model creation

We must remember that SVM also requires data normalisation.

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
  y_train <- as.factor(train_set$treatment)
  x_test <- as.matrix(test_set[,-1])
  y_test <- as.factor(test_set$treatment)
  
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
  
  # SVM model
  x_train_final <- as.data.frame(x_train_final)
  x_test_final <- as.data.frame(x_test_final)
  x_train_final$treatment <- y_train

  cost_values = c(0.001, 0.01, 0.1, 1, 5, 10, 100)
  
  # Hyperparameter tuning: selecting the optimal cost value
  tune_out = tune(svm, factor(treatment) ~ ., data = x_train_final, kernel = "radial", ranges = list(cost = cost_values), probability = TRUE)
  
  best_model = tune_out$best.model
  
  # Make predictions using the best model
  probs <- predict(best_model, newdata = x_test_final, probability=TRUE)
  probs <- attr(probs, "probabilities") # select the probabilities
  
  predictions <- ifelse(probs[,"Drought"] > 0.5, "Drought", "Control")

  all_probs[test_indices] <- probs[,"Drought"] # store predicted probabilities
  
  cv_results[test_indices, ] <- data.frame(Actual = as.character(y_test), Predicted = predictions, stringsAsFactors = FALSE) # store results
}
```

We obtain a data frame with all predictions, which allows us to create the confusion matrix:

```
table(cv_results$Predicted, cv_results$Actual)
```

From the confusion matrix, we obtain an average accuracy of approximately 0.75.

Then, we can generate the ROC curve:

```
ROCit_svm <- rocit(score = all_probs, class = as.factor(cv_results$Actual))
plot(ROCit_svm, col = c(1,"gray50"), legend = FALSE, YIndex = FALSE)
legend("bottomright", col = 1, legend = paste("SVM (AUC =",round(ROCit_svm$AUC,2),")"), lwd = 2)
```

We observe that the average AUC value of our model is approximately 0.8.

We can now compare our model's ROC curve with those of the Lasso and Random Forest models:

```
load("results/treatment/lasso_ROC.RData")
load("results/treatment/rf_ROC.RData")

ROCit_svm <- rocit(score = all_probs, class = as.factor(cv_results$Actual))
plot(ROCit_svm, col = c(1,"gray50"), legend = FALSE, YIndex = FALSE)
lines(ROCit_lasso$TPR~ROCit_lasso$FPR, col = 2, lwd = 2)
lines(ROCit_rf$TPR~ROCit_rf$FPR, col = 3, lwd = 2)
legend("bottomright", col = c(1,2,3), c("SVM","Lasso","RF"), lwd = 2)
```

We save the ROC curve.

```
save(ROCit_svm, file = "results/treatment/svm_ROC.RData")
```