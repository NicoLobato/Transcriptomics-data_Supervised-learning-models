# Random Forest

## Loading libraries

```
library("knitr")
# library("caret")       # for general model fitting
library("rpart")       # for fitting decision trees
library("ipred")  
library("foreach")     # for parallel processing with for loops
library("ggplot2")
library("rpart.plot")
library("doParallel")  # for parallel backend to foreach
library("data.table")
library("tidymodels")
library("randomForest")
library("ROCit")
library("pROC")
library("dplyr")
#library("NOISeq")
#library("limma")
```

## Reading the data

We read the prepared data in the "0.data_preparation.md":

```
load("data/prepared_data.RData")
```

This model will focus on the treatment type.

## Study of the treatment type effect

Before starting, we need to select only the columns relevant to this binary analysis, discarding Stem and Root types:

```
data <- prepared_data %>%
  select(-c(root,stem))
rm(prepared_data)
```

We check if our data is balanced:

```
table(data$treatment) # 2 classes of 12 observations each
```

As we see, our data is perfectly balanced, containing 2 classes with 12 observations each. Due to the relatively small sample size, when splitting the data, we are left with very few samples in the test set. Although the AUC results aligned with the confusion matrix, they are not reliable due to the insufficient test samples for making accurate predictions. To overcome these limitations, we will perform k-fold cross-validation with k=4, resulting in 4 folds of 6 samples each and providing a total of 24 predictions.

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
  train_set <- data[-test_indices,]  # select all rows except those in the test set
  
  train_id <- train_set$Sample_ID
  test_id <- test_set$Sample_ID
  train_set$Sample_ID <- NULL
  test_set$Sample_ID <- NULL
  
  x_train <- as.matrix(train_set[,-1])
  y_train <- as.factor(train_set$treatment)
  x_test <- as.matrix(test_set[,-1])
  y_test <- as.factor(test_set$treatment)
  
  tuned_rf <- tuneRF(x_train, y_train, mtryStart = sqrt(ncol(x_train)), ntreeTry=50, stepFactor=2, improve=0.05, trace=TRUE, plot=TRUE)
                   
  num_row <- which(tuned_rf[,2] == min(tuned_rf[,2])) # the first row that minimizes OOB error
  num_vars <- tuned_rf[num_row,1]  # 'mtry' value corresponding to min OOB error
  
  rf.treatment <- randomForest(x = x_train, y = y_train, ntree = 100, mtry = 20,
                             importance = TRUE)
                             
  print(rf.treatment)
  
  # Make predictions
  
  probs <- predict(rf.treatment, newdata = x_test, type = "prob")
  
  predictions <- ifelse(probs[,"Drought"] > 0.5, "Drought", "Control")

  all_probs[test_indices] <- probs[, "Drought"] # store predicted probabilities
  
  cv_results[test_indices, ] <- data.frame(Actual = as.character(y_test), Predicted = predictions, stringsAsFactors = FALSE) # store results
}
```

Therefore, we obtain a data frame with all predictions, which allows us to create the confusion matrix:

```
table(cv_results$Predicted, cv_results$Actual) # TP = 10 ; FP = 2; FN= 2; TN = 10
```

We can generate the ROC curve:

```
ROCit_rf <- rocit(score = all_probs, class = as.factor(cv_results$Actual))
plot(ROCit_rf, col = c(1,"gray50"), legend = FALSE, YIndex = FALSE)
legend("bottomright", col = 1, legend = paste("RF (AUC =",round(ROCit_rf$AUC,2),")"), lwd = 2)
```

We observe that our model's accuracy is 0.88.

We can now compare our model's ROC curve with Lasso model's:

```
load("data/treatment/lasso_ROC.RData")

ROCit_rf <- rocit(score=probs[,"Drought"],class=y_test)
plot(ROCit_rf, col = c(1,"gray50"), legend = FALSE, YIndex = FALSE)
lines(treatmentROClasso~treatmentROClasso,col = 2, lwd = 2)
legend("bottomright", col = c(1,2), c("RF","Lasso"), lwd = 2)
```

We save the ROC curve for next comparisons.

```
save(ROCit_rf, file = "data/treatment/rf_ROC.RData")
```
