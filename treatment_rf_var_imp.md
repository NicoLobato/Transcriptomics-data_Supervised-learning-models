# Variable Importance from Random Forest model for Treatment Type

## Loading libraries

```
library("tidyverse")
library("caret")
library("NOISeq")
library("limma")
library("randomForest")
```

## Reading the data

We load the prepared data:

```
load("data/prepared_data.RData")
```

We discard the Scion and Rootstock types:

```
data <- prepared_data %>%
  select(-c(scion,rootstock))
rm(prepared_data)
```

We perform k-fold cross-validation with k=4, resulting in 4 folds of 6 samples each. To study the importance of the variables, we will combine the important variables identified in each fold.

```
n <- nrow(data)
k <- 4  # number of folds
folds <- createFolds(data$treatment, k = k, list = TRUE, returnTrain = FALSE) # create stratified folds indices
importance_list <- list() # list to store variable importance
```

## Model creation

We implement the model including the study of variable importance.

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
  
  # RF model
  rf.treatment <- randomForest(x = x_train_final, y = y_train, ntree = 500, mtry = sqrt(ncol(x_train_final)), importance = TRUE) # we use this standard value of 'mtry' for classification problems
  
  # Variable importance
  imp <- as.data.frame(importance(rf.treatment)) %>%
  select(MeanDecreaseGini) %>%
  rownames_to_column(var = "Variable") %>%
  filter(MeanDecreaseGini > 0)  # select variables with 'MeanDecreaseGini' greater than 0
  
  importance_list[[fold]] <- imp  # store importance from the current fold
}
```

We obtain the most important variables from each fold, which we will combine for representation.

```
importance_combined <- do.call(rbind, importance_list) %>%
  group_by(Variable) %>% # group the same variables
  summarise(
    Importance = sum(MeanDecreaseGini) / n(), # weighted average by the frequency of occurence
    .groups = 'drop' # remove the groupings
  ) %>%
  slice_max(order_by = Importance, n = 10) %>% # select and order the top 10 genes in descending order
  mutate(Variable = factor(Variable, levels = rev(unique(Variable))))
  
var_imp <- ggplot(importance_combined, aes(x = Importance, y = Variable)) +
  geom_bar(stat="identity") +
  labs(y = NULL) +
  theme_minimal()

var_imp
```

We can now compare the top 10 most important genes obtained through the Random Forest model with those obtained through the Logistic Regression model:

```
rf_imp_genes <- as.character(importance_combined$Variable)
save(rf_imp_genes, file = "results/treatment/rf_imp_genes.RData")

load("results/treatment/lasso_imp_genes.RData")
#load("results/treatment/rf_imp_genes.RData")

intersect(rf_imp_genes,lasso_imp_genes)
```
