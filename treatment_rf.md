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
library("NOISeq")
library("limma")
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
  select(-c(stem,root))
rm(prepared_data)
```

## Splitting the data into training and test set

We now split the data into the training and test sets in a stratified way.

Although the Random Forest model does not require data normalization, we will apply it to ensure our models are more comparable.

```
seed = 17
set.seed(seed)
n = nrow(data) # sample size

training_set <- data %>%
  group_by(treatment) %>%
  sample_frac(size = 0.75)

test_recs <- !(data$Sample_ID %in% training_set$Sample_ID)
test_set <- data[test_recs,]
rm(data)

training_id <- training_set$Sample_ID
test_id <- test_set$Sample_ID
training_set$Sample_ID <- NULL
test_set$Sample_ID <- NULL

x_train <- as.matrix(training_set[,-1])
y_train <- as.factor(training_set$treatment) 
x_test <- as.matrix(test_set[,-1])
y_test <- as.factor(test_set$treatment)

table(y_train)
table(y_test)
```

## Preprocessing

We now have to preprocess the data. We proceed in this order because we do not normally have a testing set (it does not exist until we want to apply the model), but we have it fortunately. Consequently, we begin dividing the data into two independent sets in order to not contamining with information each other.

Preprocessing includes a data normalization (essential in Logistic Regression) and a selection of genes with more variability (more significative), that gets reducing the number of variables and improving the data quality.

In the first place, we preprocess the training set:

```
# Training set preprocessing
x_train_prep <- t(x_train)
colnames(x_train_prep) <- training_id
dim(x_train_prep) # 44881    18
head(x_train_prep)
str(x_train_prep)

# We apply 3 popular normalization methods in RNA-Seq techniques: VOOM, TMM and TMM + VOOM
# VOOM
x_train_prep_norm_voom <- voom(x_train_prep)
# TMM
x_train_prep_norm_tmm <- tmm(x_train_prep)
# TMM + VOOM
x_train_prep_norm_tmm_voom <- voom(x_train_prep_norm_tmm)

# We convert the previous objects into data.frame to boxplot them before and after normalization:
x_train_prep_df <- as.data.frame(x_train_prep)
x_train_prep_norm_voom_df <- as.data.frame(x_train_prep_norm_voom)
x_train_prep_norm_tmm_df <- as.data.frame(x_train_prep_norm_tmm)
x_train_prep_norm_tmm_voom_df <- as.data.frame(x_train_prep_norm_tmm_voom)

# We plot and compare them:
boxplot(x_train_prep_df, outline=FALSE, main="Before normalization", xaxt="n")
boxplot(x_train_prep_norm_voom_df, outline=FALSE, main="After normalization with VOOM", xaxt="n")
boxplot(x_train_prep_norm_tmm_df, outline=FALSE, main="After normalization with TMM", xaxt="n")
boxplot(x_train_prep_norm_tmm_voom_df, outline=FALSE, main="After normalization with TMM + VOOM", xaxt="n")

# We select TMM + VOOM normalization because is the most suitable for our data.

# We now choose the 1500 genes with more variance
varianza <- apply(x_train_prep_norm_tmm_voom_df, 1, stats::var)
varianza <- sort(varianza, decreasing=TRUE) # decreasing = TRUE to select the 1500 genes with more variance
milquinientosgenes <- varianza[1:1500]
genes_train <- names(milquinientosgenes)
x_train_final <- x_train_prep_norm_tmm_voom_df[genes_train,]
dim(x_train_final) # 1500   18
```

In that way, we select the 1500 genes more significative of training set.

```
# Test set preprocessing
x_test_prep <- t(x_test)
colnames(x_test_prep) <- test_id
dim(x_test_prep) # 44881     6
head(x_test_prep)
str(x_test_prep)

# We apply 3 popular normalization methods in RNA-Seq techniques: VOOM, TMM and TMM + VOOM
# VOOM
x_test_prep_norm_voom <- voom(x_test_prep)
# TMM
x_test_prep_norm_tmm <- tmm(x_test_prep)
# TMM + VOOM
x_test_prep_norm_tmm_voom <- voom(x_test_prep_norm_tmm)

# We convert the previous objects into data.frame to boxplot before and after normalization:
x_test_prep_df <- as.data.frame(x_test_prep)
x_test_prep_norm_voom_df <- as.data.frame(x_test_prep_norm_voom)
x_test_prep_norm_tmm_df <- as.data.frame(x_test_prep_norm_tmm)
x_test_prep_norm_tmm_voom_df <- as.data.frame(x_test_prep_norm_tmm_voom)

# We plot and compare them:
boxplot(x_test_prep_df, outline=FALSE, main="Before normalization", xaxt="n")
boxplot(x_test_prep_norm_voom_df, outline=FALSE, main="After normalization with VOOM", xaxt="n")
boxplot(x_test_prep_norm_tmm_df, outline=FALSE, main="After normalization with TMM", xaxt="n")
boxplot(x_test_prep_norm_tmm_voom_df, outline=FALSE, main="After normalization with TMM + VOOM", xaxt="n")

# The most suitable normalization method is TMM + VOOM again.

# We now choose the 1500 genes that were selected in the training set before:
x_test_final <- x_test_prep_norm_tmm_voom_df[genes_train,]
sum(is.na(x_test_final)) # there are not missing values
dim(x_test_final) # 1500    6
```

Once the data has been normalized and reduced, we can proceed with the model.

```
x_train_final <- t(x_train_final) # transponing to have the genes by columns
x_test_final <- t(x_test_final)
```

## Tuning hyperparameters

Recall that Random Forest not only does Bagging, which involves bootstrapping or sampling with replacement of observations to reduce variance, but also does it with the features. Therefore, we use the function ‘tuneRF’ to choose the best value of the hyperparameter 'mtry', which corresponds to the number of variables to be sampled in each tree. This function includes important parameters such as 'mtryStart' (initial value of 'mtry', which is usually the square root of the number of variables for classification problems), 'ntreeTry' (number of trees used in the tuning), 'stepFactor' (value that increases or decreases 'mtry') and 'improve' (minimum improvement value that the OOB error must exceed to continue the search for the optimal 'mtry').

```
tuned_rf <- tuneRF(x_train_final, y_train, mtryStart = sqrt(ncol(x_train_final)), ntreeTry=50,
                   stepFactor=2, improve=0.01, trace=TRUE, plot=TRUE)
head(tuned_rf)
# We observe that the lowest OOB error is obtained using 38.73 variables in each tree.

num_row <- which(tuned_rf[,2] == min(tuned_rf[,2])) # the first row that minimizes OOB error
num_vars <- tuned_rf[num_row,1]  # 'mtry' value corresponding to min OOB error
print(num_vars)
```

## Fitting the model

WE fit the model using the 'mtry' value that minimises the OOB error:

```
rf.treatment <- randomForest(x = x_train_final, y = y_train, ntree = 100, mtry = num_vars,
                             importance = TRUE)
rf.treatment
```

We obtain a relatively good confusion matrix with an error rate of 11.11%.

## Making predictions

We are now ready to evaluate our tuned RF model on the test data:

```
probs <- predict(rf.treatment, newdata = x_test_final, type = "prob")
predictions <- ifelse(probs[,"Drought"] > 0.5, "Drought", "Control")
table(predictions, y_test) # TP = 3 ; FP = 0; FN= 0; TN = 3
```

We can implement ROC curves:

```
ROCit_rf <- rocit(score=probs[,"Drought"],class=y_test)
plot(ROCit_rf, col = c(1,"gray50"), legend = FALSE, YIndex = FALSE)
legend("bottomright", col = 1, "RF", lwd = 2)

roc(y_test, probs[,"Drought"]) %>% print(auc(.)) # AUC = accuracy = 1
```

Examining the confusion matrix, we saw that our model could predict the 100% of the cases correctly, that we also see in the ROC curve. Moreover, we can determine the AUC value (Area Under Curve) where we see that the accuracy of our Lasso-penalised Logistic Regression model is 1.

However, we can observe again that results are not realistic due to the relatively small sample size. In the following script, we will develop a new Random Forest model to address this limitation.
