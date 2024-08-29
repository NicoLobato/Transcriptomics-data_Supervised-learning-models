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

As we can see, our data is perfectly balanced, containing 2 classes with 12 observations each. Therefore, it is not necessary to perform stratified splitting into training and testing sets, although we will still do so.

## Splitting the data into training and test set

We now split the data into training and test sets in a stratified manner, despite having balanced data.

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

We now need to preprocess the data. We proceed in this order because we do not normally have a testing set (it does not exist until we want to apply the model), but fortunately, we have it. Consequently, we begin by dividing the data into two independent sets to avoid contamination.

Preprocessing includes data normalisation (essential in Logistic Regression) and selecting genes with higher variability (more significant), which reduces the number of variables and improves data quality.

Firstly, we preprocess the training set:

```
# Training set preprocessing
x_train_prep <- t(x_train)
colnames(x_train_prep) <- training_id
dim(x_train_prep) # 44881    18
head(x_train_prep)
str(x_train_prep)

# We apply 3 popular normalisation methods in RNA-Seq techniques: VOOM, TMM and TMM + VOOM
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

In this way, we select the 1500 most significant genes from the training set.

```
# Test set preprocessing
x_test_prep <- t(x_test)
colnames(x_test_prep) <- test_id
dim(x_test_prep) # 44881     6
head(x_test_prep)
str(x_test_prep)

# VOOM
x_test_prep_norm_voom <- voom(x_test_prep)
# TMM
x_test_prep_norm_tmm <- tmm(x_test_prep)
# TMM + VOOM
x_test_prep_norm_tmm_voom <- voom(x_test_prep_norm_tmm)

x_test_prep_df <- as.data.frame(x_test_prep)
x_test_prep_norm_voom_df <- as.data.frame(x_test_prep_norm_voom)
x_test_prep_norm_tmm_df <- as.data.frame(x_test_prep_norm_tmm)
x_test_prep_norm_tmm_voom_df <- as.data.frame(x_test_prep_norm_tmm_voom)

boxplot(x_test_prep_df, outline=FALSE, main="Before normalization", xaxt="n")
boxplot(x_test_prep_norm_voom_df, outline=FALSE, main="After normalization with VOOM", xaxt="n")
boxplot(x_test_prep_norm_tmm_df, outline=FALSE, main="After normalization with TMM", xaxt="n")
boxplot(x_test_prep_norm_tmm_voom_df, outline=FALSE, main="After normalization with TMM + VOOM", xaxt="n")

# The most suitable normalization method is TMM + VOOM again.

# We now select the 1500 genes that were chosen in the training set before:
x_test_final <- x_test_prep_norm_tmm_voom_df[genes_train,]
sum(is.na(x_test_final)) # there are not missing values
dim(x_test_final) # 1500    6
```

Once the data has been normalised and reduced, we can proceed with the model.

```
x_train_final <- t(x_train_final) # transpose to have the genes by columns
x_test_final <- t(x_test_final)
```

We then apply k-fold cross-validation to select the best value for the regularization parameter lambda (model tuning).

```
cv.lasso <- cv.glmnet(x_train_final, y_train, family = "binomial", type.measure = "class",
                      nfolds = 4, alpha = 1)
plot(cv.lasso)
print(cv.lasso$lambda.min) # the value of lambda that minimises the classification error: 0.3468004
```

Therefore, this is the best lambda value (that minimises the classification error) that will be used to fit the final model.

## Fit the final model on the training data

We can now fit the final model (tuned with cross-validation) to the training set:

```
model <- glmnet(x_train_final, y_train, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.min)
tmp_coeffs <- coef(model)
df <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
kable(df)
```

## Make predictions

With the solutions from the trained model, we can now make predictions:

```
probs <- predict(model, newx = x_test_final, type = "response")
predictions <- ifelse(probs > 0.5, "Drought", "Control")
table(predictions, y_test) # TP = 3 ; FP = 1; FN= 0; TN = 2
```

We can implement ROC curves, which can help us evaluate different classifiers over all possible decision thresholds:

```
ROCit_lasso <- rocit(score=probs[,1],class=y_test)
plot(ROCit_lasso, col = c(1,"gray50"), legend = FALSE, YIndex = FALSE)
legend("bottomright", col = 1, legend = paste("Lasso (AUC =",round(ROCit_lasso$AUC,2),")"), lwd = 2) # AUC = accuracy = 1
```

We can observe that the results obtained from the confusion matrix do not align with the AUC curve (Area Under the Curve). This may be due to the relatively small sample size, such that the split of the data results in a test set with insufficient samples, leading to discrepancies in the model’s predictions. Therefore, in the following script, we will develop a new Logistic Regression model to address this limitation.
