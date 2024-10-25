# Supervised Learning Models for Transcriptomics Data

This repository contains **Machine Learning (ML) algorithms** for analyzing **transcriptomics data** from *Pinus pinaster* Ait. The analysis focuses on identifying the most suitable algorithm for predicting gene biomarkers related to drought resistance. The dataset includes expression levels of thousands of genes across samples (note: the data are available on request from the corresponding author because the data are part of an ongoing study).

Using known sample labels, we will implement supervised learning models, specifically **Lasso-penalized Logistic Regression**, **Random Forest (RF)**, and **Support Vector Machines (SVM)**, to predict the following categories: **Treatment type** (Drought/Control), **Scion type** (Tolerant/Sensitive), and **Rootstock type** (Tolerant/Sensitive). If the models effectively predict these labels, variable importance analyses will be performed to identify reliable biomarkers for drought tolerance (only Logistic Regression and RF yielded trustworthy results). The repository is organized into the following sections:

- Section 0: ![Data Preparation](data_preparation.md)![ and Exploratory Data Analysis](treatment_lasso.md)
- Section 1: ![Lasso-penalized Logistic Regression for Treatment Type](treatment_lasso_kfoldcv.md)
- Section 1.1: ![Variable Importance Analysis of Lasso-penalized Logistic Regression for Treatment Type](treatment_lasso_var_imp.md)
- Section 2: ![Random Forest for Treatment Type](treatment_rf_kfoldcv.md)
- Section 2.1: ![Variable Importance Analysis of Random Forest for Treatment Type](treatment_rf_var_imp.md)
- Section 3: ![Support Vector Machines for Treatment Type](treatment_svm_kfoldcv.md)
- Section 4: ![Lasso-penalized Logistic Regression for Scion Type](scion_lasso_kfoldcv.md)
- Section 5: ![Random Forest for Scion Type](scion_rf_kfoldcv.md)
- Section 6: ![Support Vector Machines for Scion Type](scion_svm_kfoldcv.md)
- Section 7: ![Lasso-penalized Logistic Regression for Rootstock Type](rootstock_lasso_kfoldcv.md)
- Section 7.1: ![Variable Importance Analysis of Lasso-penalized Logistic Regression for Rootstock Type](rootstock_lasso_var_imp.md)
- Section 8: ![Random Forest for Rootstock Type](rootstock_rf_kfoldcv.md)
- Section 8.1: ![Variable Importance Analysis of Random Forest for Rootstock Type](rootstock_rf_var_imp.md)
- Section 9: ![Support Vector Machines for Rootstock Type](rootstock_svm_kfoldcv.md)
