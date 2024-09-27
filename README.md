# Supervised Learning Models for Transcriptomics Data

This repository contains **Machine Learning (ML) algorithms** for analysing **transcriptomics data**. The analysis aims to determine the most suitable algorithm for our data, focusing on predictive power. The data comprises thousands of genes expression levels in *Pinus pinaster* samples, with the goal of identifying biomarkers for drought resistance (note: the data cannot be published yet as it is part of an upcoming article). Since the labels of the samples are known, supervised learning models will be used. The main objective is to predict the classes of the samples, which requires the implementation of classification algorithms, including **Lasso-penalised Logistic Regression**, **Random Forest** and **Support Vector Machines (SVM)**.

Samples are characterised by three traits to predict: **Treatment type** (Control/Drought), **Stem type** (Tolerant/Sensitive) and **Root type** (Tolerant/Sensitive). Consequently, the repository is divided into the following sections:

- Section 0: ![Data Preparation](data_preparation.md)
- Section 1: ![Lasso-penalised Logistic Regression for **treatmen** type](treatment_lasso_kfoldcv.md)
- Section 2: ![Random Forest for **treatment** type](treatment_rf_kfoldcv.md)
- Section 3: ![SVM for **treatment** type](treatment_svm_kfoldcv.md)
- Section 4: ![Lasso-penalised Logistic Regression for **stem** type](stem_lasso_kfoldcv.md)
- Section 5: ![Random Forest for **stem** type](stem_rf_kfoldcv.md)
- Section 6: ![SVM for **stem** type](stem_svm_kfoldcv.md)
- Section 7: ![Lasso-penalised Logistic Regression for **root** type](root_lasso_kfoldcv.md)
- Section 8: ![Random Forest for **root** type](root_rf_kfoldcv.md)
- Section 9: ![SVM for **root** type](root_svm_kfoldcv.md)
