# Supervised Learning Models for Transcriptomics Data

This repository contains **Machine Learning (ML) algorithms** for analysing **transcriptomics data**. The analysis aims to determine the most suitable algorithm for our data, focusing on predictive power. The data comprises thousands of genes expression levels in *Pinus pinaster* samples, with the goal of identifying biomarkers for drought resistance (note: the data cannot be published yet as it is part of an upcoming article). Since the labels of the samples are known, supervised learning models will be used. The main objective is to predict the classes of the samples, which requires the implementation of classification algorithms, including **Lasso-penalised Logistic Regression**, **Random Forest** and **Support Vector Machine (SVM)**.

Samples are characterised by three traits to predict: **Treatment type** (Control/Drought), **Root type** (Tolerant/Sensitive) and **Stem type** (Tolerant/Sensitive). Consequently, the repository is divided into the following sections:

- Section 0: ![Data Preparation](data_preparation.md)

## Treatment type
- Section 1: ![Lasso-penalised Logistic Regression](treatment_lasso_kfold_cv.md)
- Section 2: ![Random Forest](treatment_rf_kfold_cv.md)
- Section 3: SVM

## Root type
- Section 1: ![Lasso-penalised Logistic Regression](root_lasso_kfold_cv.md)
- Section 2: Random Forest
- Section 3: SVM

- ## Stem type
- Section 1: Lasso-penalised Logistic Regression
- Section 2: Random Forest
- Section 3: SVM
