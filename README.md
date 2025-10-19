# Quantile Regression-based GWAS (QR-GWAS)

This repository contains R scripts implementing the **Quantile Regression GWAS** described in  
**Kim, Goo, Park & Park (2025)** – *"Enhancing Polygenic Risk Prediction by Modeling Quantile-Specific Genetic Effects"*.

## Overview
Traditional GWAS estimates the average SNP effect on a phenotype (mean-based OLS).
QR-GWAS generalizes this by modeling SNP effects across multiple conditional quantiles (τ = 0.1–0.9) of the phenotype distribution.
From these quantile-specific effects, we derive Quantile Polygenic Risk Scores (QPRS) that capture heterogeneous genetic influences across the outcome distribution.

## Main Steps
1. Run quantile regression for each τ
2. Compute residuals and rank-based test statistics
3. Combine p-values via ACAT
4. Save genome-wide summary statistics
5. Perform LD clumping and p-value thresholding (C+T)
6. Calculate QPRS with the results from (4)-(5).

## Citation
> Park, M., Kim, S., Goo, T., Park, T. (2025).  
> *Enhancing Polygenic Risk Prediction by Modeling Quantile-Specific Genetic Effects.*
