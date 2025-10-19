# Quantile Regression-based GWAS (QR-GWAS)

This repository contains R scripts implementing the **Quantile Regression GWAS** described in  
**Kim, Goo, Park & Park (2025)** – *"Enhancing Polygenic Risk Prediction by Modeling Quantile-Specific Genetic Effects"*.

## Overview
The pipeline estimates SNP-specific genetic effects across multiple phenotype quantiles (τ = 0.1–0.9), using rank-based quantile regression.  
It then aggregates multiple quantile p-values via the **Aggregated Cauchy Association Test (ACAT)** to yield robust association statistics.

## Directory Structure
- `scripts/` — Main analysis scripts and helper functions  
- `data/` — Input covariates, phenotype, and genotype files  
- `results/` — Output RDS files (GWAS summary results)

## Main Steps
1. Run quantile regression for each τ
2. Compute residuals and rank-based test statistics
3. Combine p-values via ACAT
4. Save genome-wide summary statistics

## Citation
> Park, M., Kim, S., Goo, T., Park, T. (2025).  
> *Enhancing Polygenic Risk Prediction by Modeling Quantile-Specific Genetic Effects.*
