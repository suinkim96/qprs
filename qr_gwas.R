# -------------------------------
# QR-based GWAS Main Script
# -------------------------------

rm(list=ls())
library(genio)
library(tidyverse)
library(quantreg)
library(data.table)
library(BEDMatrix)
source("functions.R")

cov <- 'your covariates data' 
pheno <- 'your phenotype data'

data_all <- cov %>%
  left_join(pheno, by = "IID") %>%
  select(-FID) %>%
  setDT()

geno <- BEDMatrix("data/1KG_hapmap_final2.bed", simple_names = TRUE)
qntl <- seq(0.1, 0.9, 0.1)
pheno.name <- 'Phenotype'

# Null model residuals
fml.null <- paste0(pheno.name, " ~ Sex + Age + BMI + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
res_df <- map_dfr(qntl, function(tau) {
  fit <- rq(as.formula(fml.null), data = data_all, tau = tau)
  tibble(
    Quantile = tau,
    IID = data_all$IID,
    Residual = data_all[[pheno.name]] - predict(fit, data_all)
  )
})

# GWAS loop
set.seed(2024)
bim_id <- colnames(geno)
result <- data.table()

for (idx in 1:ncol(geno)){
  df <- cbind(data_all, GENOTYPE = geno[, idx])
  fml.qr <- as.formula(paste0(pheno.name, " ~ Sex + Age + BMI + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + GENOTYPE"))
  
  beta <- map_dbl(qntl, ~ rq(fml.qr, data = df, tau = ., method = "fn")$coefficients["GENOTYPE"])
  lm_fit <- lm(fml.qr, data = df)
  lm_summary <- summary(lm_fit)
  
  p <- qrank.test(df.qr = df, fit.residual = res_df %>% filter(IID %in% df$IID), tau = qntl)
  result[, bim_id[idx] := c(lm_summary$coefficients["GENOTYPE", 4], p$cauchy.pvalue, p$quantile.pvalue, lm_summary$coefficients["GENOTYPE", 1], beta)]

  if (idx %% 100 == 0) {
    cat("Processed SNP:", idx, "\n")
    saveRDS(result, "qr_gwas.RDS")
  }
}
