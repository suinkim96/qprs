# -------------------------------
# Quantile Regression GWAS Helpers
# -------------------------------

# Aggregated Cauchy Association Test (ACAT)
acat <- function(pvals) {
  pvals <- pvals[!is.na(pvals)]
  if (length(pvals) == 0) return(NA)
  cauchy <- tanpi(0.5 - pvals)
  stat <- mean(cauchy)
  return(pcauchy(stat, lower.tail = FALSE))
}

# Rank-based Quantile Test
qrank.test <- function(df.qr, fit.residual, tau = c(0.25, 0.5, 0.75)) {
  ltau <- length(tau)
  VN <- outer(tau, tau, Vectorize(function(a,b) min(a,b) - a*b))
  xstar <- lm(GENOTYPE ~ Sex + Age + PC1 + PC2 + PC3 + PC4 +
              PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = df.qr)$residuals

  SN <- sapply(tau, function(t) {
    dual <- (fit.residual %>% filter(Quantile == t) %>% pull(Residual)) > 0
    ranks <- dual - (1 - t)
    as.numeric(t(xstar) %*% ranks)
  })

  VN2 <- matrix(outer(VN, t(xstar) %*% xstar, "*"), nrow = ltau)
  pvals <- pchisq(SN^2 / diag(VN2), df = 1, lower.tail = FALSE)
  list(
    tau = tau,
    quantile.pvalue = pvals,
    cauchy.pvalue = acat(pvals)
  )
}
