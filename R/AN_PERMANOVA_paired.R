
# Multivariate Repeated Measurements with adonis(): 
# https://thebiobucket.blogspot.com/2011/04/repeat-measure-adonis-lately-i-had-to.html



adonis_blocks <- function(x, condition, permutation = 999, method = "jaccard", blocks, ...){
  formula <- as.formula("x ~ condition")
  fit <- adonis(formula = formula, permutations = 1, method = method, ...)
  pop <- rep(NA, permutation + 1)
  pop[1] <- fit$aov.tab[1, 5]
  ctrl <- how(blocks = blocks, within = Within(type = "series", mirror = FALSE))
  nobs <- nrow(x)
  formula_rand <- update(formula, . ~ .[idx])
  for(i in 2:(permutation + 1)) {
    idx <- shuffle(nobs, control = ctrl)
    fit_rand <- adonis(formula = formula_rand, permutations = 1, method = method, ...)
    pop[i] <- fit_rand$aov.tab[1, 5]
  }
  pval <- sum(pop >= pop[1]) / (permutation + 1)
  fit$aov.tab["condition", "Pr(>F)"] <- pval
  return(fit)
}




