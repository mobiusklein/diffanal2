do.limma.lmFit <- function(exprs, pheno, model.formula){
  design <- model.matrix(model.formula, pheno)
  fit <- lmFit(exprs, design)
  fit.ebayes <- eBayes(fit)
  return(fit.ebayes)
}