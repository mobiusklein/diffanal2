
# Converts a list of terms into a formula object for building models with
formula.from.list <- function(predictor, confounders = list()){
  # Attach the predictor and tilde first
  relation <- paste("~", predictor)
  # Combine any confounders into a chain of additions 
  trailing.terms <- do.call(paste, c(confounders, sep = " + "))
  # If there are trailing terms to merge
  if(length(trailing.terms) > 0){
    # Attach the predictor and the confounder string
    relation <- paste(relation, trailing.terms, sep=" + ")
  }
  full.formula <- as.formula(relation)
  return(full.formula)
}
