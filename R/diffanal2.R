require(plyr)
require(genefilter)
require(Biobase)
require(limma)
require(data.table)

# Column names with aggregate meaning will have tokens joined by the JOIN.CHR
# e.g. cls1-cls2-p.value means the p.value comparing cls1 and cls2.
JOIN.CHR = "-"
SAMPLE.ID.COL <- '.rownames'

#' 
#' @param exprs numeric matrix of gene expression data
#' @param pheno.model.frame data.frame. Each column corresponds to a phenotypic predictor or confounder
#' @param model.formula should have no RHS. Corresponds to the columns in \code{pheno.model.frame}
#' @param ngenes controls the number of genes to report
#' @param alpha significance threshold value
#' @param p.adjust.method name of p value adjustment method to use. \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"} 
#' @param one.vs.all 
#'   
diffanal2 <- function(exprs, pheno.model.frame, model.formula = NULL, predictor = NULL, confouders = list(), strategy = c('t.test'), ngenes = 250, alpha = 0.05, p.adjust.method = 'fdr', one.vs.all = F, fold.change.transform.fn = unlog, unique=F){
  # Gets parameter list, for later writing to file
  args <- formals()
  args$exprs <- substitute(exprs)
  args$pheno.model.frame <- substitute(pheno.model.frame)
  args$model.formula <- substitute(model.formula)
  #print(args)
  
  ### 
  # Validate model specification 
  if(is.null(model.formula) && is.null(predictor)){
    stop("Either `model.formula` or `predictor` must have a value")
  }
  else if(!is.null(model.formula) && !is.null(predictor)){
    stop("Only one of `model.formula` or `predictor` may have a value, but not both.")
  }
  else if(!is.null(predictor)){
    model.formula <- formula.from.list(predictor, confounders)
  }
  
  model.formula.terms <- terms(model.formula)    
  model.formula.cols <- row.names(attr(model.formula.terms, 'factors'))
  if(!all(model.formula.cols %in% colnames(pheno.model.frame))){
    missing.terms <- !(model.formula.cols %in% colnames(pheno.model.frame))
    missing.terms <- as.list(model.formula.cols[missing.terms])
    err.msg <- sprintf("term %s is missing from the phenotype model frame.", missing.terms)
    stop(err.msg)
  }
  
  cleaned <- clean.samples(exprs, pheno.model.frame, model.formula.cols)
  exprs <- cleaned[["exprs"]]
  pheno.model.frame <- cleaned[["pheno"]]
  
  pheno.model.frame <- prepare.pheno.frame(pheno.model.frame, model.formula.cols)
  
  results <- NULL
  if(strategy == 't.test'){
    results <- diffanal.t.test(exprs=exprs, pheno.model.frame=pheno.model.frame, 
                               model.formula=model.formula, alpha=alpha, 
                               p.adjust.method=p.adjust.method, one.vs.all=one.vs.all, 
                               fold.change.transform.fn=fold.change.transform.fn)
  }
  
  results.table <- do.call(cbind, results)
  # Include probe names in a separate column as row.names are not preserved by all operations
  results.table <- cbind(row.names(results.table), results.table)
  colnames(results.table) <- c("Probe", unlist(lapply(results, colnames)))
  # Sort by the first p.value column. This may not be general enough
  results.table <- results.table[order(results.table[,2]),]
  results.table <- results.table[1:ngenes,]
  results.table <- diffanal.results.table(results.table, ngenes=ngenes, 
                                        unique=unique, column.classes=lapply(results, colnames),
                                        call=args)
  return(results.table)
}
