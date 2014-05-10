#' @import plyr
#' @import genefilter
#' @import data.table
#' @import Biobase
#' @import limma

require(plyr)
require(genefilter)
require(Biobase)
require(limma)
require(data.table)


#' @title do.diffanal2
#' 
#' @description Factory function which performs general input cleaning and verification before passing it along to the specified strategy.
#' 
#' @param exprs A numeric \code{matrix} of gene expression data
#' @param pheno.model.frame A \code{data.frame}. Each column corresponds to a phenotypic predictor or confounder
#' @param model.formula A \code{formula} object. Should have no RHS. Corresponds to the columns in \code{pheno.model.frame}. Exclusive with the (\code{predictor}, \code{confounder}) arguments.
#' @param ngenes An integer. Controls the number of genes to report for each test class.
#' @param alpha The significance threshold value to require to report a gene
#' @param p.adjust.method The name of p-value adjustment method to use. \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"} 
#' @param one.vs.all A flag indicating whether or not to use a one-vs-all design. This is most useful for initial exploration. 
#' @param summary.transform A function to be applied to \code{exprs} to transform it before calculating summary statistics
#' @export   
do.diffanal2 <- function(exprs, pheno.model.frame, model.formula = NULL, predictor = NULL, confounders = list(), strategy = c('t.test'), ngenes = 250, alpha = 1, p.adjust.method = 'fdr', one.vs.all = T, summary.transform = unlog, unique=F, 
                      ...){
  # Gets parameter list, for later writing to file
  args <- formals()
  args$exprs <- substitute(exprs)
  args$pheno.model.frame <- substitute(pheno.model.frame)
  args$model.formula <- substitute(model.formula)
  
  dot.args <- list(...)
  
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

  model.formula.terms <- get.model.formula.terms(model.formula)
  
  if(!all(model.formula.terms %in% colnames(pheno.model.frame))){
    missing.terms <- !(model.formula.terms %in% colnames(pheno.model.frame))
    missing.terms <- as.list(model.formula.terms[missing.terms])
    err.msg <- sprintf("term %s is missing from the phenotype model frame.", missing.terms)
    stop(err.msg)
  }
  
  cleaned <- clean.samples(exprs, pheno.model.frame, model.formula.terms)
  exprs <- cleaned[["exprs"]]
  pheno.model.frame <- cleaned[["pheno"]]
  
  pheno.model.frame <- prepare.pheno.frame(pheno.model.frame, model.formula.terms)
  
  results <- NULL
  if(strategy == 't.test'){
    results.table <- diffanal.t.test(exprs=exprs, 
                                     pheno.model.frame=pheno.model.frame,                                      
                                     model.formula=model.formula, 
                                     alpha=alpha, 
                                     ngenes=ngenes,
                                     p.adjust.method=p.adjust.method, 
                                     one.vs.all=one.vs.all, 
                                     summary.transform=summary.transform,
                                     ...)
  } else if(strategy == 'limma'){
    results.table <- diffanal.limma.test(exprs=exprs, 
                                         pheno.model.frame=pheno.model.frame,
                                         model.formula=model.formula, 
                                         alpha=alpha, 
                                         ngenes=ngenes,
                                         p.adjust.method=p.adjust.method, 
                                         one.vs.all=one.vs.all,
                                         summary.transform=summary.transform,
                                         ...)
  } else if(strategy == 'perm.t'){
    results.table <- diffanal.perm.t.test(exprs=exprs, 
                                         pheno.model.frame=pheno.model.frame,
                                         model.formula=model.formula, 
                                         alpha=alpha, 
                                         ngenes=ngenes,
                                         p.adjust.method=p.adjust.method, 
                                         one.vs.all=one.vs.all,
                                         summary.transform=summary.transform,
                                         ...)    
  }
  else if(is.function(strategy)){
    results.table <- strategy(exprs=exprs, pheno.model.frame=pheno.model.frame, 
                        model.formula=model.formula, alpha=alpha,
                        p.adjust.method=p.adjust.method, one.vs.all=one.vs.all, 
                        summary.transform=summary.transform, 
                        # Pass along any extra arguments
                        dot.args)
  }
  
  if(is.null(results.table)) stop("The provided test strategy did not return any results.")

  attr(results.table, 'call') <- args

  return(results.table)
}
