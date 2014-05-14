#' @import plyr
#' @import genefilter
#' @import Biobase
#' @import limma

require(plyr)
require(genefilter)
require(Biobase)
require(limma)

#' @title .do.diffanal2
#' @name .do.diffanal2
#' @description Factory function which performs general input cleaning and verification before passing it along to the specified strategy.
#' @details Performs class-comparison differential gene expression using the strategy provided. Returns a table with summary statistics about each class. 
#' @param exprs A numeric \code{matrix} of gene expression data
#' @param pheno.model.frame A \code{data.frame}. Each column corresponds to a phenotypic predictor or confounder
#' @param model.formula A \code{formula} object. Should have no RHS. Corresponds to the columns in \code{pheno.model.frame}. Exclusive with the (\code{predictor}, \code{confounder}) arguments.
#' @param ngenes An integer. Controls the number of genes to report for each test class.
#' @param p.adjust.method The name of p-value adjustment method to use. \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"} 
#' @param one.vs.all A flag indicating whether or not to use a one-vs-all design. This is most useful for initial exploration. 
#' @param summary.transform A function to be applied to \code{exprs} to transform it before calculating summary statistics
#' @return A \code{diffanal.results} instance
#' @seealso \link{diffanal.results}, \link{do.diffanal2}
.do.diffanal2 <- function(exprs, pheno.model.frame, model.formula = NULL, predictor = NULL, confounders = list(), strategy = c('t.test', 'limma', 'perm.t'), ngenes = 250, p.adjust.method = 'fdr', one.vs.all = T, summary.transform = unlog, ...){
  ### 
  # Validate model specification 
  if(is.null(model.formula) && is.null(predictor)){
    stop("Either `model.formula` or `predictor` must have a value")
  }
  else if(!is.null(model.formula) && !is.null(predictor)){
    verbose("model.formula", model.formula, "predictor", predictor)
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
  args <- list()
  args$model.formula <- substitute(model.formula)
  args$strategy <- strategy
  
  cleaned <- clean.samples(exprs, pheno.model.frame, model.formula.terms)
  exprs <- cleaned[["exprs"]]
  pheno.model.frame <- cleaned[["pheno"]]
  
  pheno.model.frame <- prepare.pheno.frame(pheno.model.frame, model.formula.terms)
  
  results <- NULL
  if(strategy == 't.test'){
    results.table <- diffanal.t.test(exprs=exprs, 
                                     pheno.model.frame=pheno.model.frame,                                      
                                     model.formula=model.formula, 
                                     ngenes=ngenes,
                                     p.adjust.method=p.adjust.method, 
                                     one.vs.all=one.vs.all, 
                                     summary.transform=summary.transform,
                                     ...)
  } 
  else if(strategy == 'limma'){
    results.table <- diffanal.limma.test(exprs=exprs, 
                                         pheno.model.frame=pheno.model.frame,
                                         model.formula=model.formula, 
                                         ngenes=ngenes,
                                         p.adjust.method=p.adjust.method, 
                                         one.vs.all=one.vs.all,
                                         summary.transform=summary.transform,
                                         ...)
  } 
  else if(strategy == 'perm.t'){
    results.table <- diffanal.perm.t.test(exprs=exprs, 
                                         pheno.model.frame=pheno.model.frame,
                                         model.formula=model.formula, 
                                         ngenes=ngenes,
                                         p.adjust.method=p.adjust.method, 
                                         one.vs.all=one.vs.all,
                                         summary.transform=summary.transform,
                                         ...)    
  }
  else if(is.function(strategy)){
    results.table <- strategy(exprs=exprs, pheno.model.frame=pheno.model.frame, 
                        model.formula=model.formula,
                        p.adjust.method=p.adjust.method, 
                        one.vs.all=one.vs.all, 
                        summary.transform=summary.transform, 
                        # Pass along any extra arguments
                        ...)
  }
  
  if(is.null(results.table)) stop("The provided test strategy did not return any results.")

  attr(results.table, 'call') <- args

  return(results.table)
}

#' @title do.diffanal2
#' @export
#' @param object An R Object
#' @param ... Additional parameters to pass to the strategy
do.diffanal2 <- function(object, ...){
  UseMethod("do.diffanal2")
}

#' @rdname do.diffanal2
#' @title do.diffanal2.default
#' @description Assumes that the expression data and pheno.model.frame are given independently. Takes the parameters as is and passes them along to the general function. 
#' @param exprs A numeric \code{matrix} of gene expression data
#' @param pheno.model.frame A \code{data.frame}. Each column corresponds to a phenotypic predictor or confounder
#' @param model.formula A \code{formula} object. Should have no RHS. Corresponds to the columns in \code{pheno.model.frame}. Exclusive with the (\code{predictor}, \code{confounder}) arguments.
#' @param ngenes An integer. Controls the number of genes to report for each test class.
#' @param p.adjust.method The name of p-value adjustment method to use. \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"} 
#' @param one.vs.all A flag indicating whether or not to use a one-vs-all design. This is most useful for initial exploration. 
#' @param summary.transform A function to be applied to \code{exprs} to transform it before calculating summary statistics
#' @param ... Additional parameters to pass to the strategy
#' @S3method do.diffanal2 default
#' @method do.diffanal2 default
#' @seealso \link{df2.t.test}, \link{df2.perm.t.test}, \link{df2.limma}
do.diffanal2.default <- function(exprs, pheno.model.frame, model.formula = NULL, predictor = NULL, confounders = list(), strategy = c('t.test', 'limma', 'perm.t'), ngenes = 250, p.adjust.method = 'fdr', one.vs.all = T, summary.transform = unlog, do.F.stat = F, ...){
  
  return(.do.diffanal2(exprs, pheno.model.frame, model.formula, predictor, confounders, strategy, ngenes, p.adjust.method, one.vs.all, summary.transform, do.F.stat = do.F.stat,
                       ...))

}

#' @rdname do.diffanal2
#' @title do.diffanal2.ExpressionSet
#' @description The input data is an ExpressionSet object from Bioconductor
#' @param eSet An ExpressionSet object
#' @method do.diffanal2 ExpressionSet
#' @S3method do.diffanal2 ExpressionSet
do.diffanal2.ExpressionSet <- function(eSet, model.formula = NULL, predictor = NULL, confounders = list(), strategy = c('t.test', 'limma', 'perm.t'), ngenes = 250, p.adjust.method = 'fdr', one.vs.all = T, summary.transform = unlog, do.F.stat = F, ...){
  dots <- list(...)
  exprs <- exprs(eSet)
  dots$pheno.model.frame <- pData(eSet)
  # Must unwrap ... to pick out pheno.model.frame passed as NULL, update it, and 
  # rewrap it. 
  args <- c(list(exprs=exprs, model.formula=model.formula, predictor=predictor, 
            confounders=confounders, strategy=strategy, ngenes=ngenes, 
            p.adjust.method=p.adjust.method, one.vs.all=one.vs.all, 
            summary.transform=summary.transform), do.F.stat = do.F.stat, dots)
  
  do.call(.do.diffanal2, args)
}


#' @title diffanal2 t-Test
#' @param x A numeric \code{matrix} of gene expression data or ExpressionSet object
#' @param pheno.model.frame A \code{data.frame}. A subset of columns corresponds to a phenotypic predictor or confounder. If x is an ExpressionSet object, this parameter will be drawn from it. 
#' @param model.formula A \code{formula} object. Should have no RHS. Corresponds to the columns in \code{pheno.model.frame}. Exclusive with the (\code{predictor}, \code{confounder}) arguments.
#' @param ngenes An integer. Controls the number of genes to report for each test class.
#' @param p.adjust.method The name of p-value adjustment method to use. \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"} 
#' @param one.vs.all A flag indicating whether or not to use a one-vs-all design. This is most useful for initial exploration. 
#' @param summary.transform A function to be applied to \code{exprs} to transform it before calculating summary statistics
#' @param do.F.stat Whether or not to compute an F statistic for each gene using \code{aov}
#' @param ... Additional parameters to pass to the strategy
#' @return A \code{diffanal.results} instance
#' @export
#' @seealso \link{do.diffanal2}, \link{df2.perm.t.test}, \link{df2.limma}
df2.t.test <- function(x, pheno.model.frame = NULL, model.formula = NULL, predictor = NULL, confounders = list(), ngenes = 250, p.adjust.method = 'fdr', one.vs.all = T, summary.transform = unlog, do.F.stat = F,...){
  return( do.diffanal2(x, pheno.model.frame = pheno.model.frame, model.formula=model.formula, 
                       predictor=predictor, confounders=confounders, ngenes=ngenes, 
                       p.adjust.method=p.adjust.method, one.vs.all=one.vs.all, 
                       summary.transform=summary.transform, strategy="t.test", 
                       do.F.stat=do.F.stat, ...) )
}

#' @title diffanal2 Permutation t-Test
#' @param x A numeric \code{matrix} of gene expression data or ExpressionSet object
#' @param pheno.model.frame A \code{data.frame}. A subset of columns corresponds to a phenotypic predictor or confounder. If x is an ExpressionSet object, this parameter will be drawn from it. 
#' @param model.formula A \code{formula} object. Should have no RHS. Corresponds to the columns in \code{pheno.model.frame}. Exclusive with the (\code{predictor}, \code{confounder}) arguments.
#' @param ngenes An integer. Controls the number of genes to report for each test class.
#' @param nperm An integer. Controls the number of permutations to generate. 
#' @param exhaustive A Boolean. If true, try to compute all permutations unless there are too many. 
#' @param ncores An integer. Controls the number of processors to use, using the \code{doMC} backend. Defaults to 1 which will instead use the \code{Sequential} backend
#' @param p.adjust.method The name of p-value adjustment method to use. \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"} 
#' @param one.vs.all A flag indicating whether or not to use a one-vs-all design. This is most useful for initial exploration. 
#' @param summary.transform A function to be applied to \code{exprs} to transform it before calculating summary statistics
#' @param do.F.stat Whether or not to compute an F statistic for each gene using \code{aov}
#' @param ... Additional parameters to pass to the strategy
#' @return A \code{diffanal.results} instance
#' @export
#' @seealso \link{do.diffanal2}, \link{df2.t.test}, \link{df2.limma}
df2.perm.t.test <- function(x, pheno.model.frame = NULL, model.formula = NULL, predictor = NULL, confounders = list(), ngenes = 250, nperm = 1, exhaustive = F, ncores = 1, p.adjust.method = 'fdr', one.vs.all = T, summary.transform = unlog, do.F.stat = F,...){
  do.diffanal2(x, pheno.model.frame = pheno.model.frame, model.formula=model.formula, 
               predictor=predictor, confounders=confounders, ngenes=ngenes, 
               p.adjust.method=p.adjust.method, one.vs.all=one.vs.all, 
               summary.transform=summary.transform, strategy="perm.t", 
               do.F.stat=do.F.stat, nperm=nperm, exhaustive=exhaustive, 
               ncores=ncores,
               ...)
}



#' @title diffanal2 limma Linear Model + Moderated t-Test
#' @description Uses limma's lmFit -> contrasts.fit -> eBayes workflow to calculate the effect contributed by the predictor variable independent of the confounder variable, and uses a moderated t-test to calculate significance. 
#' @param x A numeric \code{matrix} of gene expression data or ExpressionSet object
#' @param pheno.model.frame A \code{data.frame}. A subset of columns corresponds to a phenotypic predictor or confounder. If x is an ExpressionSet object, this parameter will be drawn from it. 
#' @param model.formula A \code{formula} object. Should have no RHS. Corresponds to the columns in \code{pheno.model.frame}. Exclusive with the (\code{predictor}, \code{confounder}) arguments.
#' @param ngenes An integer. Controls the number of genes to report for each test class.
#' @param p.adjust.method The name of p-value adjustment method to use. \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"} 
#' @param one.vs.all A flag indicating whether or not to use a one-vs-all design. This is most useful for initial exploration. 
#' @param summary.transform A function to be applied to \code{exprs} to transform it before calculating summary statistics
#' @param robust A Boolean. Passed along to eBayes to test for robustness to outliers.
#' @param ... Additional parameters to pass to the strategy
#' @return A \code{diffanal.results} instance
#' @export
#' @seealso \link{do.diffanal2}, \link{df2.t.test}, \link{df2.perm.t.test}
df2.limma <- function(x, pheno.model.frame = NULL, model.formula = NULL, predictor = NULL, confounders = list(), ngenes = 250, p.adjust.method = 'fdr', one.vs.all = T, summary.transform = unlog, robust = F, ...){
  return( do.diffanal2(x, pheno.model.frame = pheno.model.frame, model.formula=model.formula, 
                       predictor=predictor, confounders=confounders, ngenes=ngenes, 
                       p.adjust.method=p.adjust.method, one.vs.all=one.vs.all, 
                       summary.transform=summary.transform, strategy="limma", robust = robust,
                       ...) )
}
