#' @import genefilter

library(genefilter)

#' Runs a t.test for each gene in each partition pairing
#' @param exprs A numeric \code{matrix} of gene expression data
#' @param pheno.model.frame A \code{data.frame}. Each column corresponds to a phenotypic predictor or confounder
#' @param model.formula A \code{formula} object. Should have no RHS. Corresponds to the columns in \code{pheno.model.frame}. Exclusive with the (\code{predictor}, \code{confounder}) arguments.
#' @param ngenes.results An integer. Controls the number of genes to report for each test class. May be NULL
#' @param p.adjust.method The name of p-value adjustment method to use. \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"} 
#' @param one.vs.all A flag indicating whether or not to use a one-vs-all design. This is most useful for initial exploration. 
#' @param summary.transform A function to be applied to \code{exprs} to transform it before calculating summary statistics
#' @return A \code{diffanal.results} instance
#' @seealso df2.t.test
diffanal.t.test <- function(exprs, pheno.model.frame, model.formula, 
                            p.adjust.method = 'fdr', one.vs.all = F, 
                            ngenes.results = 250, summary.transform = unlog,
                            do.F.stat = F, ...){
  
  # Cartesian Product of cofactors separates the specific sample names
  partitions <- compute.partitions(pheno.model.frame, model.formula)
  pairs <- compute.partition.pairs(partitions, one.vs.all = one.vs.all)
  
  # Compute t.test p.values
  results <- do.t.test(exprs, partitions, pairs, one.vs.all = one.vs.all)
  p.values <- results$p.values
  t.scores <- results$t.scores
  # Compute F statistics and F.p.values  
  F.stats = NULL
  if(do.F.stat){
    F.stats <- do.aov(exprs, pheno.model.frame, model.formula)  
  }
  
  adj.p.values <- apply(p.values, 2, p.adjust, p.adjust.method)
  colnames(adj.p.values) <- gsub('p.value', "adj.p.value", colnames(p.values))
  
  
  # Summary statistics  
  transform.exprs <- summary.transform(exprs)
  
  results <- do.summaries(transform.exprs, partitions, pairs, one.vs.all)
  
  # --------------------------------------------------------------------------
  # Construct results table from the analysis. Specialized for t-test currently.
  # Break up into helper functions to generalize. 
  # --------------------------------------------------------------------------
  
  results <- c(list(p.values=p.values, adj.p.values = adj.p.values, 
                    t.scores = t.scores), results)
  if(!is.null(F.stats)){
    results <- c(results, list(F.stats=F.stats))
  }
  
  results.table <- do.call(cbind, results)
  
  # Include probe names in a separate column as row.names are not preserved by all operations
  results.table <- cbind(row.names(results.table), results.table)
  colnames(results.table) <- c("gene", unlist(lapply(results, colnames)))
  
  results.table <- sortby(table = results.table, 
                                         column.grp = colnames(results$adj.p.value),
                                         ngenes = ngenes.results
  )
  
  results.table <- diffanal.results(results.table, column.groups=lapply(results, colnames))
  return(results.table)
}

do.t.test <- function(exprs, partitions, pairs, one.vs.all = F, do.test = T){
  results <- apply(pairs, 1, function(pairs, exprs, partitions){    
    exprs.subset <- NULL;
    if(one.vs.all){
      exprs.subset <- exprs
    } else {
      exprs.subset <- exprs[,c(partitions[[ pairs[[1]] ]], partitions[[pairs[2]]])];
    }
    grouping <- as.factor(colnames(exprs.subset) %in% partitions[[ pairs[[1]] ]])
    test.stats <- t.score(exprs.subset, grouping, do.test=do.test)
    return(test.stats)
            
  },exprs,partitions)
  output <- list()
  if(do.test){
    p.vals <- t(ldply(results, function(x)x$p.value))
    row.names(p.vals) <- row.names(exprs)
    colnames(p.vals) <- apply(pairs, 1, function(ps){paste(ps[[1]], ps[[2]], "p.value", sep=JOIN.CHR)})
    output$p.values <- p.vals
  }
  score <- t(ldply(results, function(x)x$score))
  row.names(score) <- row.names(exprs)
  colnames(score) <- apply(pairs, 1, function(ps){paste(ps[[1]], ps[[2]], "t.score", sep=JOIN.CHR)})
  output$t.scores <- score
  return(output)
}

