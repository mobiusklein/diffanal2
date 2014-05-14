#' @import genefilter

library(genefilter)

#' diffanal.t.test 
#' Runs a t.test for each gene in each partition pairing
diffanal.t.test <- function(exprs, pheno.model.frame, model.formula, 
                            p.adjust.method = 'fdr', one.vs.all = F, 
                            ngenes = 250, summary.transform = unlog,
                            do.F.stat = F, ...){
  
  # Cartesian Product of cofactors separates the specific sample names
  partitions <- compute.partitions(pheno.model.frame, model.formula)
  pairs <- compute.partition.pairs(partitions, one.vs.all = one.vs.all)
  
  # Compute t.test p.values
  scores <- do.t.test(exprs, partitions, pairs, one.vs.all = one.vs.all)
  # Compute F statistics and F.p.values  
  F.stats = NULL
  if(do.F.stat){
    F.stats <- do.aov(exprs, pheno.model.frame, model.formula)  
  }
  
  adj.scores <- apply(scores, 2, p.adjust, p.adjust.method)
  colnames(adj.scores) <- gsub('p.value', "adjusted.p.value", colnames(scores))
  
  
  # Summary statistics  
  transform.exprs <- summary.transform(exprs)
  
  results <- do.summaries(transform.exprs, partitions, pairs, one.vs.all)
  
  # --------------------------------------------------------------------------
  # Construct results table from the analysis. Specialized for t-test currently.
  # Break up into helper functions to generalize. 
  # --------------------------------------------------------------------------
  
  results <- c(list(p.value=scores, adj.p.value = adj.scores), results)
  if(!is.null(F.stats)){
    results <- c(results, list(F.stats=F.stats))
  }
  
  results.table <- do.call(cbind, results)
  
  # Include probe names in a separate column as row.names are not preserved by all operations
  results.table <- cbind(row.names(results.table), results.table)
  colnames(results.table) <- c("gene", unlist(lapply(results, colnames)))
  
  results.table <- sort.diffanal.results(results.table, colnames(results$adj.p.value))
  results.table <- diffanal.results(results.table, ngenes=ngenes, 
                                          unique=unique, column.groups=lapply(results, colnames))
  return(results.table)
}

do.t.test <- function(exprs, partitions, pairs, one.vs.all = F, statistic = F){
  results <- apply(pairs, 1, function(pairs, exprs, partitions){
    exprs.subset <- NULL;
    if(one.vs.all){
      exprs.subset <- exprs
    } else {
      exprs.subset <- exprs[,c(partitions[[pairs[1]]], partitions[[pairs[2]]])];
    }
    
    test.stats <- rowttests(exprs.subset, 
                            as.factor(colnames(exprs.subset) %in% partitions[[pairs[1]]]), 
                            tstatOnly=statistic)
    ifelse(statistic, return(test.stats["statistic"]), 
           return(test.stats["p.value"])) 
  },exprs,partitions)
  
  results <- do.call(cbind, results)
  result.type <- ifelse(statistic, 't.score', 'p.value')
  names(results) <- names(results) <- apply(pairs, 1, function(ps){paste(ps[[1]], ps[[2]], result.type, sep=JOIN.CHR)})
  
  return(results)
}

