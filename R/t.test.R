#' @import genefilter

library(genefilter)

#' diffanal.t.test 
#' Runs a t.test for each gene in each partition pairing
diffanal.t.test <- function(exprs, pheno.model.frame, model.formula, 
                            alpha = 0.05, p.adjust.method = 'fdr', one.vs.all = F, 
                            ngenes = 250, 
                            summary.transform = unlog){
  
  # Cartesian Product of cofactors separates the specific sample names
  partitions <- compute.partitions(pheno.model.frame, model.formula)
  pairs <- compute.partition.pairs(partitions, one.vs.all = one.vs.all)
  
  # Compute t.test p.values
  scores <- do.t.test(exprs, partitions, pairs, one.vs.all = one.vs.all)
  # Compute F statistics and F.p.values
  
  
  signif.comparisons <- which.significant(scores, alpha)
  if(length(signif.comparisons[signif.comparisons]) == 0){
    stop("No genes were significant")
  }
  signif.scores <- scores[signif.comparisons, ]
  F.stats <- do.aov(exprs[signif.comparisons, ], pheno.model.frame, model.formula)
  adj.scores <- apply(signif.scores, 2, p.adjust, p.adjust.method)
  colnames(adj.scores) <- gsub('p.value', "adjusted.p.value", colnames(signif.scores))
  
  
  # Summary statistics
  signif.exprs <- exprs[signif.comparisons, ]
  
  transform.exprs <- summary.transform(signif.exprs)
  
  results <- do.summaries(transform.exprs, partitions, pairs, one.vs.all)
  
  
#   means  <- do.mean(transform.exprs, partitions, one.vs.all=one.vs.all)
#   std.devs <- do.std.dev(transform.exprs, partitions)
#   fold.changes <- do.fold.change(means, pairs, one.vs.all=one.vs.all)
  
  # --------------------------------------------------------------------------
  # Construct results table from the analysis. Specialized for t-test currently.
  # Break up into helper functions to generalize. 
  # --------------------------------------------------------------------------
  
  results <- c(list(p.value=signif.scores, adj.p.value = adj.scores), results, list(F.stats=F.stats))
  
  results.table <- do.call(cbind, results)
  
  # Include probe names in a separate column as row.names are not preserved by all operations
  results.table <- cbind(row.names(results.table), results.table)
  #
  colnames(results.table) <- c("gene", unlist(lapply(results, colnames)))
  
  # Collate genes by p.value (ascending) for each class
  sorts <- lapply(names(results$p.value), function(column){
    order(results.table[,column])
  })
  # Merge the assembled gene list, retaining @ngenes from each set. At most the set will
  # contain number of comparisons
  top.genes <- Reduce(function(a, b){
    a <- unlist(a)
    b <- unlist(b)
    return(union(a, b[1:ngenes]))
  }, sorts[-1], sorts[[1]][1:ngenes])
  
  results.table <- results.table[unlist(top.genes),]
  results.table <- diffanal.results.table(results.table, ngenes=ngenes, 
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

