require(genefilter)

diffanal.t.test <- function(exprs, pheno.model.frame, model.formula, 
                            alpha = 0.05, p.adjust.method = 'fdr', one.vs.all = F, 
                            fold.change.transform.fn = unlog){
  
  # Cartesian Product of cofactors separates the specific sample names
  partitions <- compute.partitions(pheno.model.frame, model.formula)
  pairs <- partition.pairs(partitions, one.vs.all = one.vs.all)
  
  # Compute t.test p.values
  scores <- do.t.test(exprs, partitions, pairs, one.vs.all = one.vs.all)
  signif.comparisons <- which.significant(scores, alpha)
  signif.scores <- scores[signif.comparisons, ]
  adj.scores <- apply(signif.scores, 2, p.adjust, p.adjust.method)
  colnames(adj.scores) <- gsub('p.value', p.adjust.method, colnames(signif.scores))
  
  # Summary statistics
  signif.exprs <- exprs[signif.comparisons, ]
  means  <- do.mean(signif.exprs, partitions, one.vs.all=one.vs.all)
  std.devs <- do.std.dev(signif.exprs, partitions)
  fold.changes <- do.fold.change(means, pairs, fold.change.transform.fn, one.vs.all=one.vs.all)
  
  results <- list(p.value=signif.scores, adj.p.value = adj.scores, 
                  fold.changes=fold.changes, means=means, std.devs=std.devs)
  
  return(results)
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

