#' @import limma
#' @author Joshua Klein

library(limma)

#' diffanal.limma.test
diffanal.limma.test <- function(exprs, pheno.model.frame, model.formula, 
         p.adjust.method = 'fdr', one.vs.all = F, robust = F,
         ngenes = 250, summary.transform = unlog, ...){

  # Ensure there is no intercept term in the design matrix
  model.formula <- adjust.formula(model.formula)
  design <- model.matrix(model.formula, pheno.model.frame)

  # Get the columns of the design matrix that correspond to the predictor variable.
  # These variables correspond to the coefficients of the fit
  predictor = get.predictor(model.formula)
  predictor.cols <- attr(design, 'assign') == 1
  # Remove the name of the predictor variable from the column names
  colnames(design)[predictor.cols] <- simplify.names(predictor,
                                                     colnames(design)[predictor.cols])
    
  predictor.col.names <- colnames(design)[predictor.cols]
    
  # Only partition on predictor
  partitions <- compute.partitions(pheno.model.frame, get.predictor(model.formula))
  pairs <- compute.partition.pairs(partitions, one.vs.all=one.vs.all)
  
  condition.names <- as.character(unique(unlist(pairs)))
  
  contrasts.pairs <- compute.contrast.pairs(predictor.col.names, one.vs.all)
  
  # --------------------------------------------------------------------------
  # Must call limma workflow in this order. Otherwise results are meaningless.
  # No error is raised if tests are called out of order. 
  # --------------------------------------------------------------------------
  
  lm.method <- ifelse(robust, "robust", "ls")
  fit.lm <- lmFit(exprs, design, method=lm.method)
  
  # Drop any confounder or interaction terms that absorb significance from the 
  # predictor.They are not given in the contrast expression, and contrasts.fit 
  # won't know what to do with them. http://stats.stackexchange.com/a/64250
  restrict.fit <- fit.lm[, predictor.col.names]
  
  # Subsetting MArrayLM has odd behavior, must copy the old design matrix and subset
  # separately.
  restrict.fit$design <- fit.lm$design[,predictor.col.names]
  fit.contrasts <- contrasts.fit(fit=restrict.fit, contrasts=contrasts.pairs)
  
  # Carry out last model fitting step, scoring each contrast and evaluating the 
  # moderatd t-test for each gene.
  fit.ebayes <- eBayes(fit.contrasts, robust=robust)
  
  # Determine which comparisons are significant.
  # Separately compute the adjusted score here because decideTests does not return
  # the adjusted p.value computed.
  adj.scores <- apply(fit.ebayes$p.value, 2, p.adjust, p.adjust.method)

  means <- summary.transform(restrict.fit$coefficients)
  colnames(means) <- paste(colnames(means), 'mean', sep=JOIN.CHR)
  
  std.devs <- summary.transform(restrict.fit$stdev.unscaled)
  colnames(std.devs) <- paste(colnames(std.devs), 'std.dev', sep=JOIN.CHR)
  
  if(one.vs.all){
    
    # Construct the complementary design matrix with the "not" classes
    complement.design <- make.complement.design(design,                                                                   predictor.col.names)
    # Fit the "Not" model
    complement.fit <- lmFit(exprs, complement.design)
    
    # Extract summary statistics
    comp.means <- summary.transform(complement.fit$coefficients)  
    colnames(comp.means) <- paste(colnames(comp.means), 'mean', sep=JOIN.CHR)
    means <- cbind(means, comp.means)
    
    comp.std.devs <- summary.transform(complement.fit$stdev.unscaled)
    colnames(comp.std.devs) <- paste(colnames(comp.std.devs), 'std.dev', sep=JOIN.CHR)
    std.devs <- cbind(std.devs, comp.std.devs)
  }
                        
  colnames(fit.ebayes$p.value) <- paste(fix.dash.names(condition.names, 
                                                       colnames(fit.ebayes$p.value)), 
                                        'p.value', sep=JOIN.CHR)
  colnames(adj.scores) <- paste(fix.dash.names(condition.names, 
                                               colnames(fit.ebayes$p.value)), 
                                'adj.p.value', sep=JOIN.CHR)
  colnames(fit.ebayes$lods) <- paste(fix.dash.names(condition.names, 
                                                    colnames(fit.ebayes$lods)), 
                                     "lods", sep=JOIN.CHR)
  
  colnames(fit.ebayes$t) <- paste(fix.dash.names(condition.names, colnames(fit.ebayes$t)), 
                                  't.score', sep=JOIN.CHR)
  
  results <- do.summaries(summary.transform(exprs), partitions, pairs, one.vs.all, 
                          means = means, std.devs=std.devs)
  
  results <- c(list(gene = row.names(fit.ebayes), 
                  p.value=fit.ebayes$p.value, adj.p.value = adj.scores, 
                  lods=fit.ebayes$lods, t.score=fit.ebayes$t),
                  results, 
                  F.stat = list(data.frame(F.score=fit.ebayes$F, 
                                           F.p.value=fit.ebayes$F.p.value)))
  
  results.table <- do.call(cbind, results)
  colnames(results.table) <- c("gene", c(unlist(lapply(results, colnames))))
    
  # Collate genes by adj.p.value (ascending) for each class
  results.table <- sort.diffanal.results(results.table, colnames(results$adj.p.value), ngenes)
  results.table <- results.table[!is.na(results.table[,1]),]
  
  results.table <- diffanal.results.table(results.table, ngenes=ngenes, 
                                          unique=F, 
                                          # gene is not a group
                                          column.groups=lapply(results[-1], colnames))
  return(results.table)
  
}

make.complement.design <- function(design, predictor.col.names){
  pred.cols <- as.matrix(design[,predictor.col.names])
  complement.cols <- apply(pred.cols, 2, function(col){
    return(as.numeric(!(col)))
  })    
  design[,predictor.col.names] <- complement.cols
  colnames(design)[seq_along(predictor.col.names)] <- paste("Not", colnames(complement.cols), sep='.')
  return(design)
}

fix.dash.names <- function(condition.names, name.data){
  pattern <- paste("(", paste(condition.names, collapse="|"), 
                   ")-", sep='')
  cleaned.names <- gsub(pattern=pattern, replacement="\\1..", name.data)
  return(cleaned.names)
}

simplify.names <- function(predictor.name, name.data){
  cleaned.names <- gsub(pattern=predictor.name, replacement='', name.data)
  return(cleaned.names)
}



compute.contrast.pairs <- function(predictor.col.names, one.vs.all = F){
  # Special case for the one.vs.all case.
  if(one.vs.all){
    contrasts.strings <- sapply(seq_along(predictor.col.names), function(i){
      paste(predictor.col.names[i], 
            do.call(paste, as.list(c(predictor.col.names[-i], sep='-'))), sep='-')
    })
    contrasts <- makeContrasts(contrasts=contrasts.strings, levels=predictor.col.names)
    colnames(contrasts) <- sapply(predictor.col.names, function(name){
      paste(name, paste("Not.", name, sep=''), sep=JOIN.CHR)
    })
    return(contrasts)
  }
  pairs <- expand.grid(group1=predictor.col.names, group2=predictor.col.names) 
  pairs <- pairs[!pairs[,1]==pairs[,2],]
  # Iteratively drop redundant pairs, A + B vs B + A
  keep.inds <- 1
  for(i in 2:dim(pairs)[1]){
    cur.pair <- pairs[i,]  
    rev.match <- apply(pairs[keep.inds,],1,function(pr){
      pr[[1]] == cur.pair[[2]] && pr[[2]] == cur.pair[[1]]
    })
    if(!any(rev.match)){
      keep.inds = c(keep.inds, i)
    }    
  }
  pairs <- pairs[keep.inds,]
  contrasts.strings <- unlist(alply(pairs, 1, function(pair){
    pair <- lapply(pair, as.character)
    do.call(paste, c(pair, sep='-'))
  }))
  
  contrasts <- makeContrasts(contrasts=contrasts.strings, levels=predictor.col.names)
  return(contrasts)
}