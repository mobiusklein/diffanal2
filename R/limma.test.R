#' @import limma
#' @author Joshua Klein

library(limma)

#' @title Running limma through diffanl2
#' @param x A numeric \code{matrix} of gene expression data or ExpressionSet object
#' @param pheno.model.frame A \code{data.frame}. A subset of columns corresponds to a phenotypic predictor or confounder. If x is an ExpressionSet object, this parameter will be drawn from it. 
#' @param model.formula A \code{formula} object. Should have no RHS. Corresponds to the columns in \code{pheno.model.frame}. Exclusive with the (\code{predictor}, \code{confounder}) arguments.
#' @param ngenes.results An integer. Controls the number of genes to report for each test class. May be NULL
#' @param p.adjust.method The name of p-value adjustment method to use. \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"} 
#' @param one.vs.all A flag indicating whether or not to use a one-vs-all design. This is most useful for initial exploration. 
#' @param summary.transform A function to be applied to \code{exprs} to transform it before calculating summary statistics
#' @param robust A Boolean. Passed along to eBayes to test for robustness to outliers.
diffanal.limma.test <- function(exprs, pheno.model.frame, model.formula, 
                                p.adjust.method = 'fdr', one.vs.all = F, robust = F,
                                ngenes.results = 250, summary.transform = unlog, 
                                ...){

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
  
  ## Fit the contrasts model to calculate difference in expression between each
  ## group specified in contrast.pairs, and calculate log2 fold change for each 
  ## gene in each pair of conditions.
  fit.contrasts <- contrasts.fit(fit=restrict.fit, contrasts=contrasts.pairs)
  
  ## Carry out last model fitting step, scoring each contrast and evaluating the 
  ## moderatd t-test for each gene.
  fit.ebayes <- eBayes(fit.contrasts, robust=robust)
  
  ## Extract statistics that need transforming before being repoted, and 
  ## do said transformation
  adj.scores <- apply(fit.ebayes$p.value, 2, p.adjust, p.adjust.method)

  means <- summary.transform(restrict.fit$coefficients)
  colnames(means) <- paste(colnames(means), 'mean', sep=JOIN.CHR)
  
  std.devs <- summary.transform(restrict.fit$stdev.unscaled)
  colnames(std.devs) <- paste(colnames(std.devs), 'std.dev', sep=JOIN.CHR)
  
  fold.changes <- summary.transform(fit.ebayes$coefficients)
  colnames(fold.changes) <- paste(colnames(fold.changes), 'fold.change', sep=JOIN.CHR)
  
  ## Save the fits for attaching to the output later. Save them now so that
  ## they are copied and the names are preserved, the column names are later mutated
  ## and won't match limma's expectations.
  fits <- list(restricted.fit = restrict.fit, eBayes.fit = fit.ebayes, 
               contrasts = contrasts.pairs)
  
  if(one.vs.all){
    ## Fitting the complemented Model causes the values to be on a different scale from the
    ## uncomplemented model. To keep things n the same scale, average scores from all other
    ## conditions. These means are not used to compute fold change, and the fold change will
    ## be slightly different than what would result from using these means. 
    ## They are just presented for completeness.
    
    comp.means <- means / fold.changes
    colnames(comp.means) <- paste("Not.", colnames(means),sep='')
    means <- cbind(means, comp.means)

    ## This is not aggregated the same way the mean is, and so this value will
    ## not behave in the same way when separated from the contrasts
    comp.std.devs.cols <- llply(1:ncol(std.devs), function(i, std.devs){
      rowMeans(std.devs[,-i])
    }, std.devs)
    comp.std.devs <- do.call(cbind, comp.std.devs.cols)
    colnames(comp.std.devs) <- paste("Not.", colnames(std.devs), sep='')
    std.devs <- cbind(std.devs, comp.std.devs)
  }
                   
  
  # Fix column names and add descriptors to each "facet" of the MArrayLM object
  colnames(fit.ebayes$p.value) <- paste(fix.dash.names(condition.names, 
                                                       colnames(fit.ebayes$p.value)), 
                                        'p.value', sep=JOIN.CHR)
  colnames(adj.scores) <- gsub('p.value', "adj.p.value", colnames(fit.ebayes$p.value)) 
                                
  colnames(fit.ebayes$lods) <- paste(fix.dash.names(condition.names, 
                                                    colnames(fit.ebayes$lods)), 
                                     "lods", sep=JOIN.CHR)
  
  colnames(fit.ebayes$t) <- paste(fix.dash.names(condition.names, colnames(fit.ebayes$t)), 
                                  't.score', sep=JOIN.CHR)
  
  results <- do.summaries(summary.transform(exprs), partitions, pairs, one.vs.all, 
                          means = means, std.devs=std.devs, fold.changes = fold.changes)
  
  results <- c(list(gene = row.names(fit.ebayes), 
                  p.values=fit.ebayes$p.value, adj.p.values = adj.scores, 
                  lods=fit.ebayes$lods, t.scores=fit.ebayes$t),
                  results, 
                  F.stats = list(data.frame(F.score=fit.ebayes$F, 
                                           F.p.value=fit.ebayes$F.p.value)))
  
  results.table <- do.call(cbind, results)
  colnames(results.table) <- c("gene", c(unlist(lapply(results, colnames))))
    
  # Collate genes by adj.p.value (ascending) for each class
  results.table <- sortby(results.table, 
                                           colnames(results$adj.p.value), 
                                           ngenes.results)
  
  results.table <- results.table[!is.na(results.table[,1]),]
  
  results.table <- diffanal.results(results.table, ngenes=ngenes.results, 
                                          # gene is not a group
                                          column.groups=lapply(results[-1], colnames))
  
  attr(results.table, "fits") <- fits
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

#' When creating contrasts, "-" is used to denote difference. This isn't convenient when trying to access the resulting column name later. Easier to replace every instance of "-" not in a label name with JOIN.CHR 
#' @param condition.names The names of the data.frame or matrix columns
#' @param name.data The namse of the predictor variable levels to match
#' @return A character vector of column names
fix.dash.names <- function(condition.names, name.data){
  pattern <- paste("(", paste(condition.names, collapse="|"), 
                   ")-", sep='')
  cleaned.names <- gsub(pattern=pattern, replacement="\\1..", name.data)
  return(cleaned.names)
}


#' When creating a model matrix, the name of the predictor variable is prepended to its
#' levels' column names, making it harder to reference in context. Since we only care 
#' about the one predictor, we can remove the predictor variable's name from the column 
#' names, to be more consistent with the other strategies
#' @param predictor.name The name of the predictor variable column
#' @param name.data The column names to be cleaned
#' @return A character vector of column names
simplify.names <- function(predictor.name, name.data){
  cleaned.names <- gsub(pattern=predictor.name, replacement='', name.data)
  return(cleaned.names)
}


#' @title Creating contrasts for limma's \link{limma::contrasts.fit} function
#' http://permalink.gmane.org/gmane.science.biology.informatics.conductor/46235 describes the 
#' justification for why to use the arithmetic mean of the contrasted "Other" group.
#' @param predictor.col.names A vector of the levels of the predictor variable to contrast
#' @param one.vs.all A boolean flag to decide whether to create pairwise or one-vs-all contrasts
#' @return a matrix of contrasts
compute.contrast.pairs <- function(predictor.col.names, one.vs.all = F){
  # Special case for the one.vs.all case.
  if(one.vs.all){
    contrasts.strings <- sapply(seq_along(predictor.col.names), function(i){
      paste(predictor.col.names[i], ' - ((',
            do.call(paste, as.list(c(predictor.col.names[-i], sep='+'))), 
            ")/", length(predictor.col.names[-i]), ")"
            , sep='')
    })
    contrasts <- makeContrasts(contrasts=contrasts.strings, levels=predictor.col.names)
    colnames(contrasts) <- sapply(predictor.col.names, function(name){
      paste(name, paste("Not.", name, sep=''), sep=JOIN.CHR)
    })
    return(contrasts)
  }
  else{
    # Create order-independent unique pairs
    pairs <- t(combn(predictor.col.names, 2))
    # Construct "A-B" expressions for contrast
    contrasts.strings <- unlist(alply(pairs, 1, function(pair){
      pair <- lapply(pair, as.character)
      do.call(paste, c(pair, sep='-'))
    }))
    
    contrasts <- makeContrasts(contrasts=contrasts.strings, levels=predictor.col.names)
    return(contrasts)
  }
}