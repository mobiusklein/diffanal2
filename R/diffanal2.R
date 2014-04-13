require(plyr)
require(genefilter)
require(Biobase)
require(limma)

# Column names with aggregate meaning will have tokens joined by the JOIN.CHR
# e.g. cls1-cls2-p.value means the p.value comparing cls1 and cls2.
JOIN.CHR = "-"
SAMPLE.ID.COL <- '.rownames'
STRATEGY = list(t.test = function(){}, survival = function(){}, lm = function(){})

#' 
#' 
#' 
#' 
#' @param exprs numeric matrix of gene expression data
#' @param pheno.frame data.frame. Each column corresponds to a phenotypic predictor or confounder
#' @param model.formula should have no RHS. Corresponds to the columns in the \code{pheno} data.frame
#' @param ngenes controls the number of genes to report
#' @param alpha significance threshold value
#'   
diffanal2 <- function(exprs, pheno.frame, model.formula, ngenes = 250, alpha = 0.05, p.adjust.method = 'none', one.vs.all = F){
  # Gets parameter list, for later writing to file
  print(match.call())
  
  model.formula.terms <- terms(model.formula)
  # Specify without response variable
  
  model.formula.cols <- row.names(attr(model.formula.terms, 'factors'))
  
  cleaned <-clean.phenotypes(exprs, pheno.frame, model.formula.cols)
  exprs <- cleaned[["exprs"]]
  pheno.frame <- cleaned[["pheno"]]
  
  model.cols <- prepare.pheno.model.frame(pheno.frame, model.formula.cols)
  
  # Cartesian Product of cofactors separates the specific sample names
  partitions <- compute.partitions(model.cols, model.formula)
  
  pairs <- partition.pairs(partitions, one.vs.all = one.vs.all)
  
  scores <- do.t.test(exprs, partitions, pairs, one.vs.all = one.vs.all)
  
  # Adjust p.values and extract signficant rows
  scores <- apply(scores, 2, p.adjust, p.adjust.method)
  signif.comparisons <- which.significant(scores, alpha)
  
  signif.scores <- scores[signif.comparisons, ]
  signif.exprs <- exprs[signif.comparisons, ]
  means  <- do.means(signif.exprs, partitions)
  
  
  
  return(list(scores=scores[signif.comparisons,], means=means[signif.comparisons,], partitions=partitions, exprs=exprs, pheno=model.cols))
}

# Remove samples that are NA in any of the model dimensions
clean.phenotypes <- function(exprs, pheno, model.formula.cols){
  # Must be as.data.frame in order to prevent vectorizing when only one 
  # column is selected in formula. 
  na.samples <- Reduce(`|`, as.data.frame(is.na(pheno[,model.formula.cols])))
  # Drop pheno rows with NA attributes
  pheno <- pheno[!na.samples,]
  # Drop corresponding exprs columns 
  exprs <- exprs[,!na.samples]
  return(list(exprs=exprs, pheno=pheno))
}

prepare.pheno.model.frame <- function(pheno, model.formula.cols){
  model.cols <- as.data.frame(pheno[, model.formula.cols])
  names(model.cols) <- model.formula.cols
  model.cols <- name_rows(model.cols)
  model.cols <- model.cols[order(model.cols[model.formula.cols[1]]),]
  # column vectors with the name attribute breaks dlply  
  for(i in seq_along(model.cols)){
    names(model.cols[i]) <- NULL
    attr(model.cols[i], '.Names') <- NULL
  }
  return(model.cols)
}

compute.partitions <- function(model.frame, model.formula){
  model.formula.terms <- terms(model.formula)
  # Specify without response variable
  model.formula.cols <- row.names(attr(model.formula.terms, 'factors'))
  # Cartesian Product of cofactors separates the specific sample names
  partitions <- dlply(model.frame, model.formula.cols, function(set) set[,'.rownames'])
  return(partitions)
}

partition.pairs <- function(partitions, one.vs.all=T){
  # If using a one.vs.all design, create a simple partition pair frame
  # with 
  if(one.vs.all){
    pairs <- as.data.frame(cbind(names(partitions), paste("Not", names(partitions), sep='.')))
    names(pairs) <- c('group1', 'group2')
    return(pairs)
  }
  pairs <- expand.grid(group1=names(partitions), group2=names(partitions)) 
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
  return(pairs)
}

do.t.test <- function(exprs, partitions, pairs, one.vs.all = F, statistic = F){
  results <- apply(pairs, 1, function(pairs, exprs, partitions){
    exprs.subset <- NULL;
    if(one.vs.all){
      exprs.subset <- exprs
    } else {
      exprs.subset <- exprs[,c(partitions[[pairs[1]]], partitions[[pairs[2]]])];
    }
    
    test.stats <- rowttests(exprs.subset, as.factor(colnames(exprs.subset) %in% partitions[[pairs[1]]]), tstatOnly=statistic)
    if(statistic){
      return(test.stats["statistic"])
    } else{
      return(test.stats["p.value"]) 
    }    
  },exprs,partitions)
  
  results <- do.call(cbind, results)
  result.type <- ifelse(statistic, 't.score', 'p.value')
  names(results) <- names(results) <- apply(pairs, 1, function(ps){paste(ps[[1]], ps[[2]], result.type, sep=JOIN.CHR)})
  
  return(results)
}

do.means <- function(exprs, partitions){
  results <- lapply(partitions, function(partition){
    exprs.subset <- exprs[,partition];
    result <- rowMeans(exprs.subset)
    return(result)
  })
  results <- as.data.frame(do.call(cbind, results))
  names(results) <- paste(names(partitions), 'mean', sep=JOIN.CHR)
  return(results)
}

# Get only significant columns of a frame of p.values
which.significant <- function(scores, alpha = 0.05){
  pairs <- names(scores)
  results <- apply(scores, 1, function(row){
    return(any(row < alpha))
  })
  
  return(results)
}

do.limma.lmFit <- function(exprs, pheno, model.formula, ngenes=250, alpha){
  design <- model.matrix(model.formula, pheno)
  fit <- lmFit(exprs, design)
  fit.ebayes <- eBayes(fit)
  #fit.top <- topTable(fit.ebayes, number=ngenes)
  return(fit.ebayes)
}