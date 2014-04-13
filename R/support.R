
# Converts a list of terms into a formula object for building models with
formula.from.list <- function(predictor, confounders = list()){
  # Attach the predictor and tilde first
  relation <- paste("~", predictor)
  # Combine any confounders into a chain of additions 
  trailing.terms <- do.call(paste, c(confounders, sep = " + "))
  # If there are trailing terms to merge
  if(length(trailing.terms) > 0){
    # Attach the predictor and the confounder string
    relation <- paste(relation, trailing.terms, sep=" + ")
  }
  full.formula <- as.formula(relation)
  return(full.formula)
}

unlog <- function(exprs){
  unlogged.exprs = 2 ** exprs
  return(unlogged.exprs)
}

# Remove samples that are NA in any of the model dimensions
clean.samples <- function(exprs, pheno, model.formula.cols){
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

partition.pairs <- function(partitions, one.vs.all=F){
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

do.mean <- function(exprs, partitions, one.vs.all){
  results <- lapply(partitions, function(partition){
    exprs.subset <- exprs[,partition];
    result <- rowMeans(exprs.subset)
    return(result)
  })
  results <- as.data.frame(do.call(cbind, results))
  names(results) <- paste(names(partitions), 'mean', sep=JOIN.CHR)

  # one.vs.all complementary subsets must be separately specified
  if(one.vs.all){
    sample.names <- colnames(exprs)
    others <- lapply(partitions, function(partition){
      other <- !(sample.names %in% partition)
      other <- sample.names[other]
      exprs.subset <- exprs[,other];
      result <- rowMeans(exprs.subset)
      return(result)
    })
    names(others) <- paste("Not", names(results), sep='.')
    results <- cbind(results, others)
  }
  return(results)
}

do.std.dev <- function(exprs, partitions){
  results  <- lapply(partitions, function(partition){
    exprs.subset <- exprs[, partition];
    result <- rowSds(exprs.subset)
    return(result)
  })
  results <- as.data.frame(do.call(cbind, results))
  names(results) <- paste(names(partitions), 'std.dev', sep=JOIN.CHR)
  return(results)
}

do.fold.change <- function(means, pairs, transform.fn, one.vs.all){
  # Adds the suffix to each group name to match the equivalent mean
  # column in the means frame
  pair.cols <- apply(pairs, 2, paste, 'mean', sep=JOIN.CHR)
  results <- as.data.frame(apply(pair.cols, 1, function(pair){    
    result <- transform.fn(means[,pair[[1]]]) / transform.fn(means[,pair[[2]]])
    return(result)
  }))
  
  names(results) <- apply(pairs, 1, function(ps){paste(ps[[1]], ps[[2]], 
                                                       'fold.change', sep=JOIN.CHR)})
  rownames(results) <- rownames(means)
  return(results)
}

# Get only significant columns of a frame of p.values
which.significant <- function(scores, alpha = 0.05){
  results <- apply(scores, 1, function(row){
    return(any(row < alpha))
  })
  return(results)
}

