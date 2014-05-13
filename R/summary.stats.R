#' Standard function for computing mean expression levels from the expression matrix
#' 
#' @param exprs A numeric matrix of expression data
#' @param partitions A list of labeled sets of column names in \code{exprs}. These lists are created using \code{compute.partitions}
#' @param one.vs.all Indicates that the test is being ran in a one-vs-all style, so pooled not-in-class means need to be computed as well. 
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
#' Compute the standard deviation of expression for each gene in each partition
#' 
#' @param exprs A numeric matrix of expression data
#' @param partitions A list of labeled sets of column names in \code{exprs}. These lists are created using \code{compute.partitions}
do.std.dev <- function(exprs, partitions, one.vs.all){
  results  <- lapply(partitions, function(partition){
    exprs.subset <- exprs[, partition];
    result <- rowSds(exprs.subset)
    return(result)
  })
  results <- as.data.frame(do.call(cbind, results))
  names(results) <- paste(names(partitions), 'std.dev', sep=JOIN.CHR)
  # one.vs.all complementary subsets must be separately specified
  if(one.vs.all){
    sample.names <- colnames(exprs)
    others <- lapply(partitions, function(partition){
      other <- !(sample.names %in% partition)
      other <- sample.names[other]
      exprs.subset <- exprs[,other];
      result <- rowSds(exprs.subset)
      return(result)
    })
    names(others) <- paste("Not", names(results), sep='.')
    results <- cbind(results, others)
  }
  return(results)
}

#' Compute the fold-change (mean ratio) between each pair of conditions.
#' 
#' @param means The output data.frame from \code{do.means}
#' @param pairs A data.frame of non-redundant partition name pairs to compare
do.fold.change <- function(means, pairs){
  # Adds the suffix to each group name to match the equivalent mean
  # column in the means frame
  pair.cols <- apply(pairs, 2, paste, 'mean', sep=JOIN.CHR)
  results <- as.data.frame(apply(pair.cols, 1, function(pair){
    result <- means[,pair[[1]]] / means[,pair[[2]]]
    return(result)
  }))
  col.names <- (apply(pairs, 1, function(ps){paste(ps[[1]], ps[[2]], 
                                                   'fold.change', sep=JOIN.CHR)}))
  names(results) <- col.names
  rownames(results) <- rownames(means)
  return(results)
}

#' Compute the median for each gene in each partition
#' @param exprs The gene expression data matrix
#' @param partitions The set of self-contained conditions to sample id mappings
#' @param one.vs.all Determine whether or not to create complementary Not-group data
#' @export
do.median <- function(exprs, partitions, one.vs.all){
  results <- lapply(partitions, function(partition){
    exprs.subset <- exprs[,partition];
    result <- rowMedians(exprs.subset)
    return(result)
  })
  results <- as.data.frame(do.call(cbind, results))
  names(results) <- paste(names(partitions), 'median', sep=JOIN.CHR)
  
  # one.vs.all complementary subsets must be separately specified
  if(one.vs.all){
    sample.names <- colnames(exprs)
    others <- lapply(partitions, function(partition){
      other <- !(sample.names %in% partition)
      other <- sample.names[other]
      exprs.subset <- exprs[,other];
      result <- rowMedians(exprs.subset)
      return(result)
    })
    names(others) <- paste("Not", names(results), sep=JOIN.CHR)
    results <- cbind(results, others)
  }
  return(results)
}

#' Compute the Median Absolute Deviation (MAD) for each row in a matrix or data.frame. 
#' @param exprs The matrix or data.frame to operate on
#' @param centers Each row's median. If NULL, will compute with \code{rowMedians(exprs)}
#' @param constant As parameter in \code{mad}
#' @param na.rm As parameter in \code{mad}
#' @param low As parameter in \code{mad}
#' @param high As parameter in \code{mad}
rowMads <- function(exprs, centers = NULL, constant = 1.4826, na.rm = F, low = F, high = F){
  if(is.null(centers)){
    centers = rowMedians(exprs)
    verbose("Center missing, recomputing")
  }
  mads <- sapply(seq_along(exprs[,1]), function(i){
    data <- exprs[i,]
    center <- centers[i]
    return( mad(data, center, constant, na.rm, low, high) )
  })
  names(mads) <- row.names(exprs)
  return( mads )
}

#' Compute the mad for each gene in each partition.
#' @param exprs A numeric matrix of gene expression data
#' @param medians The median for each row in exprs in each partition
#' @param partitions The set of self-contained conditions to sample id mappings
#' @param one.vs.all Boolean whether to compute complementary group statistics
#' @param constant see \code{stats::mad}
#' @param low see \code{stats::mad} 
#' @param high see \code{stats::mad}
#' @seealso \code{stats::mad}
do.mad <- function(exprs, medians, partitions, one.vs.all, constant=1.4826, low=F, high=F){
  results <- lapply(names(partitions), function(part.name){
    exprs.subset <- exprs[,partitions[[part.name]]];
    centers <- medians[,paste(part.name, 'median', sep=JOIN.CHR)]
    mads <- rowMads(exprs.subset, centers, constant, F, low, high)
    return(mads)
  })
  results <- as.data.frame(do.call(cbind, results))
  names(results) <- paste(names(partitions), 'mad', sep=JOIN.CHR)
  if(one.vs.all){
    sample.names <- colnames(exprs)
    others <- lapply(names(partitions), function(part.name){
      other <- !(sample.names %in% partitions[[part.name]])
      other <- sample.names[other]
      part.name <- paste("Not", part.name, sep=JOIN.CHR)
      exprs.subset <- exprs[,other];
      centers <- medians[,paste(part.name, 'median', sep=JOIN.CHR)]
      mads <- rowMads(exprs.subset, centers, constant, F, low, high)
      return(mads)
    })
    names(others) <- paste("Not", names(results), sep=JOIN.CHR)
    results <- cbind(results, others)
  }
  return(results)
}

# Computes the mean and standard deviation of the predictor's influence
# on each gene. Computes a separate value for each level of the predictor.
do.lm <- function(exprs, pheno.frame, model.formula, one.vs.all = F){
  data <- cbind(exprs[1,],pheno.frame)
  names(data) <- c('..data..', names(pheno.frame))
  data.formula <- as.formula(paste(c('..data..', as.character(model.formula)), collapse=' '))
  coefs <- coef(summary(lm(data.formula, data=data)))
  coefs <- t(coefs)
  colnames(coefs) <- simplify.names(get.predictor(model.formula), colnames(coefs))
  means <- coefs[1,]
  std.devs <- coefs[2,]
  return(list(means=means, std.devs=std.devs))
}

#' Computes F-score and significance of F-score
#' @param exprs A numeric matrix of gene expression data
#' @param pheno.frame A data.frame containing predictor and confounder labels for each sample
#' @param model.formula A formula object defining the predictor and confounder labels
do.aov <- function(exprs, pheno.frame, model.formula){
  verbose("Computing F statistic for entire model")
  model.formula <- adjust.formula(model.formula)
  time <- system.time(results <- alply(exprs, 1, function(row){
    data <- cbind(row,pheno.frame)
    names(data) <- c('..data..', names(pheno.frame))
    data.formula <- as.formula(paste(c('..data..', as.character(model.formula)), 
                                     collapse=' '))
    res <- unlist(summary(aov(data.formula, data=data)))
    F.score <- res[7]
    F.p.value <- res[9]
    return(c(F.score, F.p.value))
  }))
  print(time)
  results <- as.data.frame(do.call(rbind, results))
  row.names(results) <- row.names(exprs)
  colnames(results) <- c("F.score", "F.p.value")
  return(results)
}

#' Calculate all summary statistics
#' @param exprs A numeric matrix of gene expression data
#' @param partitions The set of self-contained conditions to sample id mappings
#' @param pairs A non-redundant pairing of partitions to compare
#' @param one.vs.all Boolean flag, indicating whether to compute complementary class statistics
#' @param means A data.frame containing the mean for each gene for each partition. Optional
#' @param std.devs A data.frame containing the standard deviation for each gene for each partition. Optional
#' @param medians A data.frame containing the median for each gene for each partition. Optional
#' @param mads A data.frame containing the mad for each gene for each partition. Optional
#' @param fold.changes A data.frame containing the fold change for each gene for each partition pair. Optional
#' @param parallel Boolean flag, indicating whether to try to use the registered parallel backend
#' @return A \code{list} of \code{data.frame}s, one for each summary statistic group computed.
do.summaries <- function(exprs, partitions, pairs, one.vs.all, means=NULL, std.devs=NULL, medians=NULL, mads=NULL, fold.changes = NULL, parallel=F,...){
  verbose("Calculating summary statistics")
  results <- list()
  summaries <- list()
  
  # Check if each simple summary statistic
  if(is.null(means)){
    summaries <- c(summaries, means=do.mean)
  } else {
    results <- c(results, list(means = means))
  }
  # Fast enough to recompute means, no need to make second order
  if(is.null(std.devs)){
    summaries <- c(summaries, std.devs=do.std.dev)
  } else {
    results <- c(results, list(std.devs=std.devs))
  }
  if(is.null(medians)){
    summaries <- c(summaries, medians=do.median)
  } else {
    results <- c(results, list(medians=medians))
  }
  
  hold <- llply(summaries, function(fn, exprs, partitions, one.vs.all){
    fn(exprs, partitions, one.vs.all)
  }, exprs, partitions, one.vs.all, .parallel=parallel)
  
  for(i in names(hold)){
    results[[i]] = hold[[i]]
  }
  
  summaries <- list()
  
  # Check each second order summary statistic
  if(is.null(fold.changes)){
    summaries <- c(summaries, fold.changes=function(exprs, partitions, pairs, within.summaries, one.vs.all){
      return(do.fold.change(within.summaries[['means']], pairs))
    })
  } else {
    results <- c(results, list(fold.changes=fold.changes))
  }
  if(is.null(mads)){
    summaries <- c(summaries, mads=function(exprs, partitions, pairs, within.summaries, one.vs.all){
      return(do.mad(exprs,within.summaries[["medians"]],partitions, one.vs.all))
    })
  } else {
    results <- c(results, list(mads=mads))
  }
  
  hold <- llply(summaries, function(fn, exprs, partitions, pairs, within.summaries, one.vs.all){
    fn(exprs, partitions, pairs, within.summaries, one.vs.all)
  }, exprs, partitions, pairs, results, one.vs.all, .parallel=parallel)
  
  for(i in names(hold)){
    results[[i]] = hold[[i]]
  }  
  
  return(results)
}
