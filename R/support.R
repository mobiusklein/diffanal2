#' @import plyr
#' @import genefilter

# Column names with aggregate meaning will have tokens joined by the JOIN.CHR
# e.g. cls1-cls2-p.value means the p.value comparing cls1 and cls2.
JOIN.CHR = ".."
SAMPLE.ID.COL <- '.rownames'


#' Print each argument to stdout if the verbose option or .diffanal.verbose option
#' has been set to true. Set these options using options(...).
verbose <- function(...){
  opt <- getOption("verbose", F) || getOption(".diffanal.verbose", F)
  if(!opt) invisible(NULL)
  msgs <- list(...)
  msgs <- do.call(paste, c(msgs, '\n'))
  cat(msgs)
}

#' Converts a list of terms into a formula object for building models with
formula.from.list <- function(predictor, confounders = list()){
  # Attach the predictor and tilde first
  relation <- paste("~", predictor)
  # Combine any confounders into a chain of additions
  verbose("Confounders: ", confounders)
  trailing.terms <- do.call(paste, c(confounders, sep = " + "))
  # If there are trailing terms to merge
  if(length(trailing.terms) > 0){
    # Attach the predictor and the confounder string
    relation <- paste(relation, trailing.terms, sep=" + ")
  }
  full.formula <- as.formula(relation)
  return(full.formula)
}

get.model.formula.terms <- function(model.formula){
  model.formula.terms <- terms(model.formula)    
  col.names <- row.names(attr(model.formula.terms, 'factors'))
  return(col.names)
}

#' The predictor is always stored in the first position of the formula's
#' terms object. This holds even if a 0 is prepended.
#' @param model.formula A \code{formula} object
#' @return The first non-zero term in model.formula
get.predictor <- function(model.formula){
  form.terms <- get.model.formula.terms(model.formula)
  return(form.terms[1])
}

#' Insure that an inercept term is not included.
#' @param model.formula A \code{formula} object
#' @return A copy of model.formula where the predictor has been offset by "0 +"
adjust.formula <- function(model.formula){
  form <- as.character(model.formula)
  form <- strsplit(form, "~")[[2]]
  form <- as.formula(paste("~ 0 +", form))
  return(form)
}

#' Takes the input matrix out of log2 space
unlog <- function(exprs){
  unlogged.exprs = 2 ** exprs
  return(unlogged.exprs)
}

#' Remove samples that are NA in any of the model dimensions
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

#' Cleans up the pheno.model.frame
#' @param pheno A data.frame containing predictor and confounder labels for each sample
#' @param model.formula.cols A vector of strings that label which columns from \code{pheno} to keep
prepare.pheno.frame <- function(pheno, model.formula.cols){
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

#' Create a class to sample name mapping from the model.frame.
#' @param pheno.model.frame A data.frame containing predictor and confounder labels for each sample
#' @param model.formula A formula describing the predictor and confounder labels
compute.partitions <- function(pheno.model.frame, model.formula){
  if(is.formula(model.formula)) model.formula.cols <- get.model.formula.terms(model.formula)
  else model.formula.cols <- model.formula
  # Cartesian Product of cofactors separates the specific sample names
  partitions <- dlply(pheno.model.frame, model.formula.cols, function(set) set[,'.rownames'])
  return(partitions)
}

#' Computes the complement of each partition
#' 
#' @param pheno.frame A model.frame that describes all of the columns on which a model is to be fit plus the .rownames column to label samples
#' @param partitions A list of each partition subset
complementary.indices <- function(pheno.frame, partitions){
  sample.names <- pheno.frame[,SAMPLE.ID.COL]
  complements <- lapply(partitions, function(partition){
    other <- !(sample.names %in% partition)
    other <- sample.names[other]
    return(other)
  })
  names(complements) <- paste('Not', names(complements), sep='.')
  return(complements)
}

#' Creates unique pairwise combinations of different conditions, 
#' or for a one-vs-all design, just the pairing of each condition 
#' with its complement.
#' 
#' @param partitions A list of each partition subset
#' @param one.vs.all A flag indicating whether to produce one-vs-all complementary sets
compute.partition.pairs <- function(partitions, one.vs.all=F){
  # If using a one.vs.all design, create a simple partition pair frame
  # with 
  if(one.vs.all){
    pairs <- as.data.frame(cbind(names(partitions), paste("Not", names(partitions), sep='.')))
    names(pairs) <- c('group1', 'group2')
    return(pairs)
  }
  pairs <- t(combn(names(partitions),2))
  return(pairs)
}



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
