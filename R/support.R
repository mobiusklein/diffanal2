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



