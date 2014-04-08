require(plyr)
require(genefilter)
require(Biobase)
require(limma)
require(permute)

# Column names with aggregate meaning will have tokens joined by the JOIN.CHR
# e.g. cls1-cls2-p.value means the p.value comparing cls1 and cls2.
JOIN.CHR = "-"
SAMPLE.ID.COL <- '.rownames'
STRATEGY = list(t.test = function(){}, survival = function(){}, lm = function(){})


# model.formula should have no RHS
diffanal2 <- function(exprs, pheno, model.formula, ngenes = 250, alpha = 0.05, p.adjust.method = 'none', one.vs.all = F){
  print(deparse(formals()))
  print(match.call())
  model.formula.terms <- terms(model.formula)
  # Specify without response variable
  model.formula.cols <- row.names(attr(model.formula.terms, 'factors'))
  cleaned <-clean.phenotypes(exprs, pheno, model.formula.cols)
  exprs <- cleaned[["exprs"]]
  pheno <- cleaned[["pheno"]]
  model.cols <- prepare.pheno.model.frame(pheno, model.formula.cols)
  # Cartesian Product of cofactors separates the specific sample names
  partitions <- dlply(model.cols, model.formula.cols, function(set) set[,'.rownames'])
  pairs <- partition.pairs(partitions, one.vs.all = one.vs.all)
  scores <- do.t.test(exprs, partitions, pairs, one.vs.all = one.vs.all)
  means <- do.means(exprs, partitions)
  signif.pairs <- which.significant(scores, alpha)
  
  return(list(scores=scores, means=means, partitions=partitions, exprs=exprs, pheno=pheno))
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
  #model.cols <- cbind(model.cols, .rownames = row.names(model.cols))
  model.cols<-name_rows(model.cols)
  # column vectors with the name attribute breaks dlply  
  for(i in seq_along(model.cols)){
    names(model.cols[i]) <- NULL
    attr(model.cols[i], '.Names') <- NULL
  }
  return(model.cols)
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
    
    test.stats <- rowttests(exprs.subset, as.factor(colnames(exprs.subset) %in% partitions[[pairs[1]]]))
    if(statistic){
      return(test.stats["statistic"])
    } else{
      return(test.stats["p.value"]) 
    }    
  },exprs,partitions)
  
  results <- do.call(cbind, results)
  names(results) <- names(results) <- apply(pairs, 1, function(ps){paste(ps[[1]], ps[[2]], 'p.value', sep=JOIN.CHR)})
  
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
    min <- min(row)
    # If
    if(min < alpha){
      cols <- which(row == min)
      # If there are ties, this may be longer than 1. 
      return(pairs[cols][1])
    }
    return(NA)
  })
  results <- results[!is.na(results)]
  return(results)
}

do.limma.lmFit <- function(exprs, pheno, model.formula, ngenes=250, alpha){
  design <- model.matrix(model.formula, pheno)
  fit <- lmFit(exprs, design)
  fit.ebayes <- eBayes(fit)
  #fit.top <- topTable(fit.ebayes, number=ngenes)
  return(fit.ebayes)
}