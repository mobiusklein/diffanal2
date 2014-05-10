#' @import genefilter
#' @import doMC
#' @import foreach

library(genefilter)
library(doMC)

my.levels <- function( x, sort=T )
{
  #if ( !is.numeric(x) ) stop( "x must be numeric" )
  return( if (sort) sort(unique(as.vector(x)),na.last=T) else unique(as.vector(x)) )
}

my.tabulate <- function( x,nbins=length(unique(x)) )
{
  tabulate(match(x,unique(x)),nbins=nbins)
}

many2one <- function( cls )
{
  if ( is.null(dim(cls)) || ncol(cls)==1 ) {
    cls2 <- as.numeric(match(drop(cls),sort(unique(drop(cls)))))
    levels(cls2) <- levels(cls)
    return(cls2)
  }
  # ELSE ..
  #
  levs <- sort(apply(unique(cls),1,paste,collapse="."))
  cls.new <- as.numeric( match( apply(cls,1,paste,collapse="."), levs) )
  levels(cls.new) <- levs
  cls.new
}

# Creates permutations of binary label vectors
# @export
permute.binarray <- function( cls, nperm=1, balanced=F, equalized=F, seed=NULL,
                                  exhaustive=F, control=NULL)
{
  if ( !is.null(dim(cls)) )
    if (ncol(cls)==2)
      cls <- cls[,2]
  else
    stop( "multi-dimensional class label must have 2 columns" )
  lev <- my.levels(cls)
  m  <- length(cls)
  perm <- matrix( NA, nrow=nperm, ncol=length(cls) )
  if ( !is.null(seed) ) set.seed(seed)
  
  # non-binary response variable
  #
  if ( length(lev)>2 )
  {
    if ( is.null(control) ) {
      for ( i in (1:nperm) ) {
        perm[i,] <- sample(cls,m)
      }
      return(perm)
    }
    # control for confounder
    #
    else
    {
      control <- many2one(control)
      ctl.levs <- sort(unique(control))
      
      verbose(paste("controlling for ", length(ctl.levs), "-level confounder\n", sep="" ))
      
      idx <- lapply( ctl.levs, function(z) which(control==z) )
      jar <- lapply( idx, function(z) cls[z] )
      
      nperm.tot <- 0
      for ( i in 1:length(jar) ) {
        nperm.tot <- nperm.tot + (tmp <- lmchoose(my.tabulate(jar[[i]])))
        verbose( paste("  nchoose[",ctl.levs[i],"]:   ", exp(tmp), "\n", sep="" ))
      }
      verbose( paste("nchoose[tot]:", exp(nperm.tot), "\n" ))
      
      if ( log(nperm)>nperm.tot )
        warning( "more permutation than available" )
      
      for ( j in 1:length(idx) )
        perm[,idx[[j]]] <-
        t(apply(matrix(jar[[j]],nperm,length(jar[[j]]),byrow=T),1,sample))
    }
    return(perm)
  }
  # binary response variable
  #
  n1 <- sum(cls==lev[1])
  n2 <- sum(cls==lev[2])
  i.max <- which.max(c(n1,n2))
  max.lev <- lev[i.max]
  min.lev <- lev[-i.max]
  idx.max <- (1:m)[cls==max.lev]
  idx.min <- (1:m)[cls==min.lev]
  n.max <- c(n1,n2)[i.max]
  n.min <- c(n1,n2)[-i.max]
  
  if ( balanced )
  {
    if ( exhaustive ) {
      stop( "exhaustive not implemented for balanced permutations yet" )
    }
    perm[,idx.min] <- min.lev
    if ( length(lev)>2 )
      stop( "balanced permutation requires binary class" )
    if ( n1!=n2 && !equalized )
      stop( "balanced permutation requires same-size classes (or 'equalized=T')")
    k1 <- n.min%/%2
    k2 <- n.min%%2
    
    for ( i in 1:nperm )
    {
      k <- k1 + rbinom(1,1,.5)*k2
      idx.max.perm <- if ( n1!=n2 ) sample(idx.max,n.min) else idx.max
      perm[i,idx.max.perm] <- max.lev
      idx.max.perm <- sample( idx.max.perm, k )
      idx.min.perm <- sample( idx.min, 2*k-k )
      perm[i,idx.max.perm] <- min.lev
      perm[i,idx.min.perm] <- max.lev
    }
  }
  # IF NOT controlling for confounding variables
  #
  else if ( is.null(control) )
  {
    if ( exhaustive )
    {
      if ( choose(m,n.max)<=nperm )
      {
        nperm <- round(choose(m,n.max))
        verbose( "Number of exhaustive permutations:", nperm, "\n" )
        nidx <- t(combn(m, n.max))
        perm  <- matrix( min.lev, nrow=nrow(nidx), ncol=m )
        for ( i in 1:nrow(nidx) ) {
          perm[i,nidx[i,]] <- max.lev
        }
        return( perm )
      }
      else {
        verbose("Number of exhaustive permutations is too large:",
              choose(m,n.max),"(ignored)\n")
      }
    }
    # ELSE (if not exhaustive) ..
    #
    if ( log(nperm)>lchoose(n1+n2,n1)-log(2) )
      warning( paste("Number of possible permutations less than needed (", nperm, "): ",
                     choose(n1+n2,n1)/2, sep="") )
    for ( i in (1:nperm) )
    {
      idx <- if (equalized) 
        c(idx.min, sample(idx.max,n.min))
      else
        1:m
      perm[i,idx] <- sample(cls[idx],length(idx))
    }
  }
  # ELSE (controlling for confounding covariates) ..
  #
  else {
    lev0 <- lev[1]
    lev1 <- lev[2]
    control <- many2one(control)
    levs <- sort(unique(control))
    idx <- lapply( levs, function(z){which(control==z)} )
    nctl <- tabulate( control[cls==lev1], nbins=nlevels(control) )
    
    nperm.tot <- 0
    verbose("\n\tcontrolling for ", length(levs), "-level confounder\n", sep="" )
    
    # check if enough permutations available
    #
    for ( i in 1:length(nctl) ) {
      nperm.tot <- nperm.tot + (tmp <- lchoose(length(idx[[i]]),nctl[i]))
      verbose("\t  nchoose[",levels(control)[i],"]:   ", exp(tmp), "\n", sep="" )
    }
    verbose("\tnchoose[tot]:", exp(nperm.tot), "\n" )
    
    idx[nctl==0] <- NULL
    nctl <- nctl[nctl>0]
    if (nperm.tot<log(nperm)) {
      verbose("\tNumber of possible permutations less than required, reducing.\n" )
      nperm <- exp(nperm.tot)
      perm <- perm[1:nperm,]
    }    
    for ( i in 1:length(nctl) ) idx[[i]] <- c( nctl[i], idx[[i]] )
    
    for ( i in 1:nperm )
    {
      pidx <- unlist(sapply( idx, function(z){ sample( z[-1], z[1] ) } ))
      perm[i,] <- lev0
      perm[i,pidx] <- lev1
    }
  }
  return(perm)
}


#' @param exprs gene-by-sample numeric expression matrix
#' @param perm.labels A matrix of row of permuted group mapping factors
#' @return A numeric matrix of t-scores for each gene for each permutation
perm.tscore <- function(exprs, perm.labels){
  perm.scores <- alply(perm.labels, 1, function(labels, exprs){
    rowttests(exprs, fac=as.factor(labels), tstatOnly=T)['statistic']
  }, exprs)
  perm.scores <- do.call(cbind, perm.scores)
  names(perm.scores) <- seq_along(perm.scores)
  return(perm.scores)
}


#' @param exprs .
#' @param pheno.model.frame .
#' @param model.formula .
#' @param one.vs.all .
#' @param nperm .
#' @param exhaustive . 
#' @param ncores . 
#' @param summary.transform .
#' @param alpha .
#' @param ngenes .
#' @param smoother .
#' @param alternative .
#' @return diffanal.results.table instance
diffanal.perm.t.test <- function(exprs, pheno.model.frame, model.formula, one.vs.all = F, 
                      nperm=1, exhaustive=T, ncores=1, summary.transform = unlog,
                      alpha = 0.05, ngenes = 250, smoother = 1,
                      alternative = c("two.sided", "less", "greater"),
                      
                      ...){
  model.terms <- get.model.formula.terms(model.formula)
  predictor <- model.terms[1]
  confounders <- model.terms[-1]
  alternative <- match.arg(alternative)
  if(length(model.terms) == 1){
    confounders <- NULL
  }
  parallel = F
  if(ncores > 1 && require(doMC)){
    registerDoMC(ncores)
    parallel = T
  } else {
    registerDoSEQ()
  }
  # Error check the confounders
  sapply(confounders, function(col){
    if(!length(my.levels(pheno.model.frame[,col]) == 2)){
      stop(paste(col, "is not binary. Non-binary permutation tests are not yet supported."))
    }
  })
  
  if(length(confounders) > 1){
    stop("Only one confounder is supported at this time.")
  }
  
  partitions <- compute.partitions(pheno.model.frame, predictor)
  pairs <- compute.partition.pairs(partitions, one.vs.all)
  results <- alply(pairs, 1, function(pair, one.vs.all, alternative){
    verbose("Handling perms of", pair[1], "and", pair[2], "\n")
    # Reduce phenotype labels to binary comparison
    pheno.subset <- NULL
    obs.labels <- NULL
    if(one.vs.all){
      pheno.subset <- pheno.model.frame
      obs.labels <- as.factor(pheno.subset[,predictor] == unlist(pair[1]))
    } else{ 
      pheno.subset <- rbind(pheno.model.frame[pheno.model.frame[,predictor] == pair[1],],
                            pheno.model.frame[pheno.model.frame[,predictor] == pair[2],])
      obs.labels <- as.factor(as.character(pheno.subset[,predictor]))
    } 
    
    # Mirror sample reduction in expression data
    exprs.subset <- exprs[,pheno.subset$.rownames]
    confounders.subset <- NULL
    if(!is.null(confounders)){
      confounders.subset <- as.factor(pheno.subset[,confounders])
    }
    
    obs <- rowttests(exprs.subset, obs.labels)
    obs.score <-obs[,'statistic']
    obs.p.values <- obs[,'p.value']
    
    # Build set of label permutations
    perm.labels <- permute.binarray(obs.labels, nperm=nperm, 
                                    exhaustive=exhaustive, 
                                    control = confounders.subset)
    perm.scores <- perm.tscore(exprs.subset, perm.labels)
    
    # Compute the Permutation p-value. Computes the upper tail score. 
    perm.p.values <- sapply(seq(1,nrow(perm.scores)), function(i, perm.scores, obs.score, 
                                                               alternative, smoother){
      perm.scores <- perm.scores[i,]
      obs.score <- obs.score[i]
      # Sum gets the count of boolean TRUE
      perm.p <- NULL
      if(alternative == "greater") {
        perm.p <- (sum(perm.scores > obs.score) + smoother)/(length(perm.scores) + smoother)
      }
      if(alternative == "less") {
        perm.p <- (sum(perm.scores < obs.score) + smoother)/(length(perm.scores) + smoother)
      }
      if(alternative == "two.sided") {
        perm.p <- (sum(abs(perm.scores) > obs.score) + smoother)/(length(perm.scores) + smoother)
      }
      
      return(perm.p)
    }, perm.scores, obs.score, alternative, smoother)
    names(perm.p.values) <- row.names(exprs.subset)
    max.t <- max.col(perm.scores)
    max.t <- sapply(seq_along(max.t), function(ti){
      perm.scores[ti, max.t[ti]]
    })
    names(max.t) <- row.names(exprs.subset)
    
    return(list(obs.p.values=obs.p.values,obs.t=obs.score,
         perm.p.values=perm.p.values, max.t=max.t))
  }, one.vs.all, alternative, .parallel = parallel)

  obs.p.values <- lapply(results, function(test){test[["obs.p.values"]]})
  obs.p.values <- do.call(cbind, obs.p.values)
  colnames(obs.p.values) <- apply(pairs, 1, function(row){paste(row[1], row[2], 
                                                                'obs.p.value', sep=JOIN.CHR)})
  obs.t.scores <- lapply(results, function(test){test[["obs.t"]]})
  obs.t.scores <- do.call(cbind, obs.t.scores)
  colnames(obs.t.scores) <- apply(pairs, 1, function(row){paste(row[1], row[2], 
                                                                'obs.t.score', sep=JOIN.CHR)})
  
  max.t.scores <- lapply(results, function(test){test[["max.t"]]})
  max.t.scores <- do.call(cbind, max.t.scores)
  colnames(max.t.scores) <- apply(pairs, 1, function(row){paste(row[1], row[2], 
                                                                'max.t.score', sep=JOIN.CHR)})  
  
  perm.p.values <- lapply(results, function(test){test[["perm.p.values"]]})
  perm.p.values <- do.call(cbind, perm.p.values)
  colnames(perm.p.values) <- apply(pairs, 1, function(row){paste(row[1], row[2], 
                                                              'perm.p.value', sep=JOIN.CHR)})
  
  
  signif.comparisons <- which.significant(perm.p.values, alpha)
  if(length(signif.comparisons[signif.comparisons]) == 0){
    stop("No genes were significant")
  }
  obs.p.values <- obs.p.values[signif.comparisons, ]
  obs.t.scores <- obs.t.scores[signif.comparisons, ]
  perm.p.values <- perm.p.values[signif.comparisons, ]
  max.t.scores <- max.t.scores[signif.comparisons, ]
  F.stats <- do.aov(exprs[signif.comparisons, ], 
                               pheno.model.frame, model.formula)
  
  # Summary statistics
  signif.exprs <- exprs[signif.comparisons, ]
  transform.exprs <- summary.transform(signif.exprs)
  results <- do.summaries(transform.exprs, partitions, pairs, one.vs.all)
  
  # Construct Results Table
  results <- c(list(obs.p.values=obs.p.values, obs.t.scores=obs.t.scores, 
                    perm.p.values=perm.p.values, max.t.scores=max.t.scores), 
                    results, list(F.stats=F.stats))
  print(names(results))
  results.table <- do.call(cbind, results)
  names(results.table) <- c(unlist(lapply(results, colnames)))
  results.table <- cbind(gene = row.names(results.table), results.table)
  
  results.table <- sort.table(results.table, colnames(results$perm.p.values), ngenes)
  results.table <- diffanal.results.table(results.table, ngenes=ngenes, unique=F, column.groups=lapply(results, colnames))
  
  return(results.table)
}