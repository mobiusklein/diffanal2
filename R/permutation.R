require(permute)

# data$confound is a column/vector of confounding class factors 
# n.perms is the number of permutations requested
# perm.set <- shuffleSet(nobs, control=how(plots=c(Plots(strata=data$confounds, type='free')), nperms = n.perms))
# 
# Copy the data labels to keep using
# perm.data <- data
# for each permutation i
# perm.data$pheno <- data$pheno[perm.set[i,]]
# ... do test with perm.data

# If data$confounds is not balanced, instantiating the perm.set will fail. 

build.permutation.set <- function(model.frame, n.perm = 200, exhaustive=F, 
                                  rnd.seed=Sys.time()){
  # Remove the predictor (first colum) and the row names column
  confounds <- names(model.frame)[-1]
  confounds <- confounds[confounds != '.rownames']
  
  plotSet = c()
  # For-Loop because any other form of aggrgation breaks the way 
  # Plots() collections work. 
  for(name in confounds){
    p <- (Plots(strata=model.frame[, name], type='free'))
    plotSet = append(plotSet, p)
  }
  
  set.seed(rnd.seed)
  
  # If exhaustive is some large number like 7.3e105, this function 
  # will throw an error. It might be possible to work around with
  # the maxperm parameter, but that is set to unused currently.
  control = how(plots=plotSet, nperm=n.perm, complete=exhaustive)
  pset <- shuffleSet(nobs(model.frame), control=control)
  return(pset)
}

t.score.fn <- function(exprs, model.frame, model.formula, one.vs.all = F){
  partitions <- compute.partitions(model.frame, model.formula)
  pairs <- partition.pairs(partitions, one.vs.all = one.vs.all)
  scores <- do.t.test(exprs, partitions=partitions, pairs=pairs, one.vs.all=one.vs.all, statistic=T)
  return(scores)
}

perm.test <- function(exprs, model.frame, model.formula, test.fn, n.perm = 200, exhaustive = F, rnd.seed = Sys.time(), one.vs.all = F){
  obs.score <- test.fn(exprs, model.frame, model.formula)
  perm.set <- build.permutation.set(model.frame, n.perm, exhaustive, rnd.seed)
  perm.scores <- alply(perm.set, 1, function(perm.inds, exprs, model.frame, moel.formula, test.fn, one.vs.all){
    perm.frame <- model.frame
    # Predictor variable is the first column in model.frame
    perm.frame[,1] <- model.frame[perm.inds, 1]
    # Returns a vector of scores
    results <- test.fn(exprs, perm.frame, model.formula, one.vs.all=one.vs.all)
    return(results)
  }, exprs, model.frame, model.formula, test.fn, one.vs.all)
  return(perm.scores)
}


shuffle.across.strata <- function(strata, size){
  # Duplicate strata
  strata <- as.data.frame(strata)
  out <- strata
  shifts <- sample(seq_along(strata[,1]), size)
  positions <- seq_along(out[,1])
  
  cmpr <- function(row, class, confounds){
    return(row[1] != class && all(row[-1] == confounds))
  }
  
  for(shift.ind in shifts){
    shifter <- out[shift.ind,]
    shifter.class <- shifter[,1]
    shifter.confounds <- shifter[,-1]
    print(shifter)
    valid.positions <- positions[apply(out, 1, cmpr, shifter.class, shifter.confounds)]
    #print(valid.positions)
    print(out[valid.positions,])
    new.position <- sample(valid.positions, 1)
    out[shift.ind, ] <- out[new.position, ]
    out[new.position, ] <- shifter
  }
  return(out)
}
