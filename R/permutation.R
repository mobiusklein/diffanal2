require(permute)

# data$confound is a column/vector of confounding class factors 
# n.perms is the number of permutations requested
# perm.set <- shuffleSet(nobs, control=how(plots=c(Plots(strata=data$confounds, type='free')), nperms = n.perms))
# 
# Copy the data labels to keep using
# perm.data <- data
# for each permutation i
# perm.data$pheno <- data$pheno[perm.set[i,],]
# ... do test with perm.data

# If data$confounds is not balanced, instantiating the perm.set will fail. 
