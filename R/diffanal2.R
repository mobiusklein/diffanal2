source("./logger.R")

.DiffAnal <- setRefClass("DiffAnal", fields=list(logger = "Logger", exprs = "matrix", pheno = 'data.frame'))

DiffAnal <- function(exprs = matrix(), pheno = data.frame(), verbose = F, logfile = F){
  logger <- Logger(name="DiffAnal", verbose = verbose, logfile = logfile)
  .DiffAnal(exprs = exprs, pheno = pheno, logger = logger)
}

