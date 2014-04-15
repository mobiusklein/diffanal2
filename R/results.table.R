require(RJSONIO)

diffanal.results.table <- function(obj, ngenes, unique, column.classes, call){
  obj <- as.data.table(obj)
  class(obj) <- c('diffanal.results', class(obj))
  attr(obj, 'ngenes') <- ngenes
  attr(obj, 'unique') <- unique
  attr(obj, 'column.classes') <- column.classes
  attr(obj, 'call') <- call
  return(obj)
}

to.html <- function(table, ...){
    
}