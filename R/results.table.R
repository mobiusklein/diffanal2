#' @import RJSONIO

#' @title diffanal Results Table
#' 
#' Simple wrapper on top of \code{data.frame}.
#' 
#' @param obj A data.frame or data.table object
#' @param ngenes A record of the ngenes parameter from the call. \code{attr}
#' @param unique A record of the unique paramerter from the call. \code{attr}
#' @param column.groups A list of the different subclasses of columns. \code{list("p.value", "adj.p.value", "fold.changes", "means", "std.devs")}. \code{attr}
#' @param call A symbolic record of the call parameters. Useful for merging and writing reports. \code{attr}
#' @export
diffanal.results.table <- function(obj, ngenes, unique, column.groups, call = NULL){
    
  # Attach S3 class descriptor
  class(obj) <- c('diffanal.results', class(obj))
  
  # Attach attributes
  attr(obj, 'ngenes') <- ngenes
  attr(obj, 'unique') <- F
  attr(obj, 'column.groups') <- column.groups
  attr(obj, 'call') <- call
  
  return(obj)
}

NON.CLASS.COLS <- list("F.stat", "description", "gene")

#' @export
`column.groups` <- function(x){
  UseMethod('column.groups', x)
}

#' @export
`column.groups.diffanal.results` <- function(x){
  return(attr(x, 'column.groups'))
}

#' @export
`column.groups<-` <- function(x, value){
  UseMethod("column.groups<-", x)
}

#' @export
`column.groups<-.diffanal.results` <- function(x, value){
  attr(x, 'column.groups') <- value
}

#' @export
group <- function(x, col.grp.names){
  UseMethod("group", x)
}

#' @export
group.diffanal.results <- function(x, col.grp.names){
  grp.cols <- column.groups(x)
  grp.cols.sub <- grp.cols[col.grp.names]
  x <- x[,unlist(grp.cols.sub)]
  return(x)
}

#' @export
as.data.frame.diffanal.results <- function(x, row.names=NULL, ...){
  obj<-as.data.frame.data.frame(x, row.names, list(...))
  class(obj) <- c('diffanal.results', class(obj))
  attr(obj, 'column.groups')  <- column.groups(x)
  attr(obj, 'unique') <- attr(x, 'unique')
  attr(obj, 'call') <- attr(x, 'call')
  return(obj)
}

#' @export
as.unique <- function(obj, criteria.column="p.value", score.fn = which.min){
  if(!("diffanal.results" %in% class(obj))) stop("Object must be of type `diffanal.results`")  
  if(!(criteria.column %in% names(column.groups(obj)))) stop(paste("Criteria Column", criteria.column ,"not found in", do.call(paste, as.list(c("[", names(column.groups(obj)), "]")))))
  if(attr(obj, 'unique')){
    stop("This table has already been made unique, and each gene assigned a class.")
  }
  obj <- as.data.frame(obj)
  criteria.values <- as.data.frame(group(obj, criteria.column))
  compr.classes <- unlist(strsplit(x=column.groups(obj)[[criteria.column]],
                                   split=paste(JOIN.CHR, 
                                               sub(pattern="s$", '', criteria.column), 
                                               sep='')
                                   ))
  best.cases <- apply(criteria.values, 1, score.fn)
  best.class <- compr.classes[best.cases]
  best.class.cols <- lapply(column.groups(obj), function(category){
    cols <- as.data.frame(obj[,category]);
    best.cols <- sapply(seq_along(best.cases), function(i){
      return(cols[i,best.cases[i]])
    })
    return(best.cols)
  })
  unique.table <- do.call(cbind, best.class.cols)
  unique.table <- data.frame(obj[,"gene"], best.class, unique.table)
  group.names <- names(column.groups(obj))
  column.names <- paste('class', gsub("s$", "", group.names), sep=JOIN.CHR)
  names(unique.table) <- c("gene", "class", column.names)
  groups <- list()
  for(i in seq_along(group.names)){
    groups[[ group.names[[i]] ]] <- column.names[i]
  }
  unique.table <- diffanal.results.table(unique.table, attr(obj, 'ngenes'), T, groups, 
                                         attr(obj, 'call'))
  return(unique.table)
}


as.json <- function(table, ...){
  json.obj <- list()
  json.obj$columns.groups <- (column.groups(table))
  json.obj$unique <- (attr(table, 'unique'))
  json.obj$data <- alply(table, 1, function(x)x)
  names(json.obj$data) <- NULL
  json.obj <- toJSON(json.obj)
  return(json.obj)
}

to.html <- function(table, ...){
    
}


sort.table <- function(table, columns, ngenes = 250){
  
  sorts <- lapply(columns, function(column){
    order(table[,column])
  })
  
  # Merge the assembled gene list, retaining @ngenes from each set. At most the set will
  # contain number of comparisons
  top.genes <- Reduce(function(a, b){
    a <- unlist(a)
    b <- unlist(b)
    return(union(a, b[1:ngenes]))
  }, sorts[-1], sorts[[1]][1:ngenes])
  
  #print(top.genes)
  
  table <- table[unlist(top.genes),]
  return(table)
}
