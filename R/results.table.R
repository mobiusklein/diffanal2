#' @title \code{diffanal} Results Table
#' 
#' @description A simple wrapper on top of \code{data.frame} that bundles columns into groups
#' 
#' @param obj A data.frame or data.table object
#' @param ngenes A record of the ngenes parameter from the call. \code{attr}
#' @param unique A record of the unique paramerter from the call. \code{attr}
#' @param column.groups A list of the different subclasses of columns. \code{list("p.value", "adj.p.value", "fold.changes", "means", "std.devs")}. \code{attr}
#' @param call A symbolic record of the call parameters. Useful for merging and writing reports. \code{attr}
#' @export
diffanal.results.table <- function(obj, ngenes = Inf, unique = F, column.groups, call = NULL){
    
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

#' @title column.groups
#' @param x An R Object
#' @export
`column.groups` <- function(x){
  UseMethod('column.groups', x)
}

#' @S3method column.groups diffanal.results
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

#' @S3method unique diffanal.results
#' @method unique diffanal.results
#' @export
unique.diffanal.results <- function(obj, criteria.column="p.value", score.fn = which.min){
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

# @importFrom RJSONIO toJSON
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

#' Perform row-wise comparison between a data.frame or matrix and a single vector
#' 
#' @param df A data.frame or matrix
#' @param match A vector that shares columns with \code{df}
#' @param cols Optionally, a character vector naming which columns to compare between \code{df} and \code{match}
row.eq <- function(df, match, cols){
  # Restrict to comparing only specific columns
  if(!is.null(cols)){
    df <- as.data.frame(df[,cols])
    match <- match[cols]
  }
  # row-wise compare to match
  apply(df, 1, function(row){
    all(row == match)
  })
}


#' @name  sort.diffanal.results
#' @title sort.diffanal.results
#' @description  Given an unordered results table, order it by the column group given, and report the top ngenes given in each column of the given column group
#' @param table A data.frame
#' @param column.grp A set of columns from \code{table} or a data.frame that is parallel to it.
#' @param ngenes The number of genes to return for each class.
#' @S3method sort diffanal.results
#' @method sort diffanal.results
#' @export
sort.diffanal.results <- function(table, column.grp, ngenes = 250){
  
  sorts <- lapply(column.grp, function(column){
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


#' @title which.significant
#' @description Determine which rows of a data.frame of p-values contain a significant score.
#' @param scores A data.frame or matrix of p.values
#' @param alpha p-value threshold for significance
#' @return A boolean vector labeling each row in \code{scores} as significant or not
#' @export
which.significant <- function(scores, alpha = 0.05){
  results <- apply(scores, 1, function(row){
    return(any(row < alpha))
  })
  return(results)
}

#' @title which.class.significant
#' @description Determine which rows of a data.frame of p-values contain a significant score, and which column the minimum score was found in. 
#' @param scores A data.frame or matrix of p.values
#' @param alpha p-value threshold for significance
#' @return A vector labeling each row in scores for which column is significant
which.class.significant <- function(scores, alpha){
  results <- apply(scores, 1, function(row){
    if(any(row < alpha)){
      return(which.min(row))
    }
  })
  return(results)
}