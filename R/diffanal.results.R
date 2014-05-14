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
diffanal.results <- function(obj, ngenes = Inf, unique = F, column.groups, call = NULL){
    
  # Attach S3 class descriptor
  class(obj) <- c('diffanal.results', class(obj))
  if(!('descriptives' %in% names(column.groups))){
    column.groups[['descriptives']] <- character()
  }
  column.groups[['descriptives']] <- unique(c(column.groups[['descriptives']], "gene"))
  # Attach attributes
  attr(obj, 'ngenes') <- ngenes
  attr(obj, 'unique') <- unique
  attr(obj, 'column.groups') <- column.groups
  attr(obj, 'call') <- call
  
  return(obj)
}

NON.CLASS.GRPS <- c("F.stats")

#' @rdname Column Grouping
#' @title column.groups
#' @param x An R Object
#' @export
`column.groups` <- function(x){
  UseMethod('column.groups', x)
}

#' @S3method column.groups diffanal.results
#' @method column.groups diffanal.results
`column.groups.diffanal.results` <- function(x){
  return(attr(x, 'column.groups'))
}

#' @rdname Column Grouping
#' @param x An R Object
#' @param value A character vector of column names in \code{x}
#' @export
`column.groups<-` <- function(x, value){
  UseMethod("column.groups<-", x)
}

#' @export
`column.groups<-.diffanal.results` <- function(x, value){
  attr(x, 'column.groups') <- value
}

#' @rdname Column Grouping
#' @param x An R object
#' @param col.grp.names A vector of group names, as those returned by \link{\code{column.groups}}
#' @return A subset of \code{x} containing only columns corresponding to those mapped by \code{col.grp.names}
#' @export
group <- function(x, col.grp.names){
  UseMethod("group", x)
}

#' @S3method group diffanal.results
group.diffanal.results <- function(x, col.grp.names){
  grp.cols <- column.groups(x)
  grp.cols.sub <- grp.cols[col.grp.names]
  sub.x <- as.data.frame(x[,unlist(grp.cols.sub)])
  names(sub.x) <- unlist(grp.cols.sub)
  return(sub.x)
}

#' @S3method as.data.frame diffanal.results
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
unique.diffanal.results <- function(obj, criteria.column="p.value", score.fn = which.min){
  if(!("diffanal.results" %in% class(obj))) stop("Object must be of type `diffanal.results`")  
  if(!(criteria.column %in% names(column.groups(obj)))) stop(paste("Criteria Column", criteria.column ,"not found in", do.call(paste, as.list(c("[", names(column.groups(obj)), "]")))))
  if(attr(obj, 'unique')){
    stop("This table has already been made unique, and each gene assigned a class.")
  }
  
  criteria.values <- as.data.frame(group(obj, criteria.column))
  compr.classes <- unlist(strsplit(x=column.groups(obj)[[criteria.column]],
                                   split=paste(JOIN.CHR, 
                                               sub(pattern="s$", '', criteria.column), 
                                               sep='')
                                   ))
  # Find the best class for each row for the given column group, ordered by score.fn
  best.cases <- apply(criteria.values, 1, score.fn)
  best.class <- compr.classes[best.cases]
  
  # Select the columns that are classed
  col.grps <- column.groups(obj)[!(names(column.groups(obj)) %in% 
                                     c(NON.CLASS.GRPS, "descriptives"))]
  unclassed.grps <- column.groups(obj)[NON.CLASS.GRPS]
  # Select all columns from each group belonging to that row's best class
  best.class.cols <- lapply(col.grps, function(category){
    cols <- as.data.frame(obj[,category]);
    best.cols <- sapply(seq_along(best.cases), function(i){
      return(cols[i,best.cases[i]])
    })
    return(best.cols)
  })
  
  unique.table <- do.call(cbind, best.class.cols)
  
  #Prepend descriptive columns
  unique.table <- data.frame(group(obj, "descriptives"), best.class, 
                             unique.table)
  # Add "class" descriptor to the classed columns
  group.names <- names(col.grps)
  column.names <- paste('class', gsub("s$", "", group.names), sep=JOIN.CHR)
  names(unique.table) <- c(column.groups(obj)[["descriptives"]], "class", 
                           column.names)
  # Rebuild the column.groups list using the class columns
  groups <- list()
  for(i in seq_along(column.names)){
    groups[[ group.names[[i]] ]] <- column.names[i]
  }
  
  # Re-attach unclassed columns to both the table and the column.groups
  unique.table <- cbind(unique.table, group(obj, NON.CLASS.GRPS))
  groups <- c(groups, unclassed.grps)
  
  # Repackage as diffanal.results instance
  unique.table <- diffanal.results(unique.table, attr(obj, 'ngenes'), unique=T, groups, 
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