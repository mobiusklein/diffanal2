
#' @title Filter Noisy Genes across a Dataset
#' @param dat An ExpressionSet object to be filtered
#' @param score One of "mad", "sd", or "cv". Scoring is done with the statistic given statistic.
#' @param dir One of "top" or "bottom", controlling which section of the score-space to select from.
#' @param transform Either a function that operates on a matrix or one of "none","log2","exp2","log", or "exp", which will be applied to the gene expression data prior to scoring
#' @param ngenes The number of genes to select. 
#' @param min.score A number specifiying minimum score to be accepted.
#' @export
variation.filter <- function(dat,
                             score=c("mad","sd","cv"),
                             dir=c("top","bottom"),
                             transform=c("none","log2","exp2","log","exp"),
                             ngenes=NULL,
                             min.score=NULL,
                             min.log=1,
                             rnd=4,
                             do.plot=F,
                             pch=".",
                             lgnd.coord=1,
                             do.log=NULL,
                             qnt.lev=0.5,
                             min.qnt=-Inf,
                             no.affx=F)
{
  if (is.null(ngenes) && is.null(min.score) )
    stop( "must specify either ngenes or min.score" )
  if (!is.null(ngenes) && !is.null(min.score) )
    stop( "cannot specify both ngenes and min.score" )
  if (!is.null(ngenes) && ngenes>nrow(dat) )
    stop( "ngenes is too large" )
  if (min.log<=0)
    stop( "min.log must be positive: ", min.log )
  if ( class(dat)!='ExpressionSet' )
    stop( "ExpressionSet object expected: ",class(dat) )
  
  dir <- match.arg(dir)
  score <- match.arg(score)
  score.fun <- match.fun(score)

  ##
  # Changed from original to support passing functions as well as 
  # string labels.
  # 
  # The function must take care of any preprocessing needed.
  ## 
  if(is.function(transform)){
    exprs(dat) <- transform(exprs(dat))
  } 
  else {
    transform <- match.arg(transform)  
    if (transform=="log2" || transform=="log") { # threshold before log-transformation
      verbose( "Thresholding before log-transforming .. " )
      exprs(dat)[exprs(dat)<min.log] <- min.log
      verbose( "done.\n" )
    }
    exprs(dat) <- switch(transform,
                         none=exprs(dat),
                         log2=round(log2(exprs(dat)),rnd),
                         exp2=round(2^(exprs(dat)),rnd),
                         log=round(log(exprs(dat)),rnd),
                         exp=round(exp(exprs(dat)),rnd))
    
  }
  
  if (no.affx) {
    if ( length(rm.idx <- grep("AFFX-",featureNames(dat)))>0 ) {
      verbose("Removing 'AFFX-' probes ..")
      dat <- dat[-rm.idx,,drop=FALSE]
      verbose(" done,", length(rm.idx),"removed.\n")
    }
  }
  ctr <- if (score=="mad") apply( exprs(dat), 1, median ) else rowMeans( exprs(dat) )
  SC <- SC1 <- apply( exprs(dat), 1, score.fun )

  if (min.qnt>0)
  {
    verbose( "Filtering out genes w/ ",round(100*qnt.lev,2), "-percentile < ", min.qnt, " .. ",sep="" )
    QNT <- apply(exprs(dat),1,quantile,probs=qnt.lev)
    if ( sum(QNT>=min.qnt)<2 )
      stop( "filtering by min.qnt returns less than 2 genes (try decreasing min.qnt)" )
    dat <- dat[QNT>=min.qnt,,drop=FALSE]
    verbose( "done,", nrow(dat), "genes left.\n")

    if ( !is.null(ngenes) && nrow(dat)<=ngenes ) {
      verbose("Number of genes left is less than required, no further filtering necessary")
      return(dat)
    }
    SC1 <- SC[QNT>=min.qnt]
  }
  verbose( "Variation filtering based on", score, ".. " )
  
  verbose( "done.\n" )
  
  idx <- NULL
  if ( is.null(ngenes) ) {
    verbose( "Selecting genes with", score, "<=", min.score, ".. " )
    idx <- if(dir=="top") SC1>=min.score else SC1<=min.score
    if (sum(idx)==0)
      stop( "no genes passed the filtering criteria" )
  }
  else {
    verbose( "Selecting top", ngenes, "by", score, ".. " )
    if (dir=="top") SC1 <- -SC1
    idx <- order(SC1)[1:ngenes]
  }
  dat <- dat[idx,,drop=FALSE]
  verbose( "done,", nrow(dat), "genes selected.\n" )

  if (do.plot) {
    verbose( "Creating scatter plot .. ")
    if (is.null(do.log)) {
      do.log <- if (transform=="none" || transform=="exp2" || transform=="exp" )
        "xy"
      else
        ""
    }
    SC <- abs(SC)
    plot( ctr, SC, pch=pch, col="gray", xlab=if (score=="mad") "median" else "mean", ylab=score, log=do.log)
    plot.idx <- match( featureNames(dat),names(SC) )
    points( ctr[plot.idx], SC[plot.idx],pch=pch,col="red")
    lx <- min(ctr); ly <- max(SC); xjust <- 0; yjust <- 1
    if (lgnd.coord==2) {
      lx <- max(ctr); xjust <- 1
    }
    else if (lgnd.coord==3) {
      lx <- max(ctr); ly <- min(SC); xjust <- 1; yjust <- 0
    }
    else if (lgnd.coord==4) {
      ly <- min(SC); yjust <- 0
    }
    else if (lgnd.coord!=1)
      stop( "lgnd.coord must be btw 1 and 4" )
    
    legend(lx, ly, xjust=xjust, yjust=yjust,
           legend=c("all","passing filter"), col=c("gray","red"), pch=20)
    verbose( "done.\n")
  }
  dat
}
