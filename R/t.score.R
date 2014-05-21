t.score <- function( x, cls=NULL, y=NULL, robust=F, paired=F, generalized=F,
                     do.test=F, var.equal=F, min.sd=NULL, fix=F,
                     alternative=c("two.sided","greater","less"),cls1=NULL, cls2=NULL )
{
  # INPUT:
  #    x - m x n1 matrix (genes by experiments, condition 1)
  #    y - m x n2 matrix (genes by experiments, condition 2)
  #  OR
  #    x - m x n  matrix (genes by experiments, condition 1 & 2)
  #  cls - n vector of class labels
  #
  # OUTPUT:
  #  snr - m vector (positive are upregulated for x, or for
  #        lower label -- condition 1 -- when cls is specified)
  
  # some checks on the input
  #
  alternative <- match.arg(alternative)
  
  if ( is.null(y) & is.null(cls) )
    stop( "must specify either y or cls" )
  
  if ( is.null(y) )
  {
    lev <- sort(unique(cls))
    if ( ncol(x)!=length(cls) )
      stop( "ncol(x) must be same as length(cls)" )
    if ( length(lev)>2 )
      stop( "cls must be binary" )
    y <- x[,cls==lev[2],drop=F]
    x <- x[,cls==lev[1],drop=F]
  }  
  if ( nrow(x)!=nrow(y) ) stop( "x and y must be of same length\n" )
  if ( ncol(x)<4 ) warning( "x has less than 4 observations\n" )
  if ( ncol(y)<4 ) warning( "y has less than 4 observations\n" )
  
  score <- NULL
  
  # paired score
  #
  if ( paired )
  {
    if ( generalized )
    {
      if ( (ncol(x)>ncol(y) & ncol(x)%%ncol(y)) |
             (ncol(y)>ncol(x) & ncol(y)%%ncol(x)) )
        stop( "x and y must be have column numbers multiple of each other\n" )
      
      d <- NULL
      
      if ( ncol(x)>ncol(y) )
        for ( i in 0:(ncol(x)/ncol(y)-1) ) {
          idx <- (i*ncol(y)+1):((i+1)*ncol(y))
          verbose(, "comparing x[:",idx[1],":",idx[length(idx)],
                  "] to y[1:",ncol(y),"]","\n",sep="")
          d <- cbind( d, x[,idx,drop=F]-y )
        }
      else if ( ncol(y)>ncol(x) )
        for ( i in 0:(ncol(y)/ncol(x)-1) ) {
          idx <- (i*ncol(x)+1):((i+1)*ncol(x))
          verbose(, "comparing x[1:",ncol(x),"] to y[",
                  idx[1],":",idx[length(idx)],"]","\n",sep="")
          d <- cbind( d, x-y[,idx,drop=F] )
        }
      else
        d <- x-y
    }
    else
    {
      if ( ncol(x)!=ncol(y) )
        stop( "x and y must have same number of columns\n" )
      
      d <- x-y
    }
    if ( robust ) {
      stop( "robust paired t.score not implemented yet" )
    }
    else {
      if (do.test) {
        score <- cbind(score=drop((fast.mean(d)*sqrt(length(cls)))/fast.sd(d)),
                       p.value=drop(apply(d, 1, function(z) t.test(z)$p.value)))
        rownames(score) <- rownames(d)
      }
      else {
        score <- (drop(fast.mean(d))*sqrt(ncol(d)))/drop(fast.sd(d))
        names(score) <- rownames(d)
      }
    }
    return( score )
  }
  # ELSE not paired
  #
  x.idx <- 1:ncol(x)
  
  n1 <- (ncol(x)); if (n1<2) stop( "need at least 2 obs per class" )
  n2 <- (ncol(y)); if (n2<2) stop( "need at least 2 obs per class" )
  cls <- c( rep(1,n1), rep(0,n2) ); cls <- cbind( cls, 1-cls )
  x <- cbind(x,y)
  
  if ( robust )
  {
    rnk <- t(apply(-x,1,rank))
    score <-  (drop(rnk[,1:n1] %*% rep(1,n1)) - as.double(n1) * (as.double(n1)+1)/2)
    
    if (do.test) {
      warning( "wilcox.test p-value not implemented yet, ignoring" )
    }
  }
  else
  {
    s  <- x %*% cls
    s2 <- x^2 %*% cls
    s2[,1] <- (s2[,1] - (s[,1]^2)/n1) / (n1-1) # variance in 1st class 
    s2[,2] <- (s2[,2] - (s[,2]^2)/n2) / (n2-1) # variance in 2nd class 
    s[,1] <- s[,1]/n1
    s[,2] <- s[,2]/n2
    
    if ( fix )
    {
      s2[,1] <- fix.sd(sqrt(s2[,1]), s[,1])^2
      s2[,2] <- fix.sd(sqrt(s2[,2]), s[,2])^2
    }
    if ( !is.null(min.sd) ) {
      min.var <- min.sd*min.sd
      s2[s2[,1]<min.var,1] <- min.var
      s2[s2[,2]<min.var,2] <- min.var
    }    
    stderr <- if (var.equal)
      sqrt( (((n1-1)*s2[,1] + (n2-1)*s2[,2])/(n1+n2-2)) * (1/n1+1/n2) )
    else
      sqrt( s2[,1]/n1 + s2[,2]/n2 )
    
    score <- (s[,1]-s[,2]) / stderr
    
    results <- NULL
    if ( do.test )
    {
      df <- if ( var.equal ) # degrees of freedom
        n1+n2-2
      else
        stderr^4 / ( (s2[,1]/n1)^2/(n1-1) + (s2[,2]/n2)^2/(n2-1)) # Welch approximation of df
      
      pval <- if (alternative == "less") {
        pt(score, df=df)
      }
      else if (alternative == "greater") {
        pt(score, df=df, lower.tail=F)
      }
      else {
        2 * pt(-abs(score), df=df)
      }
      results <- data.frame( score=score, p.value=pval )
    } else{
      results <- data.frame(score = score)
    }
      
  }
  return( results )
}