
## rescale
# Apply scale along each column of a matrix or data.frame or
# along a single vector. 
rescale <- function(x)
{
  if ( is.matrix(x) )
    apply( x, 2, scale )
  else if ( is.vector(x) )
    ( x-min(x) ) / ( max(x) - min(x) )
  else if ( is.data.frame(x) )
    sapply( x, scale )
  else
    scale( as.vector(x) )
}

## 
