#' Calculate Anderson Index
#'
#' Creates an index as a weighted average of outcomes, where weights are based on the inverse covariance matrix of outcomes.
#'
#' @param xmat A matrix of outcomes.
#' @param revcols An optional vector indicating columns to reverse signs for "better" outcomes.
#' @param sgroup A logical vector indicating the control group for standardization.
#' @return A list containing the index and weights.
#' @examples
#' xmat <- matrix(rnorm(20), nrow = 5)
#' Anderson_Index(xmat)
#' @export
anderson_index <- function(	xmat,
                            #wgts=rep(1, nrow(xmat)), #nrow: number of rows present in xmat --> many 1s
                            revcols = NULL,
                            sgroup = rep(TRUE, nrow(xmat))){
  X <- matStand(xmat, sgroup)
  if(length(revcols)>0){
    X[,revcols] <-  -1*X[,revcols]
  }
  i.vec <- as.matrix(rep(1,ncol(xmat)))
  #Sx <- cov.wt(X, wt=wgts)[[1]]
  #list with estimates of the weighted covariance matrix and the mean of the data
  Sx <- cov(X,use = "pairwise.complete.obs")
  #cov: covariance of x and y if these are vectors/covariances between columns of x and columns of y are computed if these are matrices
  #use = "everything" produces NAs for the index.
  #use = "all.obs" produces an error.
  #use = "complete.obs" and use = "na.or.complete": works, NAs are handled by casewise deletion.
  #use = "pairwise.complete.obs": works, covariance between each pair of variables is computed using all complete pairs of observations on those variables
  weights <- solve(t(i.vec)%*%solve(Sx)%*%i.vec)%*%t(i.vec)%*%solve(Sx)
  index <- t(solve(t(i.vec)%*%solve(Sx)%*%i.vec)%*%t(i.vec)%*%solve(Sx)%*%t(X))
  return(list(weights = weights, index = index))
}
