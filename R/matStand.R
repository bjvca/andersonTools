#' Standardize Columns of a Matrix by Control Group
#'
#' Standardizes each column of a matrix by subtracting the mean and dividing by the standard deviation of the control group.
#'
#' @param x A numeric matrix of data to be standardized.
#' @param sgroup A logical vector indicating the control group for standardization.
#' @return A matrix with standardized columns.
#' @examples
#' mat <- matrix(rnorm(20), nrow = 5)
#' sgroup <- c(TRUE, TRUE, FALSE, TRUE, FALSE)
#' matStand(mat, sgroup)
#' @export
matStand <- function(x, sgroup = rep(TRUE, nrow(x))){
  for(j in 1:ncol(x)){
    x[,j] <- (x[,j] - mean(x[sgroup,j], na.rm = TRUE)) / sd(x[sgroup,j], na.rm = TRUE)
  }
  return(x)
}
