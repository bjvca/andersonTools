% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matStand.R
\name{matStand}
\alias{matStand}
\title{Standardize Columns of a Matrix by Control Group}
\usage{
matStand(x, sgroup = rep(TRUE, nrow(x)))
}
\arguments{
\item{x}{A numeric matrix of data to be standardized.}

\item{sgroup}{A logical vector indicating the control group for standardization.}
}
\value{
A matrix with standardized columns.
}
\description{
Standardizes each column of a matrix by subtracting the mean and dividing by the standard deviation of the control group.
}
\examples{
mat <- matrix(rnorm(20), nrow = 5)
sgroup <- c(TRUE, TRUE, FALSE, TRUE, FALSE)
matStand(mat, sgroup)
}
