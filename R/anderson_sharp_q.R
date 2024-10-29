#' Anderson Sharp q-value Calculation
#'
#' Computes the Benjamini-Krieger-Yekutieli (2006) sharpened q-values for a vector of p-values.
#'
#' @param pval A numeric vector of p-values.
#' @return A numeric vector of sharpened q-values in the same order as the input.
#' @examples
#' pval <- c(0.01, 0.05, 0.02, 0.1)
#' anderson_sharp_q(pval)
#' @export
anderson_sharp_q <- function(pval) {
  # Number of p-values
  totalpvals <- length(pval)

  # Storing original sorting order
  original_sorting_order <- 1:totalpvals

  # Sorting indices and ranking
  sorted_indices <- order(pval)
  rank <- rank(pval)

  # Initial q-value
  qval <- 1
  bky06_qval <- rep(1, totalpvals)

  while (qval > 0) {
    # First Stage
    qval_adj <- qval / (1 + qval)
    fdr_temp1 <- qval_adj * rank / totalpvals
    reject_temp1 <- ifelse(fdr_temp1 >= pval, 1, 0)
    reject_rank1 <- reject_temp1 * rank
    total_rejected1 <- max(reject_rank1, na.rm = TRUE)

    # Second Stage
    qval_2st <- qval_adj * (totalpvals / (totalpvals - total_rejected1))
    fdr_temp2 <- qval_2st * rank / totalpvals
    reject_temp2 <- ifelse(fdr_temp2 >= pval, 1, 0)
    reject_rank2 <- reject_temp2 * rank
    total_rejected2 <- max(reject_rank2, na.rm = TRUE)

    # Update q-values for rejected hypotheses
    bky06_qval[rank <= total_rejected2] <- qval

    # Decrease qval for the next iteration
    qval <- qval - 0.001
  }

  # Restore original sorting order
  sorted_indices_inverse <- order(original_sorting_order)
  bky06_qval <- bky06_qval[sorted_indices_inverse]

  # Return the sharpened q-values
  return(bky06_qval)
}

# Example usage
# pval <- c(0.01, 0.05, 0.02, 0.1)
# sharpened_qvals <- anderson_sharp_q(pval)
# print(sharpened_qvals)
