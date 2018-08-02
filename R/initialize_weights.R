#===============================================================================
# initialize_weights.R
#===============================================================================

#' @title Initialize weights
#'
#' @details
#' Initialize the weights using the histogram of p_hat
#'
#' @param data_table Input data table
#' @param n_breaks Number of breaks for the histogram
#' @param spline_order Order of splines
#' @return Histogram values
#' @export
#' @seealso \code{\link{npbin}}, \code{\link{emBspl}}, \code{\link{ebBeta}}
initialize_weights <- function(data_table, n_breaks, spline_order, plot = FALSE) {
  h <- hist(
    data_table[, p_hat],
    breaks = seq(0, 1, length.out = n_breaks + spline_order - 3),
    plot = plot
  )
  h[["density"]]
}
