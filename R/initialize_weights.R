#===============================================================================
# initialize_weights.R
#===============================================================================

#' @Initialize weights
#'
#' @details
#' Initialize the weights using the histogram of p_hat
#'
#' @param data_table Input data table
#' @param n_breaks Number of breaks for the histogram
#' @param spline_order Order of splines
#' @return Histogram values
#' @export
initialize_weights <- function(data_table, n_breaks, spline_order) {
  hist(
    data_table[, p_hat],
    breaks = seq(0, 1, length.out = n_breaks + spline_order - 3),
    plot = FALSE
  )[["density"]]
}
