#===============================================================================
# symmetric.R
#===============================================================================

#' @title Calculate beta variance from shape parameters
#'
#' @param shape1 Numeric. First shape parameter.
#' @param shape2 Numeric. Second shape parameter.
#' @return Numeric. The variance
#' @export
beta_variance <- function(shape1, shape2) {
  shape1 * shape2 / ((shape1 + shape2)^2 * (shape1 + shape2 + 1))
}

#' @title Calculate beta overdispersion from shape parameters
#'
#' @param shape1 Numeric. First shape parameter.
#' @param shape2 Numeric. Second shape parameter.
#' @return Numeric. The overdispersion
#' @export
overdispersion <- function(shape1, shape2) {
  1 / (shape1 + shape2 + 1)
}


#' @title Calculate a symmetric shape parameter from the variance
#'
#' @param variance Numeric. The variance
#' @return Numeric. The shape parameter
#' @export
symmetric_beta_shape_from_variance <- function(variance) {
  1/(8*variance) - 1/2
}