#===============================================================================
# plot_estimated_null.R
#===============================================================================

#' @title Plot an estimated null distribution
#'
#' @param p_hat Numeric. Vector of observed ref allele proportions.
#' @param shape1 Numeric. First shape parameter of the null distribution.
#' @param shape2 Numeric. Second shape parameter of the null distribution.
#' @export
plot_estimated_null <- function(p_hat, shape1, shape2) {
  pal <- color_palette()
  null_x <- 1:128 / 129
  null_y <- dbeta(null_x, shape1, shape2)
  hist(
    p_hat,
    plot = TRUE,
    freq = FALSE,
    col = "lavenderblush",
    border = pal[["pink"]],
    lty = 2,
    ylim = c(0, max(null_y)),
    xlab = "Fraction Ref Allele",
    main = ""
  )
  par(new = TRUE)
  plot(
    null_x,
    null_y,
    type = "l",
    lty = 2,
    lwd = 2,
    col = pal[["blue"]],
    bty = "n",
    ylab = "",
    xlab = "",
    yaxt = "n",
    xaxt = "n"
  )
}