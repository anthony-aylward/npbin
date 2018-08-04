#===============================================================================
# plot_estimated_null.R
#===============================================================================

#' @title Plot an estimated null distribution
#'
#' @description optionally superimpose two null distributions
#'
#' @param p_hat Numeric. Vector of observed ref allele proportions.
#' @param shape1_shape2 Numeric. Shape parameters of the first null
#'   distribution.
#' @param shape3_shape4 Numeric. Shape parameters of the second null
#'   distribution.
#' @export
plot_estimated_null <- function(p_hat, shape1_shape2, shape3_shape4 = NULL) {
  pal <- color_palette()
  null_x <- 1:128 / 129
  shape1 = shape1_shape2[[1]]
  shape2 = shape1_shape2[[2]]
  null_y1 <- dbeta(null_x, shape1, shape2)
  if (!is.null(shape3_shape4)) {
    shape3 = shape3_shape4[[1]]
    shape4 = shape3_shape4[[2]]
    null_y2 <- dbeta(null_x, shape3, shape4)
  } else {
    null_y2 <- c()
  }
  hist(
    p_hat,
    plot = TRUE,
    freq = FALSE,
    col = "lavenderblush",
    border = pal[["pink"]],
    lty = 2,
    ylim = c(0, max(c(null_y1, null_y2))),
    xlab = "Fraction Ref Allele",
    main = ""
  )
  par(new = TRUE)
  plot(
    null_x,
    null_y1,
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
  if (!is.null(shape3_shape4)) {
    par(new = TRUE)
    plot(
      null_x,
      null_y1,
      type = "l",
      lty = 2,
      lwd = 2,
      col = pal[["teal"]],
      bty = "n",
      ylab = "",
      xlab = "",
      yaxt = "n",
      xaxt = "n"
    )
  }
}