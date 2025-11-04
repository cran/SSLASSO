# code has been adapted from the ncvreg package (Breheny and Huang, 2011)

#' Standardizes a design matrix
#'
#' @description The function `standard` accepts a design matrix and returns
#' a standardized version of that matrix (i.e., each column will have
#' mean 0 and mean sum of squares equal to 1). The code has been adapted
#' from the `ncvreg` package (Breheny and Huang, 2011).
#'
#' @param X A matrix (or object that can be coerced to a matrix, such as a data frame).
#'
#' @details This function centers and scales each column of `X` so that
#' \deqn{\sum_{i=1}^n x_{ij}=0}{sum(X[,j])=0}
#' and
#' \deqn{\sum_{i=1}^n x_{ij}^2 = n}{sum(X[,j]^2)=n}
#' for all j. This is usually not necessary to call directly, as `SSLASSO`
#' internally standardizes the design matrix, but inspection of the
#' standardized design matrix can sometimes be useful. This differs from
#' the base R function `scale` in two ways: (1) `scale`
#' uses the sample standard deviation `sqrt(sum(x^2)/(n-1))`, while
#' `standard` uses the root-mean-square, or population, standard
#' deviation `sqrt(mean(sum(x^2)))`, and (2) `standard` is faster.
#' The reason for using the population standard deviation is that `SSLASSO`
#' assumes that the columns of the design matrix have been scaled to have
#' norm `sqrt(n)`.
#'
#' @return A standardized matrix with attributes for center, scale, and non-singular columns.
#'
#' @importFrom stats model.matrix
#'
standard <- function(X) {
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix")
  }
  STD <- .Call("standardize", X, PACKAGE = "SSLASSO")
  dimnames(STD[[1]]) <- dimnames(X)
  ns <- which(STD[[3]] > 1e-6)
  if (length(ns) == ncol(X)) {
    val <- STD[[1]]
  } else {
    val <- STD[[1]][, ns, drop=FALSE]
  }
  attr(val, "center") <- STD[[2]]
  attr(val, "scale") <- STD[[3]]
  attr(val, "nonsingular") <- ns
  val
}