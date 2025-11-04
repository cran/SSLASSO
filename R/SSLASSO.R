#' The Spike-and-Slab LASSO
#'
#' @description
#' Spike-and-Slab LASSO is a spike-and-slab refinement of the LASSO procedure,
#' using a mixture of Laplace priors indexed by `lambda0` (spike) and `lambda1` (slab).
#'
#' The `SSLASSO` procedure fits coefficients paths for Spike-and-Slab LASSO-penalized
#' linear regression models over a grid of values for the regularization
#' parameter `lambda0`. The code has been adapted from the `ncvreg` package (Breheny and Huang, 2011).
#'
#' @param X The design matrix (n x p), without an intercept. `SSLASSO` standardizes the data by default.
#' @param y Vector of continuous responses (n x 1). The responses will be centered by default.
#' @param penalty The penalty to be applied to the model. Either "separable" (with a fixed `theta`) or "adaptive" (with a random `theta`, where `theta ~ B(a,p)`).
#' @param variance Specifies whether the error variance is "fixed" or "unknown".
#' @param lambda1 The slab penalty parameter. Must be smaller than the smallest `lambda0`.
#' @param lambda0 A sequence of spike penalty parameters. Must be monotone increasing. If not specified, a default sequence is generated.
#' @param beta.init Initial values for the coefficients. Defaults to a vector of zeros.
#' @param nlambda The number of `lambda0` values to use in the default sequence. Defaults to 100.
#' @param theta The initial mixing proportion for the spike component. Defaults to 0.5.
#' @param sigma The initial value for the error standard deviation. Defaults to 1.
#' @param a Hyperparameter for the Beta prior on theta. Defaults to 1.
#' @param b Hyperparameter for the Beta prior on theta. Defaults to the number of predictors, p.
#' @param eps Convergence tolerance. The algorithm stops when the maximum change in coefficients is less than `eps`. Defaults to 0.001.
#' @param max.iter The maximum number of iterations. Defaults to 500.
#' @param counter The number of iterations between updates of the adaptive penalty parameters. Defaults to 10.
#' @param warn A logical value indicating whether to issue a warning if the algorithm fails to converge. Defaults to `FALSE`.
#'
#' @return An object with S3 class "SSLASSO". The object contains:
#' \item{beta}{A p x L matrix of estimated coefficients, where L is the number of regularization parameter values.}
#' \item{intercept}{A vector of length L containing the intercept terms.}
#' \item{iter}{The number of iterations for each value of `lambda0`.}
#' \item{lambda0}{The sequence of `lambda0` values used.}
#' \item{lambda1}{The `lambda1` value used.}
#' \item{penalty}{The penalty type used.}
#' \item{thetas}{A vector of length L containing the hyper-parameter values `theta` (the same as `theta` for "separable" penalty).}
#' \item{sigmas}{A vector of length L containing the values `sigma` (the same as the initial `sigma` for "known" variance).}
#' \item{select}{A (p x L) binary matrix indicating which variables were selected along the solution path.}
#' \item{model}{A single model chosen after the stabilization of the regularization path.}
#' \item{n}{The number of observations.}
#'
#' @author Veronika Rockova <Veronika.Rockova@chicagobooth.edu>, Gemma Moran <gm845@stat.rutgers.edu>
#'
#' @references
#' Ročková, V., & George, E. I. (2018). The spike-and-slab lasso. Journal of the American Statistical Association, 113(521), 431-444.
#'
#' Moran, G. E., Ročková, V., & George, E. I. (2019). Variance prior forms for high-dimensional bayesian variable selection. Bayesian Analysis, 14(4), 1091-1119.
#'
#' @seealso [plot.SSLASSO()]
#'
#' @importFrom stats sd qchisq model.matrix 
#' @useDynLib SSLASSO, .registration = TRUE
#' @export
#' @examples
#' ## Linear regression, where p > n
#' 
#' library(SSLASSO)
#' 
#' p <- 100
#' n <- 50
#' 
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' beta <- c(1, 2, 3, rep(0, p-3))
#' y = X[,1] * beta[1] + X[,2] * beta[2] + X[,3] * beta[3] + rnorm(n)
#' 
#' # Oracle SSLASSO with known variance
#' 
#' result1 <- SSLASSO(X, y, penalty = "separable", theta = 3/p)
#' plot(result1)
#' 
#' # Adaptive SSLASSO with known variance
#' 
#' result2 <- SSLASSO(X, y)
#' plot(result2)
#' 
#' # Adaptive SSLASSO with unknown variance
#' 
#' result3 <- SSLASSO(X, y, variance = "unknown")
#' plot(result3)
#' 
#' 
#'
SSLASSO <- function(X,
                    y,
                    penalty = c("adaptive", "separable"),
                    variance = c("fixed", "unknown"),
                    lambda1,
                    lambda0,
                    beta.init = numeric(ncol(X)),
                    nlambda = 100,
                    theta = 0.5,
                    sigma = 1,
                    a = 1, b,
                    eps = 0.001,
                    max.iter = 500,
                    counter = 10,
                    warn = FALSE) {
  
  # Coersion
  penalty <- match.arg(penalty)
  variance <- match.arg(variance)
  
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) {
      stop("X must be a matrix or able to be coerced to a matrix")
    }
  }
  if (storage.mode(X) == "integer") {
    storage.mode(X) <- "double"
  }
  if (!is.numeric(y)) {
    tmp <- try(y <- as.numeric(y), silent=TRUE)
    if (inherits(tmp,"try-error")) {
      stop("y must numeric or able to be coerced to numeric")
    }
  }
  
  if (any(is.na(y)) | any(is.na(X))) {
    stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to SSLASSO")
  }
  
  ## Standardize
  XX <- scale(X, center = TRUE, scale = FALSE)
  # ns <- attr(XX, "nonsingular")
  p <- ncol(XX)
  
  yy <- y - mean(y)
  n <- length(yy)
  
  if (missing(lambda0)) {
    lambda0 <- seq(1, n, length = nlambda)
    lambda1 <- lambda0[1]
  } else {
    nlambda <- length(lambda0)
    if (missing(lambda1)) {
      lambda1 <- lambda0[1]
    }
  }
  
  # Lambda0 should be an increasing sequence
  
  monotone <- sum((lambda0[-1] - lambda0[-nlambda]) > 0)
  if (monotone != nlambda - 1){
    stop("lambda0 must be a monotone increasing sequence")
  }
  if (lambda1 > min(lambda0) ) {
    stop("lambda1 must be smaller than lambda0")
  }
  
  if(missing(b)) {
    b <- p
  }
  
  # get initial value for sigma
  df = 3
  sigquant = 0.9
  sigest <- sd(yy)
  qchi <- qchisq(1 - sigquant, df)
  ncp <- sigest^2 * qchi / df
  min_sigma2 <- sigest^2 / n
  
  if (variance == "unknown") {
    if (missing(sigma)) {
      sigma <- sqrt(df * ncp / (df + 2))
    }
  }
  
  ## Fit
  res <- .Call("SSL_gaussian", XX, yy, beta.init, penalty, variance, as.double(lambda1), as.numeric(lambda0),
               as.double(theta), as.double(sigma), as.double(min_sigma2), as.double(a), as.double(b),
               eps, as.integer(max.iter), as.integer(counter), PACKAGE = "SSLASSO")
  bb <- matrix(res[[1]], p, nlambda)
  iter <- res[[3]]
  thetas<-res[[4]]
  sigmas <- res[[5]]
  
  ## Warning
  if (warn & any(iter == max.iter)) {
    warning("Algorithm did not converge for the ABOVE MENTIONED values of lambda0: ", paste(lambda0[iter == max.iter], collapse = ", "))
  }
  
  if (iter[nlambda] == max.iter) {
    warning("Algorithm did not converge at the last lambda0 value.")
  }
  
  ## Unstandardize
  beta <- matrix(0, nrow = ncol(X), ncol = nlambda)
  # bbb <- bb/attr(XX, "scale")[ns]
  # beta[ns, ] <- bbb
  bbb <- bb
  beta <- bb
  
  intercept <- rep(mean(y), nlambda) - crossprod(attr(XX, "scaled:center"), bbb)
  
  ## Names
  varnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep = "") else colnames(X)
  varnames <- c(varnames)
  dimnames(beta) <- list(varnames, round(lambda0,digits=4))
  
  ## Select
  select <- apply(beta, 2, function(x){as.numeric(x!=0)})
  
  ## Model
  
  model<-(1:p)[select[,nlambda]==1]
  
  ## Output
  val <- structure(list(beta = beta,
                        intercept = intercept,
                        iter = iter,
                        lambda0 = lambda0,
                        penalty = penalty,
                        lambda1 = lambda1,
                        thetas = thetas,
                        sigmas = sigmas,
                        select = select,
                        model = model,
                        n = n),
                   class = "SSLASSO")
  
  val
}

