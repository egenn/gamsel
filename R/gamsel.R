#' Fit Regularization Path for Gaussian or Binomial Generalized Additive Model
#'
#' Using overlap grouped lasso penalties, gamsel selects whether a term in a gam is nonzero, 
#' linear, or a non-linear spline (up to a specified max df per variable). It fits the entire 
#' regularization path on a grid of values for the overall penalty lambda, both for gaussian and 
#' binomial families.
#' 
#' The sequence of models along the lambda path is fit by (block) cordinate descent. In the case of 
#' logistic regression the fitting routine may terminate before all num_lambda values of lambda have 
#' been used. This occurs when the fraction of null deviance explained by the model gets too close 
#' to 1, at which point the fit becomes numerically unstable. Each of the smooth terms is computed 
#' using an approximation to the Demmler-Reinsch smoothing spline basis for that variable, and the 
#' accompanying diagonal pernalty matrix.
#' 
#' @param x Input (predictor) matrix of dimension nobs x nvars. Each observation is a row.
#' @param y Response variable. Quantitative for family="gaussian" and with values in {0,1} for 
#' family="binomial"
#' @param num_lambda Number of lambda values to use. (Length of lambda sequence.)
#' @param lambda User-supplied lambda sequence. For best performance, leave as NULL and allow the 
#' routine to automatically select lambda. Otherwise, supply a (preferably gradually) decreasing 
#' sequence.
#' @param family Response type. "gaussian" for linear model (default). "binomial" for logistic 
#' model.
#' @param degrees An integer vector of length nvars specifying the maximum number of spline basis 
#' functions to use for each variable.
#' @param min.degree 
#' @param max.degree 
#' @param gamma Penalty mixing parameter \eqn{0 \le\gamma\le 1}.  Values \eqn{ \gamma <37   0.5}  
#' penalize linear fit less than non-linear fit. The default is \eqn{\gamma = 0.4}, which encourages 
#' a linear term over a nonlinear term.
#' @param dfs Numeric vector of length nvars specifying the maximum (end-of-path) degrees of 
#' freedom for each variable.
#' @param min.df A list of orthonormal bases for the non-linear terms for each variable. The 
#' function \link{pseudo.bases} generates these, using the parameters \code{dfs} and \code{degrees}. 
#' See the documentation for \link{pseudo.bases}.
#' @param max.df 
#' @param bases 
#' @param tol Convergence threshold for coordinate descent. The coordinate descent loop continues 
#' until the total change in objective after a pass over all variables is less than tol. 
#' Default is \code{1e-4}.
#' @param max_iter Maximum number of coordinate descent iterations over all the variables for each 
#' lambda value. Default is 2000.
#' @param traceit If \code{TRUE}, various information is printed during the fitting process.
#' @param parallel passed on to the pseudo.bases() function. Uses multiple process if available.
#' @param failsafe 
#' @param trace 
#' @param ... additional arguments passed on to \code{pseudo.bases()}
#' @return
#' An object with S3 class \code{gamsel}.
#' @author Alexandra Chouldechova and Trevor Hastie
#' @export
#' @examples
#' \dontrun{
#' data <- gendata(n = 500, p = 12, k.lin = 3, k.nonlin = 3, de g =8, sigma = 0.5)
#' attach(data)
#' bases <- pseudo.bases(X, degree = 10, df = 6)
#' # Gaussian gam
#' gamsel.out <- gamsel(X, y, bases = bases)
#' par(mfrow=c(1,2),mar=c(5,4,3,1))
#' summary(gamsel.out)
#' gamsel.cv=cv.gamsel(X,y,bases=bases)
#' par(mfrow=c(1,1))
#' plot(gamsel.cv)
#' par(mfrow=c(3,4))
#' plot(gamsel.out,newx=X,index=20)
#' # Binomial model
#' gamsel.out=gamsel(X,yb,family="binomial")
#' par(mfrow=c(1,2),mar=c(5,4,3,1))
#' summary(gamsel.out)
#' par(mfrow=c(3,4))
#' plot(gamsel.out,newx=X,index=30)
#' }

gamsel <- function(x, y,
                   num_lambda = 50,
                   lambda = NULL,
                   family = c("gaussian", "binomial"),
                   # degrees = rep(10, p),
                   degrees = NULL,
                   min.degree = 1,
                   max.degree = 8,
                   gamma = 0.4,
                   # dfs = rep(5, p),
                   dfs = NULL,
                   min.df = 1,
                   max.df = 5,
                   bases = pseudo.bases(x, degrees, dfs, parallel = parallel, ...),
                   tol = 1e-04,
                   max_iter = 2000, traceit = FALSE, parallel = FALSE,
                   failsafe = TRUE,
                   failsafe.args = list(learner = "s.GLMNET", 
                                        family = family,
                                        alpha = 1,
                                        lambda = NULL),
                   trace = 0, ...) {
  
  # delta
  n.features <- NCOL(x)
  unique_perfeat <- apply(x, 2, function(i) length(unique(i)))
  
  if (any(unique_perfeat < 4) && failsafe) {
    warning("Failsafe on: Feature with less than 4 unique values found, returning glmnet")
    learner <- failsafe.args[[1]]
    failsafe.args[[1]] <- NULL
    .failsafe.args <- c(list(x = x, y = y),
                       failsafe.args
                       )
    return(do.call(learner, .failsafe.args))
  }
  
  if (trace > 1) cat(".: Unique vals per feat:", unique_perfeat, "\n")
  
  if (is.null(degrees)) {
    degrees <- sapply(seq_len(n.features), function(i) 
      max(min.degree, min(unique_perfeat[i] - 1, max.degree)))
  }
  if (length(degrees) < n.features) degrees <- rep(degrees, n.features)[seq_len(n.features)]
  if (trace > 1) cat(".: 'degrees' set to:", degrees, "\n")
  
  if (is.null(dfs)) {
    # -2 is playing it safe to test: fix
    dfs <- sapply(seq_len(n.features), function(i) max(min.df, min(degrees[i] - 2, max.df)))
  }
  if (length(dfs) != n.features) dfs <- rep(dfs, n.features)[seq_len(n.features)]
  if (trace > 1) cat(".: 'dfs' set to:", dfs, "\n")
  # /
  
  # this.call <- match.call()
  family <- match.arg(family)
  family_id <- ifelse(family == "gaussian", 0L, 1L)
  # n <- length(y)
  p <- NCOL(x)
  ## Create U, X, D and psis from the bases
  degrees <- sapply(bases, dim)[2, ]
  U <- do.call("cbind", bases)
  X <- do.call("cbind", lapply(bases, function(x) x[, 1, drop = FALSE]))
  parms <- lapply(bases,"attr", "parms")
  getdpsi  <- function(parms) {
    d <- parms$d
    if (length(d) > 1) {
      psi <- d[2]
      d <- d/d[2]
      d[1] <- 1
    } else {
      d <- 1
      psi <- 0
    }
    list(d = d, psi = psi)
  }
  dpsi <- lapply(parms, getdpsi)
  D_mat <- unlist(lapply(dpsi,"[[","d"))
  psi <- sapply(dpsi, "[[", "psi")
  if (is.null(lambda)) lambda <- rep(-0.5, num_lambda) else {
    lambda <- as.numeric(lambda)
    lambda <- rev(unique(sort(lambda)))
  }
  num_lambda <- as.integer(length(lambda))
  degrees <- as.integer(degrees)
  max_iter <- as.integer(max_iter)
  out <- .Call("gamselFit", y, X, U, tol, degrees, D_mat, gamma, 
               psi, family_id, max_iter, lambda, num_lambda,
               as.integer(traceit))[c("intercept", "alphas", "betas", "lambdas")]
  out$degrees <- degrees
  out$parms <- parms
  out$family <- family
  out <- c(out, fracdev(U, y, out$alphas, out$betas, out$intercept, degrees, family))
  # out$call <- this.call
  class(out) <- "gamsel"
  out
}
