gamsel <- function(x, y,
                   num_lambda = 50,
                   lambda = NULL,
                   family = c("gaussian","binomial"),
                   # degrees = rep(10, p),
                   degrees = NULL,
                   min.degree = 1,
                   max.degree = 8,
                   gamma = 0.4,
                   # dfs = rep(5, p),
                   dfs = NULL,
                   min.df = 1,
                   max.df = 5,
                   bases = pseudo.bases(x,degrees, dfs, parallel = parallel, ...),
                   tol = 1e-04,
                   max_iter = 2000, traceit = FALSE, parallel = FALSE,
                   failsafe = TRUE,
                   trace = 0, ...) {
  
  # egenn
  n.features <- NCOL(x)
  
  unique_perfeat <- apply(x, 2, function(i) length(unique(i)))
  
  if (any(unique_perfeat < 4) && failsafe) {
    warning("Feature with less than 4 unique values found, returning glm")
    return(glm(y ~ data.matrix(x)))
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
  family_id <- as.integer(ifelse(family == "gaussian", 0, 1))
  # n <- length(y)
  p <- ncol(x)
  ## Create U, X, D and psis  from the bases
  degrees <- sapply(bases, dim)[2,]
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
  dpsi <- lapply(parms,getdpsi)
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
