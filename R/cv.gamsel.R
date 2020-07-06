cv.gamsel <- function(x, y,
                      lambda = NULL,
                      family = c("gaussian", "binomial"), 
                      # degrees = rep(10, p),
                      degrees = NULL,
                      min.degree = 1,
                      max.degree = 8,
                      # dfs = rep(5, p),
                      failsafe = TRUE,
                      dfs = NULL,
                      min.df = 1,
                      max.df = 5,
                      bases = pseudo.bases(x, degrees, dfs, parallel = parallel, ...),
                      type.measure = c("mse", "mae", "deviance", "class"),
                      nfolds = 10, foldid, keep = FALSE, parallel = FALSE, 
                      trace = 0,
                      ...) {
  
  family <- match.arg(family)
  p <- NCOL(x)
  y <- drop(y)
  if (missing(type.measure)) 
    type.measure <- switch(family, gaussian = "mse", binomial = "deviance")
  else type.measure <- match.arg(type.measure)
  if (any(match(type.measure, c("deviance", "class"), FALSE)) & 
      family == "gaussian") {
    warning(paste(type.measure, "not available for gaussian family; will use mse instead"))
    type.measure <- "mse"
  }
  typenames <- c(mse = "Mean-Squared Error", mae = "Mean Absolute Error", 
                deviance = "Binomial Deviance", class = "Misclassification Error")
  if (!is.null(lambda) && length(lambda) < 2) 
    stop("Need more than one value of lambda for cv.gamsel")
  nobs <- nrow(x)
  
  # delta
  n.features <- NCOL(x)
  
  unique_perfeat <- apply(x, 2, function(i) length(unique(i)))
  
  if (any(unique_perfeat < 4) && failsafe) {
    # Consider having gamsel output gamsel object even when failsafe trigerred
    warning("Feature with less than 4 unique values found")
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
  
  # y = drop(y)
  gamsel.object <- gamsel(x, y,
                          lambda = lambda,
                          bases = bases,
                          family = family,
                          failsafe = failsafe,
                         ...)
  lambda <- gamsel.object$lambda
  if (missing(foldid)) 
    foldid <- sample(rep(seq(nfolds), length = nobs))
  else nfolds <- max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist <- as.list(seq(nfolds))
  if (parallel) {
    #  if (parallel && require(foreach)) {
    outlist <- foreach(i = seq(nfolds), .packages = c("gamsel")) %dopar% 
      {
        which <- foldid == i
        bases.sub <- lapply(bases, basis.subset, subset = !which)
        gamsel(x[!which, , drop = FALSE], y[!which], 
               lambda = lambda, bases = bases.sub, family = family,...)
      }
  } else {
    for (i in seq(nfolds)) {
      which <- foldid == i
      bases.sub <- lapply(bases, basis.subset, subset = !which)
      outlist[[i]] = gamsel(x[!which, , drop = FALSE], 
                            y[!which], lambda = lambda, bases = bases.sub, family = family,
                            ...)
    }
  }
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq_len(nfolds)) {
    which <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[which, , drop = FALSE], type = "response")
    nlami <- length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] <- preds
    nlams[i] = nlami
  }
  prob_min = 1e-04
  prob_max = 1 - prob_min
  N = nobs - apply(is.na(predmat), 2, sum)
  cvraw = switch(type.measure,
                 mse = (y - predmat)^2,
                 mae = abs(y - predmat),
                 deviance = {
                   predmat = pmin(pmax(predmat, prob_min), prob_max)
                   lp = y * log(predmat) + (1 - y) * log(1 - predmat) - 2 * lp
                 },
                 class = y * (predmat < 0.5) + (1 - y) * (predmat >= 0.5))
  nz = sapply(predict(gamsel.object, type = "nonzero"), length)
  cvm = apply(cvraw, 2, mean, na.rm = TRUE, trim = 0.01)
  cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE, trim = 0.01)/(N - 
                                                                                       1))
  out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + 
               cvsd, cvlo = cvm - cvsd, nzero = nz, name = typenames[type.measure], 
             gamsel.fit = gamsel.object)
  if (keep) 
    out = c(out, list(fit.preval = predmat, foldid = foldid))
  lamin = getmin(lambda, cvm, cvsd)
  imin = match(lamin, lambda)
  names(imin) = c("index.min", "index.1se")
  obj <- c(out, as.list(c(lamin, imin)))
  class(obj) <- "cv.gamsel"
  obj
}
