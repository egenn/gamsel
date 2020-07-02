predict.gamsel <- function(object, newdata,
                           index = NULL,
                           type = c("link","response","terms","nonzero"), ...) {
  type <- match.arg(type)
  lambda <- object$lambda
  nlambda <- length(lambda)
  if (is.null(index)) index <- seq_len(nlambda) else index <- intersect(index, seq(nlambda))
  if (type == "nonzero") return(getActive(object, index, type = "nonzero"))
  if (missing(newdata)) stop("newdata is required for prediction")
  parms <- object$parms
  betas <- object$betas
  # dimb <- dim(betas)
  degrees <- object$degrees
  p <- length(degrees)
  dimx <- dim(newdata)
  if (dimx[2] != p) stop(paste("number of columns of x different from", p))
  U <- as.list(seq_len(p))
  offset <- c(1, 1 + cumsum(degrees)[-p])
  for (i in seq_len(p)) U[[i]] <- basis.gen(newdata[,i], degree = degrees[i], parms = parms[[i]])
  betas <- betas[, index, drop = FALSE]
  betas[offset, ] <- betas[offset, ] + object$alpha[, index]
  dimnames(betas) <- list(NULL, paste0("l", index))
    pred <- switch(type,
                   terms = { 
                     fitlist = as.list(seq_len(p))
                 which = rep(seq_len(p), degrees)
                 for (i in seq_len(p)) {
                    beta = betas[which == i, , drop = FALSE]
                    fitlist[[i]] = U[[i]] %*% beta
                  }
                 dd = c(dim(fitlist[[1]]), length(fitlist))
                 dn = c(dimnames(fitlist[[1]]), list(paste0("v", seq_len(p))))
                 array(do.call("cbind", fitlist), dd, dn)
              },
           {U = do.call("cbind", U)
            
              U %*% betas + rep(1, nrow(U)) %o% object$intercept[index]
            }
           )
  if (type == "response" && object$family == "binomial") {
    pred <- exp(pred)
    pred <- pred/(1 + pred)
  }
  pred
}
