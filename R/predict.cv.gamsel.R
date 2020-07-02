#' Predict from \code{cv.gamsel} object
#' 
#' @param object \code{cv.gamsel} object
#' @param index Integer: Which lambda to use. If NULL, choose lambda based on \code{which.lambda}.
#' If specified, \code{which.lambda} is not used. Default = NULL
#' @param which.lambda String: "lambda.min" or "lambda.1se". Default = "lambda.min". This will not 
#' be used if \code{index} is specified
#' @param type String: "link", "response", "terms", "nonzero"
#' @author Efstathios D. Gennatas
#' @export

predict.cv.gamsel <- function(object, newdata,
                              index = NULL,
                              which.lambda = "lambda.min",
                              type = c("link", "response", "terms", "nonzero"), ...) {
  
  if (is.null(index)) {
    index <- if (which.lambda == "lambda.min") object$index.min else object$index.1se
  }
  predict(object$gamsel.fit, newdata = newdata,
          index = index,
          type = type, ...)
  
}