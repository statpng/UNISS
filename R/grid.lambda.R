#' @name grid.lambda
#' @title Finding a Grid of Lambda
#' @description
#' Finding a grid of lambda covering all possible lambda values in the provided dataset \code{(x, y)}.
#'
#' @param x An \code{n} by \code{p} matrix of variants; Each row represents an individual, and each column a genetic variant. Can be high-dimensional data with \code{n<p}.
#' @param y An \code{n} by \code{q} matrix of phenotypes, each with being either continuous or categorical.
#' @param family A vector of distribution family, such as 'gaussian', 'binomial', 'multinomial', and etc, corresponding to phenotypes.
#' @param iter The number of iterations to be performed for finding a lambda grid. Default is 10.
#' @param alpha Mixing proportion between ridge and lasso penalties used in elastic-net penalty. Default is 0.1.
#' @param nlambda The number of lambda values to be considered at each iteration. Default is 10.
#' @param psub Subsampling proportion. When \code{psub=1.0} bootstrap resampling is applied, and otherwise subsampling without replacement is applied. Default is 1.0.
#' @param pf A vector of penalty to be imposed for each column of \code{x}. The 0 value indicates no penalty which leads the genetic variants always to be included in the model. On the other hand, the 1 value indicates that the penalty is imposed for each column as it was.
#' @param seed An integer value for reproducible outcome. Default is 1 that does not have special meaning.
#' @param ... Additional arguments passing to function \code{glmnet::glmnet}.
#'
#' @usage
#' grid.lambda(x, y, family, iter = 10, alpha = 0.1,
#'             nlambda = 10, psub = 1.0, pf = NULL, seed = 1, ...)
#'
#'
#'
#'
#' @return
#' A list of lambda sequences. Each sequence contains \code{nlambda*iter} lambda values.
#'
#'
#' @import glmnet
#'
#' @export grid.lambda
grid.lambda <- function(x, y,
                       family,
                       iter = 10,
                       alpha = 0.1,
                       nlambda = 10,
                       psub = 1.0,
                       pf = NULL,
                       seed = 1,
                       ...) {

  Call <- match.call()

  if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
  if( is.null(pf)) pf <- rep(1, ncol(x))


  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  nsub <- n * psub

  sample.replace = FALSE
  if( psub == 1.0 ) sample.replace = TRUE

  q <- ncol(y)


  lambda.array <- array( NA, dim = c(iter, nlambda, length(alpha), ncol(y), 2),
                         dimnames = list(
                           paste0("iter.", 1:iter),
                           paste0("lambda.", 1:nlambda),
                           paste0("alpha=", alpha),
                           paste0("y.", 1:ncol(y)),
                           c("lambda", "df")
                         )
  )


  seq.lambda <- as.list(1:q)

  set.seed(seed)
  for (k in 1:q) {
    FAMILY <- family[k]

    lambda.vec <- NULL
    for (ii in 1:iter) {
      for (aa in 1:length(alpha)) {
        wsub <- sample(n, nsub)
        xsub <- x[wsub, , drop = F]
        ysub <- y[wsub, , drop = F]
        if (FAMILY == "multinomial") {
          fitsub <- glmnet(
            x = xsub,
            y = y.multinom(ysub[,k,drop=T]),
            alpha = alpha[aa],
            family = FAMILY,
            nlambda = nlambda,
            type.multinomial = "grouped",
            penalty.factor = pf,
            ... )
        } else {
          fitsub <- glmnet(
            x = xsub,
            y = ysub[, k],
            alpha = alpha[aa],
            family = FAMILY,
            nlambda = nlambda,
            penalty.factor = pf,
            ...
          )
        }


        lambda.array[ii, 1:length(fitsub$df), aa, k, 1] <- fitsub$lambda
        lambda.array[ii, 1:length(fitsub$df), aa, k, 2] <- fitsub$df

        lambda.vec <- c(lambda.vec, fitsub$lambda)
      }
    }

    seq.lambda[[k]] <- lambda.vec

  }


  attr(seq.lambda, "call") <- Call

  return(seq.lambda)

}




#' @importFrom stats model.matrix
y.multinom <- function(yk){
  droplevels( model.matrix(~ . - 1, data=as.data.frame(as.factor(yk) )))
}

