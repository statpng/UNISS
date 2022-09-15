#' @name uniss
#' @title Unified Selection Score for Multi-Trait Analysis
#' @description
#' New statistical selection method for pleiotropic variants associated with both quantitative and qualitative traits, based on a penalized regression and a resampling approach. Details are refer to the paper in reference section.
#'
#' @param x An \code{n} by \code{p} matrix of variants; Each row represents an individual, and each column a genetic variant. Can be high-dimensional data with \code{n<p}.
#' @param y An \code{n} by \code{q} matrix of phenotypes, each with being either continuous or categorical.
#' @param family A vector of distribution family, such as 'gaussian', 'binomial', 'multinomial', and etc, corresponding to phenotypes.
#' @param alpha Mixing proportion between ridge and lasso penalties used in elastic-net penalty. Default is 0.1.
#' @param lambda A grid vector of the penalty parameter used in elastic-net penalty. It can be provided through our function 'grid.lambda'.
#' @param B The number of resampling replicates when calculating the selection probability. Default is 100.
#' @param psub Subsampling proportion. When \code{psub=1.0} bootstrap resampling is applied, and otherwise subsampling without replacement is applied. Default is 1.0.
#' @param pf A vector of penalty to be imposed for each column of \code{x}. The 0 value indicates no penalty which leads the genetic variants always to be included in the model. On the other hand, the 1 value indicates that the penalty is imposed for each column as it was.
#' @param std.x Logical flag for x variable standardization. If TRUE, each variant is scaled to have zero mean and unit variance. Default is TRUE.
#' @param std.y Logical flag for y variable standardization. If TRUE, each phenotype is scaled to have zero mean and unit variance. Default is TRUE.
#' @param verbose Logical indicating whether or not to print out progress indicators.
#' @param ... Additional arguments passing to function \code{glmnet::glmnet}.
#'
#'
#'
#' @references
#' Kim, K., Jun, T., Wang, S., & Sun, H. (2022). New statistical selection method for pleiotropic variants associated with both quantitative and qualitative traits. \emph{Bioinformatics}, submitted.
#'
#'
#' @return
#' A list of fitting results containing
#' \item{selscore}{a matrix of which the first four columns indicate the individual selection probabilities for phenotypes and the last column the unified selection score calculated from the adjusted selection probabilities.}
#' \item{params}{arguments used in model fitting}
#'
#'
#'
#' @examples
#'
#' set.seed(1)
#' n=100; p=200; q=4
#'
#' X <- replicate(p, rbinom(n,2,0.2)) #generate the SNP X=0,1,2
#' b <- matrix(0, p, q)
#' b[1:5,1:4] <- 1.0
#' b[6:10,1:2] <- 1.0
#'
#'
#' Z <- replicate(1, rnorm(n))
#' g <- matrix(0, 1, q)
#' g[1,1:4] <- 0.1
#'
#' x <- cbind(Z, X)
#' beta <- rbind(g, b)
#'
#' y <- x%*%beta + replicate(q, rnorm(n,1))
#' y[,2:3] <- apply(y[,2:3], 2, function(yk) ifelse(yk > median(yk), 1, 0) )
#' family <- c("gaussian","binomial","binomial","gaussian")
#' alpha <- 0.1
#' pf <- c(0, rep(1, p))
#'
#' lambda.all <- grid.lambda(x = x, y = y,
#'                                family = family,
#'                                iter = 10,
#'                                seq.alpha = alpha, nlambda = 10,
#'                                pf=pf)
#' lambda <- lapply(lambda.all, function(x) seq(median(x), max(x), length.out=10) )
#'
#'
#' fit.uniss <- uniss( x = x, y = y, family = family,
#'                     psub = 1.0,
#'                     alpha = 0.1,
#'                     lambda = lambda,
#'                     B = 20, # Default is 100
#'                     pf = pf,
#'                     std.x = TRUE,
#'                     std.y = TRUE )
#'
#' fit.cutoff <- cutoff(params = fit.uniss$params, nperm=10)
#'
#' head(fit.cutoff, 10)
#'
#'
#'
#'
#' @import glmnet
#'
#' @export uniss
uniss <- function(x,
                 y,
                 family,
                 lambda = NULL,
                 alpha = 0.1,
                 B = 100,
                 psub = 1.0,
                 pf = NULL,
                 std.x = TRUE,
                 std.y = TRUE,
                 verbose = TRUE,
                 ...) {


  Call <- match.call()


  if( length(family) != ncol(y) ) stop("The length of family should be equal to ncol(y).")
  if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
  if( is.null(alpha)) alpha <- 1:9*0.1
  if( is.null(pf)) pf <- rep(1, ncol(x))
  if( length(unique(sapply(lambda, length))) > 1 ) stop("The lengths of lambda sequences should be equal to each other.")


  x <- as.matrix(x)
  if(std.x) x <- scale(x)
  y <- as.data.frame(y, stringsAsFactors = FALSE)


  n.pf <- sum(pf==0)
  if( n.pf > 0 ){
    wh.var <- (1:ncol(x))[-(1:n.pf)]
  } else {
    wh.var <- (1:ncol(x))
  }

  n <- nrow(x);
  p <- ncol(x)-n.pf;
  q <- ncol(y)
  nsub <- n*psub;

  sample.replace = FALSE
  if( psub == 1.0 ) sample.replace = TRUE


  lambda <- lapply( lambda, sort, decreasing = FALSE )
  nlambda <- min( sapply(lambda, length) )



  beta.array <- array( 0, dim = c(p, nlambda, ncol(y)),
                            dimnames = list(
                              paste0("", 1:p),
                              paste0("", seq_len(nlambda)),
                              paste0("", seq_len(ncol(y)))
                            ) )
  names(attributes(beta.array)$dimnames) <- c("variants", "lambda", "phenotype")






  for (l in 1:B) {
    start <- proc.time()

    wsub <- sample(n, nsub, replace = sample.replace)
    xsub <- x[wsub, , drop = F]
    ysub <- y[wsub, , drop = F]

    for (k in 1:ncol(y)) {

      FAMILY <- family[k]

      type.multinomial <- NULL
      ysub2 <- ysub[, k, drop= T]

      if (FAMILY == "multinomial"){
        ysub2 <- y.multinom(ysub[, k, drop= T])
        type.multinomial <- "grouped"
      }


      fit <- glmnet( x = xsub, y = ysub2,
                     alpha = alpha,
                     lambda = lambda[[k]],
                     family = FAMILY,
                     type.multinomial = type.multinomial,
                     standardize.response = std.y,
                     pf = pf, ... )

      fitted.beta <-
        switch(as.character(is.list(fit$beta)),
               "TRUE" = fit$beta[[1]],
               "FALSE" = fit$beta)
      fitted.beta <- fitted.beta[wh.var,]


      beta.array[, , k] <- beta.array[, , k] + as.numeric(fitted.beta != 0)

    }


    end <- proc.time()

    if(l == 1 & verbose) cat("\n The expected remained time is", (end - start)[3] * (B - l), "\n")

  }




  adj.selprob <- array(NA, dim(beta.array))

  selscore <- NULL
  for( lamb in 1:nlambda ){
    adj.selprob[,lamb,] <- enet2score( beta.array[,lamb,], delta = 1.0 )$score / B
    selscore <- cbind(selscore, apply( adj.selprob[,lamb,], 1, sum ))
  }
  colnames(selscore) <- paste0("lambda.", 1:nlambda)

  unadj.selprob <- beta.array / B






  out.list <- NULL
  for( h in 1:q ){
    out.list[[h]] <- unadj.selprob[,,h]
  }
  out.list[[q+1]] <- selscore


  out.array <- simplify2array(out.list)
  out.array <- aperm(out.array, c(1,3,2))

  dimnames(out.array) <- list(paste0(1:dim(out.array)[1]),
                              c(paste0("y.", 1:q),"total"),
                              paste0("lambda.", 1:dim(out.array)[3]))



  params <- list(
    x = x,
    y = y,
    family = family,
    lambda = lambda,
    alpha = alpha,
    B = B,
    psub = psub,
    pf = pf,
    std.x = std.x,
    std.y = std.y,
    ...
  )



  out <- list(selscore=out.array, params=params)

  attr(out, "call") <- Call

  return( out )

}












### internal functions =========================================================
## Internal function for adjusting selection probabilities
##
## @param Array An array of selection frequencies.
## @param delta A denoising parameter. As \code{delta} increases, the final selection score would decrease. Default is 1.0.
## @param weight Logical flag for adjusting individual selection probabilities.
##
## @return a list of results for adjusted selection probabilities.
enet2score <- function(Array, delta = 1.0, weight = TRUE){
  p <- nrow(Array)
  q <- ncol(Array)

  Array.sort <- apply(Array, 2, sort, decreasing=TRUE)
  Array.order <- apply(Array, 2, order, decreasing=TRUE)
  Array.sum <- apply(Array, 2, sum)

  eta <- floor( min(Array.sum) * delta )


  array.rank <- function(vec){
    count <- 0
    out <- NULL
    for( lev in unique(vec) ){
      count <- count + 1
      out[which( vec == lev )] <- count
    }
    out
  }

  Array.sort.rank <- apply( Array.sort, 2, array.rank )

  J <- S.tilde <- NULL
  for( k in 1:q ){

    J[k] <- min( which( cumsum( Array.sort[,k] ) >= eta ) )
    S.tilde[k] <- cumsum( Array.sort[,k] )[ J[k] ]

    if( sum( Array.sort.rank[,k] == Array.sort.rank[J[k], k] ) > 1 ){

      J[k] <- max( which( Array.sort.rank[,k] == Array.sort.rank[J[k], k] ) )
      S.tilde[k] <- cumsum( Array.sort[,k] )[ J[k] ]

    }

  }


  # Thresholding
  out <- NULL
  for( k in 1:q ){
    out[[k]] <- Array[ , k ][ Array.order[ 1:J[k], k ], drop=F ]
  }


  s.star <- matrix( 0, nrow = p, ncol = q )
  for( k in 1:q ){
    if( weight ){
      s.star[ as.numeric( names(out[[k]]) ), k ] <- out[[k]] * eta / max(1, S.tilde[k])
    } else {
      s.star[ as.numeric( names(out[[k]]) ), k ] <- out[[k]]
    }

  }


  # -- Start Validation --
  # sapply( 1:ncol(Array), function(k){
  #   Array[ 1:20, k ] * eta / S.tilde[k]
  # }) %>% apply(1, sum) %>% {./100}
  #
  # apply( s.star, 1, sum )[1:20] / 100
  # -- End Validation --



  list( array.sum = Array.sum,
        S.tilde = S.tilde,
        J = J,
        score = s.star )

}
