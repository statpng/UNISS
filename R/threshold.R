#' @name cutoff
#' @title A Cut-Off of Unified Selection Score
#' @description
#' When the expected number of falsely selected variants is given by \eqn{\theta}, the corresponding cut-off of unified selection score is calculated. Then, we conclude that the variants with unified selection score greater than the obtained cut-off value are potentially associated with multiple phenotypes, and the number of falsely selected variants among them is expected to be less than \eqn{\theta}. In the data frame of the output, each row index represent the value of \eqn{\theta}, that is the first row is \eqn{\theta=1}, the second row is \eqn{\theta=2}, and so on.
#'
#' @param params The parameters used in \code{main} function. If you already run the \code{main} function, it then provides the list of params.
#' @param nperm The number of permutation replicates. Default is 100.
#' @param seed An integer value for reproducible outcome. Default is 1 that does not have special meaning.
#'
#' @usage
#' cutoff(params, nperm = 100, seed = 1)
#'
#'
#' @references
#' Kim, K., Koo, J., & Sun, H. (2020). An empirical threshold of selection probability for analysis of high-dimensional correlated data. \emph{Journal of Statistical Computation and Simulation}, 90(9), 1606-1617.
#'
#'
#' @return
#' A data frame containing cut-off values of unified selection score. Each row represents the expected number (\eqn{\theta}) of falsely selected variants and each column represents a phenotype.
#'
#'
#' @import glmnet
#'
#'
#' @export cutoff
cutoff <- function(params, nperm = 100, seed = 1) {

  Call <- match.call()

  n <- nrow( params$x )
  p <- sum( params$pf == 1 )
  q <- ncol( params$y )

  nalpha <- length( params$alpha )
  nlambda <- min(sapply( params$lambda, length ))


  cutoff.array <- array(0, dim=c(p, q+1, nlambda) )
  maxFD = floor(0.1*p)
  selscore.perm.array <- array(0, dim=c(maxFD, q+1, nperm, 2, nlambda),
                      dimnames = list(paste0("top", 1:maxFD),
                                      c(paste0("y", 1:q), "total"),
                                      paste0("perm", 1:nperm),
                                      c("sort", "order"),
                                      paste0("lambda.",1:nlambda)) )

  print("Permutation get started!")
  perm.params <- params

  for( perm.i in 1:nperm ){

    if( perm.i == 1 ) start <- proc.time()
    # permutation -------------------------------------------------------------
    set.seed( seed + perm.i - 1 )


    perm.params$y <- params$y[sample(n),]
    perm.params$verbose <- FALSE

    fit.perm <- do.call("uniss", perm.params)

    # selscore.perm <- apply( fit.perm$selscore, 1:2, max )
    selscore.perm <- fit.perm$selscore
    # colnames(selscore.perm) <- c(paste0("y.", 1:ncol(params$y)), "total")


    selscore.perm.array[,,perm.i,1,] <- apply( selscore.perm, 2:3, function(sck) sort(sck, decreasing=TRUE)[1:maxFD] )
    selscore.perm.array[,,perm.i,2,] <- apply( selscore.perm, 2:3, function(sck) order(sck, decreasing=TRUE)[1:maxFD] )


    for( fdfd in 1:p ){
      selscore.perm.top <- apply( selscore.perm, 2:3, function(sck) sort(sck, decreasing=TRUE)[ fdfd ] ) / nperm

      cutoff.array[fdfd, ,] <- cutoff.array[fdfd, ,] + selscore.perm.top

    }

    if( perm.i == 1 ) cat("The expected time left = ", (proc.time() - start)["elapsed"] * nperm / 60, " (min) \n" )
  }

  dimnames(cutoff.array)[2:3] <- dimnames(selscore.perm)[2:3]
  attr(cutoff.array, "call") <- Call


  return( cutoff.array )

}


