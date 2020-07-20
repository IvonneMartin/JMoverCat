#' Nodes and weights for multivariate Gauss - Hermite quadratures
#'
#' @param n the number of quadratures
#' @param mu mean vector.
#' @param sigma covariance matrix
#' @return a set of points and nodes of multivariate Gauss-Hermite quadratures corresponding to the mean vector \code{mu} and \code{sigma}
#' @examples
#' Sigma <- cbind(c(1.5,-0.6),c(-0.6,4))
#' mgauss_hermite(n = 3, mu = c(1.5,0), sigma = Sigma)

mgauss_hermite <- function(n, mu, sigma) {
  require("ecoreg")
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")

  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)

  return(list(points=pts, weights=wts))
}
