#' Linear combination of the fixed effect for each categories in multivariate count outcome.
#'
#' @param x.vec a value of binary covariate.
#' @param b vector of regression coefficients.
#' @param Des a design matrix
#' @param method either "dirmult" for Dirichlet - multinomial or "multl" for multinomial logistics mixed model.
#' @return a vector representing the result of multiplication between \code{Des} and \code{b} according to value of \code{x.vec} and the chosen method.
#' @examples
#' Des <- as.matrix(cbind(rep(1,6),rep(0:1,3),c(0,0,1,1,0,0), c(rep(0,4),rep(1,2)),c(0,0,0,1,0,0),c(rep(0,5),1)))
#' beta <- c(0,0.5,-3.5,-1.3,0.8,-0.15,-2,1)
#' FixEf(x.vec = 1,b = beta, Des = Des,method = "dirmult")

FixEf <- function(x.vec,b,Des,method ){
  gt <- c(Des %*% b)

  if (method == "dirmult"){
    if (x.vec == 0) {
      g = gt[seq(1,nrow(Des),by=2)]
    } else {
      g = gt[seq(2,nrow(Des),by=2)]}}
  if(method == "multl") {
    if (x.vec == 0) {
      g = gt[seq(1,nrow(Des),by=2)]
    } else {
      g = gt[seq(2,nrow(Des),by=2)]}}
  return(g)}
