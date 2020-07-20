#' Compute the marginal loglikelihood of the multivariate count outcome.
#'
#' @param Ct A vector of categorical count.
#' @param g A function of expected values for each category of the multivariate count outcome.
#' @param method An argument, either "dirmult" for the Dirichlet-multinomial mixed model and "multl" for the multinomial logistics mixed model.
#' @return The loglikelihood values for the marginal distribution of the multivariate count outcome.

CNF <- function(Ct,g,method){
  gs <- sum(g)
  ys <- sum(Ct)

  if (method == "dirmult") {
    val <- lgamma(ys+1) + lgamma(gs) - lgamma(ys+gs) +
      sum(lgamma(Ct+g) - lgamma(g) - lgamma(Ct+1))
  }

  if (method == "multl"){
    prob <- g/gs
    val <- lgamma(ys + 1) + sum(Ct * log(prob) -
                                  lgamma(Ct + 1))
  }
  return(val)}
