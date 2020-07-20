#' A loglikelihood value of a subject given a random effect \code{z}
#'
#' @param z A vector of random effect .
#' @param b A vector of fixed effect parameters.
#' @param theta an overdispersion parameter for the Dirichlet - multinomial mixed model.
#' @param Sigma The matrix of covariance of the random effect.
#' @param ev the standard deviation for residuals in continous outcome.
#' @param data a vector of observation for a subject.
#' @param method An argument of either "dirmult" for Dirichlet-multinomial mixed model or "multl" for multinomial logistics mixed model.

F <- function(z,b,theta,Sigma,ev,data,method){
  Yt <- data

  if (method == "dirmult"){
    eta1 <- FixEf(Yt[7], b[1:6],Des,method)
    eta2 <- FixEf(Yt[8], b[1:6],Des,method)
    eta1 <- (1/theta)*exp(eta1 + c(z[1:3]))
    eta2 <- (1/theta)*exp(eta2 + c(z[1:3]))
    val <- CNF(Yt[1:3],eta1,method) + CNF(Yt[4:6],eta2,method) +
      dnorm(Yt[9],mean=c(1,Yt[7])%*%b[7:8] + z[4],sd = ev,log=TRUE)+
      dnorm(Yt[10],mean=c(1,Yt[8])%*%b[7:8] + z[4],sd = ev,log=TRUE)+
      dmvnorm(z, mean = rep(0,4), sigma = Sigma, log=TRUE)
    return(val)
  }

  if (method == "multl"){
    eta1 <- FixEf(Yt[7], b[1:4],Des,method)
    eta2 <- FixEf(Yt[8], b[1:4],Des,method)
    eta1 <- exp(eta1 + z[1:2])
    eta2 <- exp(eta2 + z[1:2])
    Eta1 <- c(1,eta1)
    Eta2 <- c(1,eta2)
    val <- CNF(Yt[1:3],Eta1,method) + CNF(Yt[4:6],Eta2,method) +
      dnorm(Yt[9],mean=c(1,Yt[7])%*%b[5:6] + z[3],sd = ev,log=TRUE)+
      dnorm(Yt[10],mean=c(1,Yt[8])%*%b[5:6] + z[3],sd = ev,log=TRUE)+
      dmvnorm(z, mean = rep(0,3), sigma = Sigma, log=TRUE)
  }
  return(val)
}
