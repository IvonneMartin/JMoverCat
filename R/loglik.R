#' The core function to compute the total likelihood of the joint model given the selected marginal distribution for the multivariate count outcome and the dataset
#'
#' @param params A vector of fixed effect and standard deviations of the random effects.
#' @param method An argument of either "dirmult" or "multl" defining the marginal distribution of the multivariate count outcome.
#' @param data A complete observation for all subjects including the multivariate count outcome, the covariates and the continuous outcomes for two time points.
#' @return The total loglikelihood values for the given parameters values.
#' @examples
#' # params should be in the format of params = {betas_C, betas_Y,uC's, uS's,theta or rho,uY,ev}, where
#' # The following example is for the multinomial logistics mixed model
#' betas_C <- c(0.5,3.5,1.3,0.8,0.15)
#' betas_Y <- c(2.3,0.1)
#' uC<- c(1,0.8)
#' uS <- c(0.5,0.6)
#' uY <- 0.9; ev <- 0.7; rho <- 0.1;
#' params <- c(betas_C,betas_Y, uC,uS,rho,uY,ev)
#' ll <- loglik(params,method = "multl",data)

loglik <- function(params,method,data){
  require(ecoreg)
  require(mvtnorm)
  require(Matrix)

  var.Y <- exp(params[(length(params)-1):length(params)])
  uY <- var.Y[1]
  ev <- var.Y[2]
  if (method == "dirmult"){
    var.Sigma <- exp(params[(length(params) - 8): (length(params) - 2)])
    uC <- var.Sigma[1:3]
    uS <- var.Sigma[4:6]
    theta <- var.Sigma[7]
    beta <- c(0,params[1:(length(params) - 9)])

    Sigma1 <- matrix(0,nrow=4,ncol=4)
    diag(Sigma1) <- c(uC^2 + uS^2,sum(uS^2) + uY^2)
    Sigma1[4,1] <- uS[1]^2
    Sigma1[4,2] <- uS[2]^2
    Sigma1[4,3] <- uS[3]^2
    Sigma <- as.matrix(Matrix::forceSymmetric(Sigma1,uplo = "L"))

    pts <- mgauss_hermite(3, mu=rep(0,3), sigma=Sigma)
    xGH <- pts$points
    wGH <- pts$weights

    index<-seq(1,nrow(data),by=1)
    w <- wGH
    x <- as.matrix(xGH )
    f <- match.fun(F)
    e1 <- sapply(index,integralComputation,b = beta,theta = theta,
                 Sigma = Sigma,ev,dataset = data,x = x,w = w,f = f,method = method)

    res <- sum(log(e1))
    return(res)
  }

  if (method == "multl"){
    rho <- cos(params[length(params) - 2])
    var.Sigma <- exp(params[(length(params) - 6): (length(params) - 3)])
    uC <- var.Sigma[1:2]
    uS <- var.Sigma[3:4]
    beta <- params[1:(length(params) - 7)]

    Sigma1 <- matrix(0,nrow=3,ncol=3)
    diag(Sigma1) <- c(uC^2 + uS^2,sum(uS^2) + uY^2)
    Sigma1[2,1] <- rho*uC[1]*uC[2]
    Sigma1[3,1] <- uS[1]^2
    Sigma1[3,2] <- uS[2]^2
    Sigma <- as.matrix(Matrix::forceSymmetric(Sigma1,uplo = "L"))

    pts <- mgauss_hermite(3, mu=rep(0,3), sigma=Sigma)

    xGH <- pts$points
    wGH <- pts$weights

    index<-seq(1,nrow(data),by=1)
    w <- wGH
    x <- as.matrix(xGH )
    f <- match.fun(Ffun)
    e1 <- sapply(index,integralComputation,b = beta,theta = 1,
                 Sigma = Sigma,ev,dataset = dat1,x = x,w = w,f = f,method = method)

    res <- sum(log(e1))
  }
  return(res)}
