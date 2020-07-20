#' Generate observation for joint model between multivariate count outcome and continuous outcome
#' The marginal distribution for the multivariate count outcome could be either Dirichlet - multinomial mixed model or multinomial logistic mixed model
#' @param N A number of subjects.
#' @param method An argument defining the marginal distribution of the multivariate count outcome, either "dirmult" or "multl".
#' @param var.level An argument of either "high" or "low" as defined in the simulation setting
#' @return A dataset in wideformat consisting the observations of the multivariate count outcome for two time points, the covariate and the continuous outcome.
#' @examples
#' DataGenerate(100,method = "dirmult",var.level = "high")



DataGenerate <- function(N,method,var.level){
  require(HMP)
  require(mvtnorm)
  require(Matrix)

  Q <- 3
  S <- 2000
  uY <- 0.9
  eps <- 0.7

  if (method == "dirmult") {
    th <- 0.1
    uC <- c(1,1,0.9)
    if(var.level == "high") {
      uS <- c(1,1,0.9) }
    if (var.level =="low") {
      uS <- c(1,0.5,0.6)}

    Sigma1 <- matrix(0,nrow = 4,ncol = 4)
    diag(Sigma1) <- c(uC^2 + uS^2, sum(uS^2) + uY^2)
    Sigma1[4,1] <- uS[1]^2
    Sigma1[4,2] <- uS[2]^2
    Sigma1[4,3] <- uS[3]^2
    Sigma <- as.matrix(Matrix::forceSymmetric(Sigma1,uplo = "L"))

    Des <- as.matrix(cbind(rep(1,6),rep(0:1,3),c(0,0,1,1,0,0),
                           c(rep(0,4),rep(1,2)),c(0,0,0,1,0,0),c(rep(0,5),1)))
    beta <- c(0,0.5,-3.5,-1.3,0.8,-0.15,-2,1)
    Umc <- rmvnorm(N,mean=rep(0,4),sigma=Sigma)
    Eps <- rnorm(2*N, mean = 0, sd = eps)
    X <- cbind(rbinom(N,size=1,prob = 0.5),rbinom(N,size = 1,prob = 0.5))

    XB1 <- t(sapply(X[,1],FixEf,b = beta[1:6],Des,method))
    XB.tilde1 <- (1/th)*exp(cbind(XB1+Umc[,1:3]))
    C1 <- t(apply(XB.tilde1,1,Dirichlet.multinomial,Nrs = S))
    Y1 <- cbind(1,X[,1])%*%beta[7:8] + Umc[,4] + Eps[1:N]


    XB2 <- t(sapply(X[,2],FixEf,b = beta[1:6],Des,method))
    XB.tilde2 <- (1/th)*exp(cbind(XB2+Umc[,1:3]))
    C2 <- t(apply(XB.tilde2,1,Dirichlet.multinomial,Nrs = S))
    Y2 <- cbind(1,X[,2])%*%beta[7:8] + Umc[,4] + Eps[(N+1):(2*N)]

    data1 <-  cbind(C1,C2,X,Y1,Y2)
  }

  if(method == "multl"){

    rho <- 0.1
    uC <- c(1,0.8)
    if(var.level == "high") {
      uS <- c(1,0.9)
    }
    if (var.level == "low") {
      uS <- c(0.5,0.6)
    }

    Sigma1 <- matrix(0,nrow = 3,ncol = 3)
    diag(Sigma1) <- c(uC^2 + uS^2, sum(uS^2) + uY^2)
    Sigma1[2,1] <- rho*uC[1]*uC[2]
    Sigma1[3,1] <- uS[1]^2
    Sigma1[3,2] <- uS[2]^2
    Sigma <- as.matrix(Matrix::forceSymmetric(Sigma1,uplo = "L"))

    Des <- cbind(c(rep(1,2),rep(0,2)),c(rep(0,2),rep(1,2)),
                 c(0,1,0,0),c(0,0,0,1))
    beta <- c(-3.5,-1.3,0.8,-0.15,-2,0.1)
    Umc <- rmvnorm(N,mean = rep(0,3),sigma = Sigma)
    Eps <- rnorm(2*N,mean = 0,sd = eps)

    eta1 <- t(sapply(X[,1],FixEf,b = beta[1:4],Des,method))
    eta2 <- t(sapply(X[,2],FixEf,b = beta[1:4],Des,method))

    Eta1 <- cbind(1, exp(eta1 + Umc[,c(1,2)]))
    SumEta1 <- apply(Eta1,1,sum)
    Eta1 <- Eta1/SumEta1

    Eta2 <- cbind(1, exp(eta2 + Umc[,c(1:2)]))
    SumEta2 <- apply(Eta2,1,sum)
    Eta2 <- Eta2/SumEta2

    C1 <- t(apply(Eta1,1,rmultinom,n=1,size = S))
    C2 <- t(apply(Eta2,1,rmultinom,n=1,size = S))

    Y1 <- cbind(1,X[,1])%*%beta[5:6] + Umc[,3] + Eps[1:N]
    Y2 <- cbind(1,X[,2])%*%beta[5:6] + Umc[,3] + Eps[(N+1):(2*N)]

    data1 <- cbind(C1,C2,X[,1],X[,2],Y1,Y2)
  }
  colnames(data1) <- c("C11","C21","C31","C12","C22","C32","X1","X2","Y1","Y2")
  return(data1)
}
