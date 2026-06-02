#' simILD
#'
#' Generates a dataset in either long format (an individual is represented by at least
#' two rows) or short format (an individual is represented by at most two rows), with
#' random right censoring, for parametric and semi-parametric illness--death models
#' under constrained and unconstrained formulations.
#'
#' @param n Number of individuals in the dataset.
#'
#' @param censure \emph{Rate} parameter of the exponential distribution used to generate
#' random right censoring.
#'
#' @param lambda Vector of scale parameters for the baseline intensity functions
#' corresponding to the considered distributions (Weibull, Gompertz, and generalized
#' log-logistic). The columns respectively represent the generation parameters for the
#' transitions $1 -> 2$, $1 -> 3$, and $2 -> 3$. For the piecewise constant exponential
#' distribution, \code{lambda} is a matrix whose rows represent the parameters associated
#' with each interval. This parameter is not used for the piecewise linear approximation.
#'
#' @param alpha Vector of shape parameters for the baseline intensity functions
#' corresponding to the considered distributions (Weibull, Gompertz, and generalized
#' log-logistic). The columns respectively represent the generation parameters for the
#' transitions $1 -> 2$, $1 -> 3$, and $2 -> 3$. This parameter is not used for the
#' piecewise constant exponential distribution nor for the piecewise linear approximation.
#'
#' @param rho Vector of location parameters for the generalized log-logistic distribution.
#' The columns respectively represent the generation parameters for the transitions
#' $1 -> 2$, $1 -> 3$, and $2 -> 3$.
#'
#' @param beta Matrix of covariate coefficients of dimension
#' (\emph{number of covariates} $*$ 3).
#'
#' @param rates Matrix of intercept and slope parameters for the piecewise linear
#' approximation. The columns respectively represent the generation parameters for the
#' transitions $1 -> 2$, $1 -> 3$, and $2 -> 3$, while the rows correspond to the
#' parameters associated with each interval.
#'
#' @param covariable Number of covariates to simulate. Covariates are generated from
#' a standard normal distribution. The default value is 0.
#'
#' @param distribution Name of the model used for data generation.
#' Possible values are: Weibull (\code{"weibull"}), Gompertz
#' (\code{"gompertz"}), generalized log-logistic (\code{"gll"}),
#' piecewise constant exponential (\code{"pwexp"}), and piecewise
#' linear approximation (\code{"pla"}). The default model is Weibull.
#'
#' @param lat_long Indicates whether the dataset is returned in long format
#' (\code{lat_long = TRUE}) or short format (\code{lat_long = FALSE}).
#' The default value is \code{TRUE}.
#'
#' @param breakpoints Vector of cut points used for the piecewise constant
#' exponential distribution and the piecewise linear approximation.
#'
#' @return The function returns a dataset containing the variables:
#' \code{id} (individual identifier), \code{start} (entry time into state $k$),
#' \code{stop} (exit time from state $k$), \code{status} (indicates whether the
#' transition $k -> l$ is observed), \code{from} (starting state), \code{to}
#' (destination state), as well as the generated covariates when they are included
#' in the model.
#' @examples
#' library(flexsurv)
#'
#' simILD(
#'   n = 1000,
#'   censure = 0.001,
#'   lambda = c(0.5, 0.8, 0.8 * exp(-0.4)),
#'   alpha = c(2, 1.5, 1.5),
#'   beta = matrix(c(0.4, 0.2, 0.2), 1, 3, byrow = TRUE),
#'   covariable = 1,
#'   distribution = "weibull",
#'   lat_long = TRUE
#' )
#'
#'library(PWEXP)
#'
#' simILD(
#'   n = 1000,
#'   censure = 0.2,
#'   lambda = t(cbind(c(0.18, 0.2),
#'                     c(0.15, 0.2),
#'                     c(0.2, 0.25))),
#'   breakpoints = 3,
#'   beta = matrix(c(0.4, 2, 1.5), 1, 3, byrow = TRUE),
#'   covariable = 1,
#'   distribution = "pwexp",
#'   lat_long = TRUE
#' )
#'
#' @export
#'
#' @import survival
#' @import truncdist
#' @import eha
#' @import PWEXP
#' @importFrom stats rexp rnorm qnorm optim pchisq uniroot runif
#' @import dplyr
#' @importFrom flexsurv pweibullPH rweibullPH dweibullPH qweibullPH
simILD <- function(n, censure, lambda = NULL, alpha = NULL, rho = NULL,
                   beta = NULL, rates = NULL, covariable = 0, distribution = "weibull",
                   lat_long = TRUE, breakpoints =NULL){

  lois_dispo <- c("weibull", "gompertz", "gll", "pwexp","pla")

  if(distribution %in% c("pwexp")){
    stopifnot(all(as.vector(lambda) > 0))
  }else if(!(distribution %in% c("pla"))){
    stopifnot( length(lambda) == 3, length(alpha) == 3, all(lambda > 0), nrow(beta) == covariable)
  }
  if(distribution %in% c("weibull","gll")){
    stopifnot(all(alpha > 0))
  }
  if(distribution == "gll"){
    stopifnot(length(rho) == 3,all(rho >0))
  }
  loi <- match.arg(tolower(distribution), lois_dispo)
  nom_data <- paste0("sim_", loi)
  if(loi == "gll"){
    do.call(nom_data, list(N = n, lamc = censure, lambda = lambda, alpha = alpha, Betax = beta,
                           rho= rho, cov = covariable, lat_long = lat_long))
  }else if(loi == "pwexp"){
    do.call(nom_data, list(N = n, lamc = censure, lambda = lambda,
                           pts = breakpoints, Betax = beta,
                           cov = covariable, lat_long = lat_long))
  }else if(loi == "pla"){
    do.call(nom_data, list(N = n, lamc = censure, breakpoints = breakpoints,
                          rates = rates, Betax = beta, cov = covariable, lat_long = lat_long))
  }else if(loi == "weibull" | loi == "gompertz"){
    do.call(nom_data, list(N = n, lamc = censure, lambda = lambda, alpha = alpha,
                           Betax = beta, cov = covariable, lat_long = lat_long))
  }

}




#' @noRd
qtruncWPH = function(p, a = - Inf ,b = Inf, ...){
  if ( a >= b )
    stop( "argument a is greater than or equal to b" )
  G.a <- pweibullPH( a, ... )
  G.b <- pweibullPH( b, ... )
  if ( G.a == G.b ) {
    stop( "Trunction interval is not inside the domain of the density function" )
  }
  G = qweibullPH(pweibullPH(a, ...) - p*(pweibullPH(a, ...) - pweibullPH(b, ...)),  ...)
  return(G)
}
#' @noRd
rtruncWPH = function(n, a = - Inf, b = Inf , ...){
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtruncWPH(u, a = a, b = b, ...)
  return(x)
}
#' @noRd
sim_weibull <- function(N,lamc, lambda, alpha, Betax = NULL, cov = 0, lat_long = TRUE){
  outdata <- NULL
  lam12 <- lambda[1]
  lam13 <- lambda[2]
  lam23 <- lambda[3]
  alp12 <- alpha[1]
  alp13 <- alpha[2]
  alp23 <- alpha[3]
  for (i  in 1:N) {
    X <- NULL
    if(cov != 0){
      ux <- 1
      x <- rnorm(cov)
    }else{
      ux <- 0
      x <- 0
    }
    t12 <- flexsurv::rweibullPH(1, shape = alp12, scale = lam12*exp(ux*sum(Betax[,1]*x)))
    t13 <- flexsurv::rweibullPH(1,shape = alp13, scale = lam13*exp(ux*sum(Betax[,2]*x)))
    CC = rexp(1,rate = lamc)
    t1 <- min(t12,t13)
    s1<- ifelse(t12 < t13,2,3)
    if((t12 < t13) && (t12 < CC)){
      t23 <- rtruncWPH(1, a = t1, b = Inf, shape = alp23, scale = lam23*exp(ux*sum(Betax[,3]*x)))
      s23 <- 3
      start <- c(0,0,t1)
      stop <- c(t1,t1,min(t23, CC))
      status <- c(1,0, ifelse(t23 < CC, 1,0))
      from <- c(1,1,s1)
      to <- c(s1,s23,s23)
      for.etm <- c(1,0,1)
      if (cov != 0){
        X <- matrix(rep(x,3),3,cov, byrow = TRUE)
      }
    }else{
      start <- rep(0,2)
      stop <- rep(min(t1,CC),2)
      status <- c(ifelse(t1 < CC, 1,0),0)
      from <- c(1,1)
      to <- c(s1,2)
      for.etm <- c(1,0)
      if(cov != 0){
        X <- matrix(rep(x,2),2,cov, byrow = TRUE)
      }
    }
    id <- rep(i, length(start))
    if(isTRUE(lat_long)){
      if(cov!=0){
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm,X))
      }else{
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm))
      }
    }else{
      if(cov!=0){
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm,X))
        outdata <- outdata[outdata$for.etm == 1, ]
      }else{
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm))
        outdata <- outdata[outdata$for.etm == 1, ]
      }
    }
  }
  outdata <- outdata[order(outdata$id, outdata$from, outdata$to),]
  return(data.frame(outdata))
}
#' @noRd
sim_gompertz <- function(N,lamc, lambda, alpha, Betax = NULL, cov = 0, lat_long = TRUE){
  a12 <- lambda[1]
  a13 <- lambda[2]
  a23 <- lambda[3]
  r12 <- alpha[1]
  r13 <- alpha[2]
  r23 <- alpha[3]
  outdata <- NULL
  for (i  in 1:N) {
    X <- NULL
    if(cov != 0){
      ux <- 1
      x <- rnorm(cov)
    }else{
      ux <- 0
      x <- 0
    }
    t12 <- eha::rgompertz(1, shape = a12*exp(ux*sum(Betax[,1]*x)), rate = r12, param = "rate")
    t13 <- eha::rgompertz(1, shape = a13*exp(ux*sum(Betax[,2]*x)), rate = r13, param = "rate")
    CC <- rexp(1, rate = lamc)
    t1 <- min(t12,t13)
    s1<- ifelse(t12 < t13,2,3)
    if((t12 < t13) && (t12 < CC)){
      t23 <- rtrunc(1,spec = "gompertz", a = t1, b = Inf, shape = a23*exp(ux*sum(Betax[,3]*x)),
                    rate = r23, param = "rate")
      s23 <- 3
      start <- c(0,0,t1)
      stop <- c(t1,t1,min(t23, CC))
      status <- c(1,0, ifelse(t23 < CC, 1,0))
      from <- c(1,1,s1)
      to <- c(s1,s23,s23)
      for.etm <- c(1,0,1)
      if (cov != 0){
        X <- matrix(rep(x,3),3,cov, byrow = TRUE)
      }
    }else{
      start <- rep(0,2)
      stop <- rep(min(t1,CC),2)
      status <- c(ifelse(t1 < CC, 1,0),0)
      from <- c(1,1)
      to <- c(s1,2)
      for.etm <- c(1,0)
      if(cov != 0){
        X <- matrix(rep(x,2),2,cov, byrow = TRUE)
      }
    }
    id <- rep(i, length(start))
    if(isTRUE(lat_long)){
      if(cov!=0){
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm,X))
      }else{
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm))
      }
    }else{
      if(cov!=0){
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm,X))
        outdata <- outdata[outdata$for.etm == 1, ]
      }else{
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm))
        outdata <- outdata[outdata$for.etm == 1, ]
      }
    }
  }
  outdata <- outdata[order(outdata$id, outdata$from, outdata$to),]
  return(data.frame(outdata))
}
#' @noRd
pGLLPH<-function(x,lambda, alpha, rho){
  (1-((1+rho*(x)^alpha)^(-lambda/rho)))
}
#' @noRd
dGLLPH<-function(x,lambda, alpha, rho){
  ((alpha*lambda)*(x^(alpha-1)))/(1+(rho*(x)^alpha))^((lambda/rho)+1)
}
#' @noRd
hGLLPH<-function(x,lambda, alpha, rho){
  ((alpha*lambda)*(x^(alpha-1)))/(1+(rho*(x)^alpha))
}
#' @noRd
HGLLPH<- function(x,lambda, alpha, rho){
  cdf <- (1-((1+(rho*(x)^alpha))^(-(lambda/rho))))
  return(-log(1-cdf))
}
#' @noRd
qGLLPH <- function(lambda,alpha,rho, u)
{
  num1 <- rho/lambda
  num2 <- (1/(1-u))^num1
  overall <- (num2-1)/rho
  qfinal <- (overall)^(1/alpha)
  return(qfinal)
}
#' @noRd
rGLLPH<- function(n,lambda, alpha, rho){
  u= runif(n)
  s <-  qGLLPH(lambda,alpha,rho,u)
  return(s)
}
#' @noRd
qtruncGLLPH = function(p, a = - Inf ,b = Inf, ...){
  if ( a >= b )
    stop( "argument a is greater than or equal to b" )
  G.a <- pGLLPH( a, ... )
  G.b <- pGLLPH( b, ... )
  if ( G.a == G.b ) {
    stop( "Trunction interval is not inside the domain of the density function" )
  }
  G = qGLLPH(pGLLPH(a, ...) - p*(pGLLPH(a, ...) - pGLLPH(b, ...)),  ...)
  return(G)
}
#' @noRd
rtruncGLLPH = function(n, a = - Inf, b = Inf , ...){
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtruncGLLPH(u, a = a, b = b, ...)
  return(x)
}
#' @noRd
sim_gll <- function(N,lamc, lambda, alpha, rho, Betax = NULL, cov = 0, lat_long = TRUE){
  k12 <- lambda[1]
  k13 <- lambda[2]
  k23 <- lambda[3]
  a12 <- alpha[1]
  a13 <- alpha[2]
  a23 <- alpha[3]
  e12 <- rho[1]
  e13 <- rho[2]
  e23 <- rho[3]
  outdata <- NULL
  for (i  in 1:N) {
    X <- NULL
    if(cov != 0){
      ux <- 1
      x <- rnorm(cov)
    }else{
      ux <- 0
      x <- 0
    }
    t12 <- rGLLPH(1, lambda  = k12*exp(ux*sum(Betax[,1]*x)), alpha = a12, rho = e12)
    t13 <- rGLLPH(1, lambda = k13*exp(ux*sum(Betax[,2]*x)), alpha = a13, rho = e13)
    CC <- rexp(1, rate = lamc)
    t1 <- min(t12,t13)
    s1<- ifelse(t12 < t13,2,3)
    if((t12 < t13) && (t12 < CC)){
      t23 <- rtruncGLLPH(1, a = t1, b = Inf,lambda = k23*exp(ux*sum(Betax[,3]*x)), alpha = a23, rho = e23)
      s23 <- 3
      start <- c(0,0,t1)
      stop <- c(t1,t1,min(t23, CC))
      status <- c(1,0, ifelse(t23 < CC, 1,0))
      from <- c(1,1,s1)
      to <- c(s1,s23,s23)
      for.etm <- c(1,0,1)
      if (cov != 0){
        X <- matrix(rep(x,3),3,cov, byrow = TRUE)
      }
    }
    else{
      start <- rep(0,2)
      stop <- rep(min(t1,CC),2)
      status <- c(ifelse(t1 < CC, 1,0),0)
      from <- c(1,1)
      to <- c(s1,2)
      for.etm <- c(1,0)
      if(cov != 0){
        X <- matrix(rep(x,2),2,cov, byrow = TRUE)
      }
    }
    id <- rep(i, length(start))
    if(isTRUE(lat_long)){
      if(cov!=0){
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm,X))
      }else{
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm))
      }
    }else{
      if(cov!=0){
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm,X))
        outdata <- outdata[outdata$for.etm == 1, ]
      }else{
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm))
        outdata <- outdata[outdata$for.etm == 1, ]
      }
    }
  }
  outdata <- outdata[order(outdata$id, outdata$from, outdata$to),]
  return(data.frame(outdata))
}
#' @noRd
indi_vec = function(x, beta){
  p = sapply(x, function (x){ifelse(x >= beta[1] & beta[2] > x, 1,0)})
  return(p)
}
#' @noRd
int1 = function(t, beta){
  if(t >  beta[2]){
    p = beta[2] - beta[1]
  }else if(beta[1] <= t & t < beta[2] ){
    p = t - beta[1]
  }else{
    p = 0
  }
  return(p)
}
#' @noRd
int1_vec = function(t, beta){
  p = sapply(t, int1, beta = beta )
  return(p)
}
#' @noRd
int2 = function(s,t, beta){
  if(s <= beta[1] & beta[2] < t ){
    p = beta[2] - beta[1]
  }else if(beta[1] <= s & s < beta[2]  &  beta[2] < t){
    p = beta[2] - s
  }else if(beta[1] <= s & t < beta[2]){
    p = t - s
  }else if(s <= beta[1] & beta[1] < t  &  t < beta[2]){
    p = t -  beta[1]
  }else{
    p = 0
  }
  return(p)
}
#' @noRd
int2_vec = function(s, t, beta){
  p = NULL
  for(i in 1:length(t)){
    if(s[i] <= beta[1] & beta[2] < t[i] ){
      p = c(p, beta[2] - beta[1])
    }else if(beta[1] <= s[i] & s[i] < beta[2]  &  beta[2] < t[i]){
      p = c(p, beta[2] - s[i])
    }else if(beta[1] <= s[i] & t[i] < beta[2]){
      p = c(p, t[i] - s[i])
    }else if(s[i] <= beta[1] & beta[1] < t[i]  &  t[i] < beta[2]){
      p =  c( p, t[i] -  beta[1])
    }else{
      p = c(p,0)
    }
  }
  return(p)
}
#' @noRd
sim_pwexp <- function(N,lamc, lambda, pts, Betax = NULL, cov = 0, lat_long = TRUE){
  lam12 <- as.vector(lambda[1,])
  lam13 <- as.vector(lambda[2,])
  lam23 <- as.vector(lambda[3,])
  outdata <- NULL
  for (i  in 1:N) {
    X <- NULL
    if(cov != 0){
      ux <- 1
      x <- rnorm(cov, 0, 1)
    }else{
      ux <- 0
      x <- 0
    }
    t12 <- PWEXP::rpwexp(1, rate= lam12*exp(ux*sum(Betax[,1]*x)), breakpoint = pts)
    t13 <- PWEXP::rpwexp(1, rate= lam13*exp(ux*sum(Betax[,2]*x)), breakpoint = pts)
    CC <- rexp(1, rate = lamc)
    t1 <- min(t12,t13)
    s1 <- ifelse(t12 < t13,2,3)
    if((t12 < t13) && (t12 < CC)){
      t23 <- rtrunc(1,spec = "pwexp", a = t1, b = 100000, rate =lam23*exp(ux*sum(Betax[,3]*x)),
                  breakpoint = pts)
      s23 <- 3
      start <- c(0,0,t1)
      stop <- c(t1,t1,min(t23, CC))
      status <- c(1,0, ifelse(t23 < CC, 1,0))
      from <- c(1,1,s1)
      to <- c(s1,s23,s23)
      for.etm <- c(1,0,1)
      if (cov != 0){
        X <- matrix(rep(x,3),3,cov, byrow = TRUE)
      }
    }else{
      start <- rep(0,2)
      stop <- rep(min(t1,CC),2)
      status <- c(ifelse(t1 < CC, 1,0),0)
      from <- c(1,1)
      to <- c(s1,2)
      for.etm <- c(1,0)
      if(cov != 0){
        X <- matrix(rep(x,2),2,cov, byrow = TRUE)
      }
    }
    id <- rep(i, length(start))
    if(isTRUE(lat_long)){
      if(cov!=0){
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm,X))
      }else{
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm))
      }
    }else{
      if(cov!=0){
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm,X))
        outdata <- outdata[outdata$for.etm == 1, ]
      }else{
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm))
        outdata <- outdata[outdata$for.etm == 1, ]
      }
    }
  }
  outdata <- outdata[order(outdata$id, outdata$from, outdata$to),]
  return(data.frame(outdata))
}
#' @noRd
Hpla <- function (t, breakpoints, rates){
  lent <- length(t)
  cumhaz <- vector(mode = "numeric", length = lent)
  for (i in 1:lent) {
    int <- findInterval(t[i], breakpoints, all.inside = TRUE)
    z1 <- t[i] - breakpoints[1:int]
    z2 <- (t[i])^2-(breakpoints[1:int])^2
    exposure1 <- c(-diff(z1), z1[int])
    exposure2 <- c(-diff(z2), z2[int])
    kk <- as.vector(rates[1:int,])*c(exposure1,0.5*exposure2)
    cumhaz[i] <- sum(kk)
  }
  return(cumhaz)
}
#' @noRd
HHpla <- function (t, s, breakpoints, rates){
  lent <- length(t)
  cumhaz <- vector(mode = "numeric", length = lent)
  for (i in 1:lent) {
    breakp <- breakpoints
    int <- findInterval(t[i], breakpoints, all.inside = TRUE)
    int2 <- findInterval(s[i], breakpoints, all.inside = TRUE)
    breakp[int2] <- s[i]
    breakp[int+1] <- t[i]
    breakp <- breakp[int2:(int+1)]
    z1 <- t[i] - breakp
    z2 <- (t[i])^2-(breakp)^2
    exposure1 <- c(-diff(z1), z1[length(z1)])[-length(z1)]
    exposure2 <- c(-diff(z2), z2[length(z2)])[-length(z1)]
    kk <- as.vector(rates[int2:int,])*c(exposure1,0.5*exposure2)
    cumhaz[i] <- sum(kk)
  }
  return(cumhaz)
}
#' @noRd
hpla <- function(t, breakpoints, rates) {
  N <- nrow(rates)
  pp <- numeric(length(t))
  indices <- findInterval(t, breakpoints,all.inside = TRUE )
  for (i in 1:length(t)) {
    if (indices[i] > 0 && indices[i] <= N) {
      p <- sum(rates[indices[i], ] * c(1, t[i]))
      pp[i] <- p
    } else {
      pp[i] <- 0
    }
  }
  return(pp)
}
#' @noRd
ppla <- function (t, breakpoints, rates){
  lent <- length(t)
  cumhaz <- vector(mode = "numeric", length = lent)
  for (i in 1:lent) {
    int <- findInterval(t[i], breakpoints, all.inside = TRUE)
    z1 <- t[i] - breakpoints[1:int]
    z2 <- (t[i])^2-(breakpoints[1:int])^2
    exposure1 <- c(-diff(z1), z1[int])
    exposure2 <- c(-diff(z2), z2[int])
    kk <- as.vector(rates[1:int,])*c(exposure1,0.5*exposure2)
    cumhaz[i] <- sum(kk)
  }
  return(1- exp(-cumhaz))
}
#' @noRd
pw_root <- function(t, breakpoints,rates, uu){
  aa <-  Hpla(t, breakpoints,rates) + log(1-uu)
  return(aa)
}
#' @noRd
gte <- function(x, rat, bre){
  if(x==1){
    return(rat[,x]*diff(bre^x))
  }else{
    return(0.5*rat[,x]*diff(bre^x))
  }
}
#' @noRd
ff <- function(rates, breakpoints){
  sapply(1:(ncol(rates)), gte, rat=rates, bre=breakpoints)
}
#' @noRd
qpla <- function(breakpoints,rates, u){
  interval = c(breakpoints[1], breakpoints[length(breakpoints)])
  V4 <- cumsum(apply(ff(rates, breakpoints),1,sum))
  PP = NULL
  for (i in 1:length(u)) {
    index <- findInterval(-log(1-u[i]),V4)
    if(length(breakpoints)-1==index){
      interval <- c(breakpoints[length(breakpoints)-1],100+0.01)
    }else if(index==0){
      interval <- c(breakpoints[1],100+0.01)
    }else{
      interval <- c(breakpoints[index+1],100+0.01)
    }
    PP = c(PP,uniroot(f=pw_root,interval=interval,breakpoints,rates,uu=u[i])$root)
  }
  return(PP)
}
#' @noRd
rpla <-function(n,breakpoints,rates){
  u= runif(n)
  sim <- qpla(breakpoints,rates,u)
  return(sim)
}
#' @noRd
qtruncpla <- function(p, aa = - Inf ,bb = Inf, ...){
  if ( aa >= bb )
    stop( "argument a is greater than or equal to b" )
  G.a <- ppla( aa, ... )
  G.b <- ppla( bb, ... )
  if ( G.a == G.b ) {
    stop( "Trunction interval is not inside the domain of the density function" )
  }
  G = qpla(..., ppla(aa, ...) - p*(ppla(aa, ...) - ppla(bb, ...)))
  return(G)
}
#' @noRd
rtruncpla = function(n, aa = - Inf, bb = Inf , ...){
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtruncpla(u, aa = aa, bb = bb, ...)
  return(x)
}
#' @noRd
sim_pla <- function(N,lamc, breakpoints, rates,Betax = NULL, cov = 0, lat_long = TRUE){
  rates12 = rates[[1]]
  rates13 = rates[[2]]
  rates23 = rates[[3]]
  outdata <- NULL
  for (i  in 1:N) {
    X <- NULL
    if(cov != 0){
      ux <- 1
      x <- rnorm(cov, 0, 1)
    }else{
      ux <- 0
      x <- 0
    }
    t12 <- rpla(1, breakpoints = breakpoints, rates12*exp(ux*sum(Betax[,1]*x)))
    t13 <- rpla(1, breakpoints = breakpoints, rates13*exp(ux*sum(Betax[,2]*x)))
    CC <- rexp(1,rate = lamc)
    t1 <- min(t12,t13)
    s1<- ifelse(t12 < t13,2,3)
    if((t12 < t13) && (t12 < CC)){
      t23 <- rtruncpla(1, aa = t1, bb = 1000000, breakpoints = breakpoints
                       , rates23*exp(ux*sum(Betax[,3]*x)))
      s23 <- 3
      start <- c(0,0,t1)
      stop <- c(t1,t1,min(t23, CC))
      status <- c(1,0, ifelse(t23 < CC, 1,0))
      from <- c(1,1,s1)
      to <- c(s1,s23,s23)
      for.etm <- c(1,0,1)
      if (cov != 0){
        X <- matrix(rep(x,3),3,cov, byrow = TRUE)
      }
    }
    else{
      start <- rep(0,2)
      stop <- rep(min(t1,CC),2)
      status <- c(ifelse(t1 < CC, 1,0),0)
      from <- c(1,1)
      to <- c(s1,2)
      for.etm <- c(1,0)
      if(cov != 0){
        X <- matrix(rep(x,2),2,cov, byrow = TRUE)
      }
    }
    id <- rep(i, length(start))
    if(isTRUE(lat_long)){
      if(cov!=0){
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm,X))
      }else{
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm))
      }
    }else{
      if(cov!=0){
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm,X))
        outdata <- outdata[outdata$for.etm == 1, ]
      }else{
        outdata <- rbind(outdata, data.frame(id, start, stop, status, from, to, for.etm))
        outdata <- outdata[outdata$for.etm == 1, ]
      }
    }
  }
  outdata <- outdata[order(outdata$id, outdata$from, outdata$to),]
  return(data.frame(outdata))
}


