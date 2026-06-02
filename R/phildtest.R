#' phildtest
#'
#' Performs a likelihood ratio test between the constrained and unconstrained
#' illness--death models, considering Weibull, Gompertz, and generalized
#' log-logistic distributions.
#'
#' @param lambda Vector of initial values for the scale parameters of the
#' baseline intensity functions for the considered distributions (Weibull,
#' Gompertz, and generalized log-logistic).
#'
#' @param alpha Vector of initial values for the shape parameters of the
#' considered distributions (Weibull, Gompertz, and generalized log-logistic).
#'
#' @param rho Vector of initial values for the location parameters of the
#' generalized log-logistic distribution.
#'
#' @param beta Initial value of the proportionality parameter.
#' This argument is used only for estimating the constrained
#' illness--death model.
#'
#' @param betax Matrix of initial values for covariate coefficients.
#'
#' @param distribution Name of the model used. Possible values are:
#' Weibull (\code{"weibull"}), Gompertz (\code{"gompertz"}), and
#' generalized log-logistic (\code{"gll"}). The default model is Weibull.
#'
#' @param dataset Dataset used for estimation.
#'
#' @param NomCovariable Vector containing the names of the covariates
#' included in the dataset.
#'
#' @return The function returns the likelihood ratio test statistic
#' along with the associated p-value.
#' @seealso phildII
#'
#' @export
phildtest <- function(lambda= NULL, alpha= NULL, rho = NULL,beta = NULL,
                      betax = NULL, distribution="weibull", dataset = NULL, NomCovariable= NULL){

  loi <- tolower(distribution)
  cat("\n", rep("=", 100), "\n", sep = "")
  if(loi == "gll"){
    cat(">> Step 1: Optimization of the restricted model (H0)...")
    incomplet <-  phild1(lambda = lambda[1:2], alpha = alpha[1:2], rho = rho[1:2],beta = beta, betax = betax,
                         distribution = distribution, nom_loglik = "logincomplet",
                         dataset = dataset,NomCovariable = NomCovariable)

    parInComplet <- incomplet[,1]
    fn_nameincomplet <- paste0("logincomplet", "_", loi)
    fn_nameincomplet <- get(fn_nameincomplet, mode = "function")
    H0G <- do.call(fn_nameincomplet, list(par =  parInComplet, data = dataset, NomCovariable=NomCovariable))

    cat(">> Step 2: Optimization of the full model (H1)...")
    complet <- phild1(lambda = lambda, alpha = alpha, rho = rho, beta = beta , betax = betax,
                      distribution = distribution, nom_loglik = "logcomplet",
                      dataset = dataset, NomCovariable = NomCovariable)
    parComplet <- complet[,1]
    fn_namecomplet <- paste0("logcomplet", "_", loi)
    fn_namecomplet <- get(fn_namecomplet, mode = "function")
    H1G <- do.call(fn_namecomplet, list(par= parComplet, data = dataset, NomCovariable=NomCovariable))
  }else if(loi == "weibull" | loi == "gompertz"){

    cat(">> Step 1: Optimization of the restricted model (H0)...")
    incomplet <-  phild1(lambda = lambda[1:2], alpha = alpha[1:2], beta = beta, betax = betax,
                         distribution = distribution, nom_loglik = "logincomplet",
                         dataset = dataset,NomCovariable = NomCovariable)
    parInComplet <- incomplet[,1]
    fn_nameincomplet <- paste0("logincomplet", "_", loi)
    fn_nameincomplet <- get(fn_nameincomplet, mode = "function")
    H0G <- do.call(fn_nameincomplet, list(par =  parInComplet, data = dataset, NomCovariable=NomCovariable))
    print(H0G)
    cat(">> Step 2: Optimization of the full model (H1)...")
    complet <- phild1(lambda = lambda, alpha = alpha, beta = beta , betax = betax,
                      distribution = distribution, nom_loglik = "logcomplet",
                      dataset = dataset, NomCovariable = NomCovariable)

    parComplet <- complet[,1]
    fn_namecomplet <- paste0("logcomplet", "_", loi)
    fn_namecomplet <- get(fn_namecomplet, mode = "function")
    H1G <- do.call(fn_namecomplet, list(par= parComplet, data = dataset, NomCovariable=NomCovariable))
    print(H1G)
  }

  cat(">> Step 3: Likelihood Ratio Test calculation...")
  ddl <- nrow(complet) - nrow(incomplet)
  DG <- 2*(H0G-H1G)
  pvalG <- 1 - pchisq(DG, ddl)
  cat("\n >>> Test completed successfully.\n")
  cat("\n", rep("=", 100), "\n", sep = "")
  structure(

    list(
      statistic = c(LRT = DG),
      p.value = as.numeric(pvalG),
      parameter = c(df = ddl),
      method = paste("Likelihood Ratio Test for Proportional Hazards assumption (", distribution, ")", sep=""),
      data.name = deparse(substitute(dataset)),
      alternative = "two-sided"
    ),
    class = "htest"
  )
}
#' @noRd
logincomplet_weibull <- function(par, data, NomCovariable=NULL){
  l12 <- par[1]
  l13 <- par[2]
  a12 <- par[3]
  a13 <- par[4]
  beta<- par[5]
  data <- data.frame(data)
  if(!is.null(NomCovariable)){
    X <- as.matrix(data[,NomCovariable])
    ncov <- ncol(X)
    data$lin12 <- as.vector(X %*% par[5+(1:ncov)])
    data$lin13 <- as.vector(X %*% par[5+((1+ncov):(2*ncov))])
  }else{
    data$lin12 <- data$lin13 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*(log(a13*l13*(stop)^(a13 - 1)) +  lin13)
                                                           - l13*((stop)^a13)*exp(lin13) ))
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(status*(log(a12*l12*(stop)^(a12-1)  ) + lin12)
                                                           - l12*((stop)^a12)*exp(lin12)))
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(status*(log(a13*l13*exp(beta)*(stop)^(a13 - 1) ) + lin13)
                                                           -l13*exp(beta)*(stop^a13 - start^a13 )*exp(lin13)))

  L <- L12 + L23 + L13
  return(-L)
}
#' @noRd
logincomplet_gompertz <- function(par,data, NomCovariable=NULL){
  a12 <- par[1]
  a13 <- par[2]
  r12 <- par[3]
  r13 <- par[4]
  beta <- par[5]
  data <- data.frame(data)
  if(!is.null(NomCovariable)){
    X <- as.matrix(data[, NomCovariable])
    ncov <- ncol(X)
    data$lin12 <- as.vector(X %*% par[5+(1:ncov)])
    data$lin13 <- as.vector(X %*% par[5+((1+ncov):(2*ncov))])
  }else{
    data$lin12 <- data$lin13 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*(log(a13*exp(r13*stop) ) + lin13 )
                                                           -(a13/r13)*(exp(r13*stop) - 1)*exp(lin13)) )
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(status*(log(a12*exp(r12*stop)) + lin12 )
                                                           -(a12/r12)*(exp(r12*stop) - 1)*exp(lin12)) )
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(status*(log(a13*exp(beta)*exp(r13*stop)) + lin13 )
                                                           - ((a13*exp(beta))/r13)*( exp(r13*stop) - exp(r13*start) )*exp(lin13) ))

  L <- L12 + L23 + L13
  return(-L)
}
#' @noRd
logincomplet_gll <- function(par,data, NomCovariable=NULL){
  l12 <- par[1]
  l13 <- par[2]
  a12 <- par[3]
  a13 <- par[4]
  r12 <- par[5]
  r13 <- par[6]
  beta <- par[7]
  data <- data.frame(data)
  if(!is.null(NomCovariable)){
    X <- as.matrix(data[,NomCovariable])
    ncov <- ncol(X)
    data$lin12 <- as.vector(X %*% par[7+(1:ncov)])
    data$lin13 <- as.vector(X %*% par[7+((1+ncov):(2*ncov))])
  }else{
    data$lin12 <- data$lin13 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*( log(a13*l13*(stop)^(a13 -1))
                                                                    - log(1 + r13*(stop)^a13) + lin13)
                                                           - (l13/r13)*log(1 + r13*(stop)^a13)*exp(lin13)  ))
  L12 <- with(data[(data$from == 1)&(data$to ==2), ], sum(status*( log(a12*l12*(stop)^(a12 -1))
                                                                   - log(1 + r12*(stop)^a12) + lin12)
                                                          - (l12/r12)*log(1 + r12*(stop)^a12)*exp(lin12)  ))
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum( status*( log(a13*l13*exp(beta)*(stop)^(a13 -1))
                                                                     - log(1 + r13*(stop)^a13) + lin13)
                                                            - (l13*exp(beta)/r13)* ( log(1 + r13*(stop)^a13)
                                                                                     - log(1 + r13*(start)^a13))*exp(lin13)  ))
  L <- L12 + L23 + L13
  return(-L)
}
#' @noRd
logcomplet_weibull <- function(par, data, NomCovariable = NULL){
  l12 <- par[1]
  l13 <- par[2]
  l23 <- par[3]
  a12 <- par[4]
  a13 <- par[5]
  a23 <- par[6]
  data <- data.frame(data)
  if(is.null(NomCovariable)){
    data$lin12 <- data$lin13 <- data$lin23 <- 0
  }else{
    X <- as.matrix(data[,NomCovariable])
    ncov <- ncol(X)
    data$lin12 <- as.vector(X %*% par[6+(1:ncov)])
    data$lin13 <- as.vector(X %*% par[6+((1+ncov):(2*ncov))])
    data$lin23 <- as.vector(X %*% par[(2*ncov+7):length(par)])
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*(log(a13*l13*(stop)^(a13 - 1)) + lin13)
                                                           - l13*((stop)^a13)*exp(lin13)))
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(status*(log(a12*l12*(stop)^(a12-1)) + lin12)
                                                           -  l12*((stop)^a12)*exp(lin12)))
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(status*(log(a23*l23*(stop)^(a23 - 1) ) + lin23)
                                                           - l23*(stop^a23 - start^a23 )*exp(lin23)))
  L <- L12 + L23 + L13
  return(-L)
}
#' @noRd
logcomplet_gompertz <- function(par,data, NomCovariable=NULL){
  a12 <- par[1]
  a13 <- par[2]
  a23 <- par[3]
  r12 <- par[4]
  r13 <- par[5]
  r23 <- par[6]
  data <- data.frame(data)
  if(!is.null(NomCovariable)){
    X <- as.matrix(data[,NomCovariable])
    ncov <- ncol(X)
    data$lin12 <- as.vector(X %*% par[6+(1:ncov)])
    data$lin13 <- as.vector(X %*% par[6+((1+ncov):(2*ncov))])
    data$lin23 <- as.vector(X %*% par[(2*ncov+7):length(par)])
  }else{
    data$lin12 <- data$lin13 <- data$lin23 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*(log(a13*exp(r13*stop) ) + lin13 )
                                                           -(a13/r13)*(exp(r13*stop) - 1)*exp(lin13) ) )
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(status*(log(a12*exp(r12*stop)) + lin12 )
                                                           -(a12/r12)*(exp(r12*stop) - 1)*exp(lin12)) )
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(status*(log(a23*exp(r23*stop)) + lin23  )
                                                           - (a23/r23)*( exp(r23*stop) - exp(r23*start) )*exp(lin23)  ))
  L <- L12 + L23 + L13
  return(-L)
}
#' @noRd
logcomplet_gll <- function(par,data, NomCovariable=NULL){
  l12 <- par[1]
  l13 <- par[2]
  l23 <- par[3]
  a12 <- par[4]
  a13 <- par[5]
  a23 <- par[6]
  r12 <- par[7]
  r13 <- par[8]
  r23 <- par[9]
  data <- data.frame(data)
  if(!is.null(NomCovariable)){
    X <- as.matrix(data[,NomCovariable])
    ncov <- ncol(X)
    data$lin12 <- as.vector(X %*% par[9+(1:ncov)])
    data$lin13 <- as.vector(X %*% par[9+((1+ncov):(2*ncov))])
    data$lin23 <- as.vector(X %*% par[(2*ncov+10):length(par)])
  }else{
    data$lin12 <- data$lin13 <- data$lin23 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*( log(a13*l13*(stop)^(a13 -1))
                                                                    - log(1 + r13*(stop)^a13) + lin13 )
                                                           - (l13/r13)*log(1 + r13*(stop)^a13)*exp(lin13) ))
  L12 <- with(data[(data$from == 1)&(data$to ==2), ], sum(status*( log(a12*l12*(stop)^(a12 -1))
                                                                   - log(1 + r12*(stop)^a12) + lin12 )
                                                          - (l12/r12)*log(1 + r12*(stop)^a12)*exp(lin12)  ))
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum( status*( log(a23*l23*(stop)^(a23 -1))
                                                                     - log(1 + r23*(stop)^a23) + lin23 )
                                                            - (l23/r23)* ( log(1 + r23*(stop)^a23)
                                                                           - log(1 + r23*(start)^a23))*exp(lin23) ))
  L <- L12 + L23 + L13
  return(-L)
}
#' @noRd
phild1 <- function(lambda= NULL, alpha = NULL, beta = NULL, rho = NULL, betax = NULL,
                   distribution= "weibull", nom_loglik = "logcomplet",dataset = NULL,NomCovariable= NULL){

  if(is.null(betax)){
    betavector <- NULL
  }else{
    if(nom_loglik == "logcomplet"){
      suffixes <- c("12", "13", "23")
      betavector <- as.vector(betax)
      beta_names <- paste0("beta", rep(1:nrow(betax), length(suffixes)), rep(suffixes, each = nrow(betax)))
    }
    if(nom_loglik == "logincomplet"){
      suffixes <- c("12", "13")
      betax <- matrix(betax[,1:2], nrow =nrow(betax), ncol = ncol(betax)-1)
      betavector <- as.vector(betax)
      beta_names <- paste0("beta", rep(1:nrow(betax), length(suffixes)), rep(suffixes, each = nrow(betax)))
    }
  }
  loi <- tolower(distribution)
  stopifnot(all(lambda > 0))
  if(loi %in% c("weibull","gll")){
    stopifnot(all(alpha > 0) | all(rho > 0) )
  }
  if(nom_loglik == "logcomplet"){
    stopifnot( length(lambda) == 3, length(alpha) == 3)
    if(is.null(rho)){
      parinit <- c(lambda, alpha,betavector)
      parametres_names <- c(paste0("lambda", suffixes),paste0("alpha", suffixes),
                            beta_names)
    }else{
      parinit <- c(lambda, alpha, rho, betavector)
      parametres_names <- c(paste0("lambda", suffixes),paste0("alpha", suffixes),
                            paste0("rho", suffixes),beta_names)
    }
  }
  if(nom_loglik == "logincomplet"){
    stopifnot(length(lambda) == 2, length(alpha) == 2, !is.null(beta))
    suffixes2 <- c("12", "13")
    if(is.null(rho)){
      parinit <- c(lambda, alpha, beta, betavector)
      parametres_names <- c(paste0("lambda", suffixes2),paste0("alpha", suffixes2), "prop", beta_names)
    }else{
      parinit <- c(lambda, alpha,rho, beta, betavector)
      parametres_names <- c(paste0("lambda", suffixes2),paste0("alpha", suffixes2),
                            paste0("rho", suffixes2),"prop", beta_names)
    }
  }
  npar <- length(parinit)

  if(nom_loglik == "logincomplet"){
    Low <- c(rep(0.00001, npar-length(betavector)-1),-Inf, rep(-Inf, length(betavector)))
  }else{
    Low <- c(rep(0.00001, npar-length(betavector)), rep(-Inf, length(betavector)))
  }


  if(loi == "gompertz"){
    npar2 <- (length(alpha)+length(lambda))%/%2
    Low <- c(rep(0.00001, npar2), rep(-Inf, npar-npar2))
  }

  Upe <- rep(Inf, npar)
  fn_name <- paste0(nom_loglik, "_", loi)
  # vérification d'existence
  if (!exists(fn_name, mode = "function")) {
    stop(sprintf("La loi '%s' n'existe pas. Fonction attendue : %s()",loi,fn_name))
  }
  # récupération de la fonction
  fn_objective <- get(fn_name, mode = "function")
  optim_log <- optim(
    par = parinit,
    fn = fn_objective,
    method = "L-BFGS-B",
    data = dataset,
    NomCovariable = NomCovariable,
    lower = Low,
    upper = Upe,
    control = list(maxit = 1000, ndeps = rep(1e-4, npar)),
    hessian = TRUE
  )
  cat("\n >>> Optimization diagnostics:\n")
  cat("Convergence status (0 = success): ", optim_log$convergence, "\n")
  cat("\n")
  parametres_estimes <- optim_log$par
  variance_estime <- tryCatch({
    diag(solve(optim_log$hessian))
  }, error = function(e) {
    return(NULL)
  })
  if (!is.null(variance_estime)) {
    varetst <- sqrt(variance_estime)
    tetst <- parametres_estimes / varetst
    zts <- qnorm(1-0.05/2)
    OptimResult <- data.frame(
      Estimate = parametres_estimes,
      Se = varetst,
      Z_stat   = tetst,
      ICL   = parametres_estimes - varetst*zts,
      ICU   = parametres_estimes + varetst*zts,
      row.names = parametres_names
    )

  } else {
    message("Calcul de la variance impossible (matrice singulière).")
    OptimResult <- data.frame(
      Estimate = parametres_estimes,
      Se = NA,
      Z_stat   = NA,
      row.names = parametres_names
    )
  }
  return( OptimResult)
}


