#' phild
#'
#' Estimates the parameters of constrained and unconstrained illness--death models
#' using the maximum likelihood method, considering Weibull, Gompertz, and
#' generalized log-logistic distributions.
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
#' @param nom_loglik Specifies the illness--death model to be estimated:
#' \code{"logcomplet"} for the unconstrained model and
#' \code{"logincomplet"} for the constrained model.
#' The default is the unconstrained model.
#'
#' @param dataset Dataset used for estimation.
#'
#' @param NomCovariable Vector containing the names of the covariates
#' included in the dataset.
#'
#' @return The function returns the maximum likelihood estimators of the
#' transition intensity functions, as well as their variances and
#' 95% confidence intervals.
#' @examples
#' library(survival)
#' library(eha)
#'
#' data(heart)
#'
#' dataset <- transILD(
#'   heart,
#'   identifiant = "id",
#'   observation = "event",
#'   indicateur = "transplant",
#'   tempsUn = "start",
#'   tempsDeux = "stop",
#'   format_long = FALSE
#' )
#'
#' M <- phild(
#'   lambda = rep(0.001, 3),
#'   alpha = rep(0.001, 3),
#'   betax = matrix(
#'     c(rep(-0.001, 3),
#'       rep(-0.001, 3),
#'       rep(-0.001, 3)),
#'     3, 3,
#'     byrow = TRUE
#'   ),
#'   distribution = "weibull",
#'   nom_loglik = "logcomplet",
#'   dataset = dataset,
#'   NomCovariable = c("age", "year", "surgery")
#' )
#'
#' M
#'
#' phreg(
#'   Surv(start, stop, status) ~ age + year + surgery,
#'   data = dataset[(dataset$from == 1 & dataset$to == 3), ],
#'   dist = "weibull"
#' )
#'
#' M$`Transition 13`
#'
#' @seealso phildII
#'
#' @export
phild <- function(lambda= NULL, alpha = NULL, beta = NULL, rho = NULL, betax = NULL,
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
    stop(sprintf("The '%s' law does not exist. Expected function : %s()",loi,fn_name))
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
    message("Cannot calculate variance (singular matrix).")
    OptimResult <- data.frame(
      Estimate = parametres_estimes,
      Se = NA,
      Z_stat   = NA,
      row.names = parametres_names
    )
  }
  tableaux_separés <- lapply(suffixes, function(s) {
    indices <- grep(paste0(s, "$"), rownames(OptimResult))
    sub_tab <- OptimResult[indices, , drop = FALSE]
    current_names <- rownames(sub_tab)
    is_beta <- grepl("^beta", current_names)
    if (any(is_beta)) {
      current_names[is_beta] <- paste0(NomCovariable, s)
    }
    clean_names <- gsub(paste0(s, "$"), "", current_names)
    rownames(sub_tab) <- clean_names
    sub_tab
  })

  names(tableaux_separés) <- paste0("Transition ", suffixes)
  if ("prop" %in% rownames(OptimResult)) {
    tab_prop <- OptimResult["prop", , drop = FALSE]
    rownames(tab_prop) <- "prop"
    tableaux_separés[["Proportionnal coefficient"]] <- tab_prop
  }
  return(tableaux_separés)
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


