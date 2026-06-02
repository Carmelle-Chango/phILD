#' phildtestII
#'
#' Performs a likelihood ratio test between the constrained and unconstrained
#' illness--death models, considering the piecewise constant exponential model
#' and the piecewise linear approximation.
#'
#' @param Lambda Matrix of initial parameter values of dimension $(2 * 3)$.
#'
#' @param beta Initial value of the proportionality parameter.
#' This argument is used only for estimating the constrained
#' illness--death model.
#'
#' @param betax Matrix of initial values for covariate coefficients.
#'
#' @param distribution Name of the model used. Possible values are:
#' piecewise constant exponential (\code{"pwexp"}) and
#' piecewise linear approximation (\code{"pla"}).
#' The default model is \code{"pwexp"}.
#'
#' @param breakpoints Vector of cut points for the piecewise linear approximation
#' or matrix of intervals for the piecewise constant exponential model.
#'
#' @param dataset Dataset used for estimation.
#'
#' @param NomCovariable Vector containing the names of the covariates
#' included in the dataset.
#'
#' @return The function returns the likelihood ratio test statistic
#' along with the associated p-value.
#' @seealso phildtest
#'
#' @export
phildtestII <- function(Lambda,beta = NULL,betax = NULL,
                        distribution="pwexp", dataset = NULL,
                        breakpoints, NomCovariable= NULL){

  loi <- tolower(distribution)
  cat("\n", rep("=", 100), "\n", sep = "")
  if(loi == "pla"){
    cat(">> Step 1: Optimization of the restricted model (H0)...")
    incomplet <-  phildII1(Lambda= Lambda, betax = betax, beta=beta, distribution=loi,
                           nom_loglik = "logincomplet",dataset = dataset,
                           NomCovariable = NomCovariable, breakpoints = breakpoints)

    parInComplet <- incomplet[,1]
    fn_nameincomplet <- paste0("logincomplet", "_", loi)
    fn_nameincomplet <- get(fn_nameincomplet, mode = "function")
    H0G <- do.call(fn_nameincomplet, list(par =  parInComplet, data = dataset,
                                          NomCovariable = NomCovariable,
                                          breakpoints=breakpoints))

    cat(">> Step 2: Optimization of the full model (H1)...")
    complet <-   phildII1(Lambda= Lambda, betax = betax, distribution=loi,
                          nom_loglik = "logcomplet",dataset = dataset,
                          NomCovariable = NomCovariable, breakpoints = breakpoints)

    parComplet <- complet[,1]
    fn_namecomplet <- paste0("logcomplet", "_", loi)
    fn_namecomplet <- get(fn_namecomplet, mode = "function")
    H1G <- do.call(fn_namecomplet, list(par= parComplet, data = dataset,
                                        NomCovariable = NomCovariable,
                                        breakpoints=breakpoints))
  }else if(loi == "pwexp"){
    fn_name <- "totalLoglik"
    fn_name <- get(fn_name, mode = "function")

    cat(">> Step 1: Optimization of the restricted model (H0)...")

    incomplet <- phildII1(Lambda= Lambda, betax = betax, beta=beta, distribution= loi,
                          nom_loglik = "logincomplet", dataset = dataset,
                          NomCovariable = NomCovariable, breakpoints = breakpoints)

    parInComplet <- incomplet[,1]
    H0G <- do.call(fn_name, list(par =  parInComplet, data = dataset,nom_loglik = "logincomplet",
                                 NomCovariable = NomCovariable, nbeta = 2*nrow(betax) + 1,
                                 breakpoints=breakpoints))
    cat(">> Step 2: Optimization of the full model (H1)...")

    complet <- phildII1(Lambda= Lambda, betax = betax, distribution= loi,
                        nom_loglik = "logcomplet", dataset = dataset,
                        NomCovariable = NomCovariable, breakpoints = breakpoints)

    parComplet <- complet[,1]
    H1G <- do.call(fn_name, list(par= parComplet, data = dataset,nom_loglik = "logcomplet",
                                 NomCovariable = NomCovariable,nbeta= 3*nrow(betax),
                                 breakpoints=breakpoints))
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
}
#' @noRd
logincomplet_pwexp <- function(par, data, breakpoints, NomCovariable=NULL){
  q12 <- par[1]
  q13 <- par[2]
  beta <- par[3]
  data <- data.frame(data)
  if(!is.null(NomCovariable)){
    X <- as.matrix(data[, NomCovariable])
    ncov <- ncol(X)
    sequence <- seq(from = 4, length.out = 2*ncov)
    data$lin12 <- as.vector(X %*% par[sequence[1:ncov]])
    data$lin13 <- as.vector(X %*% par[sequence[(ncov+1):(2*ncov)]])
  }else{
    data$lin12 <- data$lin13 <- data$lin23 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*indi_vec(stop, breakpoints)*( log(q13) + lin13)
                                                           - q13*int1_vec(stop, breakpoints)*exp(lin13)  ))
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(status*indi_vec(stop, breakpoints)*(log(q12) + lin12 )
                                                           - q12*int1_vec(stop, breakpoints)*exp(lin12) ))
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(status*indi_vec(stop, breakpoints)*( log(q13*exp(beta))  + lin13 )
                                                           - q13*exp(beta)*int2_vec(start, stop, breakpoints)*exp(lin13) ))
  L <- L12 + L23 + L13
  return(-L)
}
#' @noRd
logincomplet_pla <- function(par, breakpoints, data , NomCovariable=NULL){
  if(!is.null(NomCovariable)){
    X <- as.matrix(data[,NomCovariable])
    ncov <- ncol(X)
  }
  ind <- length(breakpoints)-1
  indices_auto <- c(1:(ind-1), (ind-1))
  lenpar <- length(par)
  par12 <- par[1:ind]
  par13 <- par[(ind+1):(2*ind)]
  beta <- par[(2*ind+1)]
  c <- par[(2*ind+2):(lenpar)]
  ind_c <- length(c)
  groupes <- split(1:ind_c, ceiling(seq_along(1:ind_c)/ncov))
  par23 <- par13*exp(beta)
  b12 <- pmax(diff(par12)/diff(breakpoints[1:ind]), 0)
  b13 <- pmax(diff(par13)/diff(breakpoints[1:ind]),0)
  b23 <- pmax(diff(par23)/diff(breakpoints[1:ind]),0)
  a12 <- par12[-1] - b12*breakpoints[-c(1,ind+1)]
  a13 <- par13[-1] - b13*breakpoints[-c(1,ind+1)]
  a23 <- par23[-1] - b23*breakpoints[-c(1,ind+1)]
  rates12 <- cbind(a12, b12)[indices_auto,]
  rates13 <- cbind(a13, b13)[indices_auto,]
  rates23 <- cbind(a23, b23)[indices_auto,]
  data <- data.frame(data)
  if(!is.null(NomCovariable)){
    X <- as.matrix(data[,NomCovariable])
    ncov <- ncol(X)
    data$lin12 <- as.vector(X %*% c[groupes[[1]]])
    data$lin13 <- as.vector(X %*% c[groupes[[2]]])
  }else{
    data$lin12 <- data$lin13 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*(log(hpla(stop,breakpoints,rates13)) + lin13)
                                                           - Hpla(stop,breakpoints, rates13)*exp(lin13)))
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(status*(log(hpla(stop,breakpoints,rates12)) + lin12 )
                                                           - Hpla(stop,breakpoints, rates12)*exp(lin12 )) )
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(status*(log(hpla(stop,breakpoints,rates23)) + lin13 )
                                                           - HHpla(stop, start,breakpoints, rates23)*exp(lin13 )))
  L <- L12 + L23 + L13
  return(-L)
}
#' @noRd
logcomplet_pwexp <- function(par, data, breakpoints, NomCovariable=NULL){
  q12 <- par[1]
  q13 <- par[2]
  q23 <- par[3]
  data <- data.frame(data)
  if(!is.null(NomCovariable)){
    X <- as.matrix(data[, NomCovariable])
    ncov <- ncol(X)
    sequence <- seq(from = 4, length.out = 3*ncov)
    data$lin12 <- as.vector(X %*% par[sequence[1:ncov]])
    data$lin13 <- as.vector(X %*% par[sequence[(ncov+1):(2*ncov)]])
    data$lin23 <- as.vector(X %*% par[sequence[(2*ncov+1):(3*ncov)]])
  }else{
    data$lin12 <- data$lin13 <- data$lin23 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*indi_vec(stop, breakpoints)*( log(q13) + lin13)
                                                           - q13*int1_vec(stop, breakpoints)*exp(lin13)  ))
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(status*indi_vec(stop, breakpoints)*(log(q12) + lin12 )
                                                           - q12*int1_vec(stop, breakpoints)*exp(lin12) ))
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(status*indi_vec(stop, breakpoints)*( log(q23)  + lin23 )
                                                           - q23*int2_vec(start, stop, breakpoints)*exp(lin23) ))
  L <- L12 + L23 + L13

  return(-L)
}
#' @noRd
logcomplet_pla <- function(par, breakpoints, data = NULL, NomCovariable=NULL){
  if(!is.null(NomCovariable)){
    X <- as.matrix(data[,NomCovariable])
    ncov <- ncol(X)
  }
  ind <- length(breakpoints)-1
  indices_auto <- c(1:(ind-1), (ind-1))
  lenpar <- length(par)
  par12 <- par[1:ind]
  par13 <- par[(ind+1):(2*ind)]
  par23 <- par[(2*ind+1):(3*ind)]
  c <- par[(3*ind+1):(lenpar)]
  ind_c <- length(c)
  groupes <- split(1:ind_c, ceiling(seq_along(1:ind_c)/ncov))
  b12 <- pmax(diff(par12)/diff(breakpoints[1:ind]), 0)
  b13 <- pmax(diff(par13)/diff(breakpoints[1:ind]),0)
  b23 <- pmax(diff(par23)/diff(breakpoints[1:ind]),0)
  a12 <- par12[-1] - b12*breakpoints[-c(1,ind+1)]
  a13 <- par13[-1] - b13*breakpoints[-c(1,ind+1)]
  a23 <- par23[-1] - b23*breakpoints[-c(1,ind+1)]
  rates12 <- cbind(a12, b12)[indices_auto,]
  rates13 <- cbind(a13, b13)[indices_auto,]
  rates23 <- cbind(a23, b23)[indices_auto,]
  data <- data.frame(data)

  if(!is.null(NomCovariable)){
    X <- as.matrix(data[,NomCovariable])
    ncov <- ncol(X)
    data$lin12 <- as.vector(X %*% c[groupes[[1]]])
    data$lin13 <- as.vector(X %*% c[groupes[[2]]])
    data$lin23 <- as.vector(X %*% c[groupes[[3]]])
  }else{
    data$lin12 <- data$lin13 <- data$lin23 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*(log(hpla(stop,breakpoints,rates13)) + lin13)
                                                           - Hpla(stop,breakpoints, rates13)*exp(lin13)))
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(status*(log(hpla(stop,breakpoints,rates12)) + lin12 )
                                                           - Hpla(stop,breakpoints, rates12)*exp(lin12 )) )
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(status*(log(hpla(stop,breakpoints,rates23)) + lin23 )
                                                           - HHpla(stop, start,breakpoints, rates23)*exp(lin23 )))
  L <- L12 + L23 + L13
  return(-L)
}
#' @noRd
totalLoglik <- function(par, breakpoints, data = NULL, nom_loglik = "logcomplet",
                        NomCovariable = NULL, nbeta =NULL ){
  lambda_vec <- NULL
  n_L <- length(par) - nbeta
  lambda_vec <- par[1:n_L]
  n_int <- nrow(breakpoints)
  lambda_mat <- matrix(lambda_vec, nrow = n_int, byrow = FALSE)
  beta_uniques <- par[(n_L + 1):length(par)]
  if(nom_loglik == "logincomplet"){
    beta <- par[n_L+1]
    beta_uniques <- par[(n_L + 2):length(par)]
  }else if(nom_loglik == "logcomplet"){
    beta <- NULL
  }
  fn_name <- paste0(nom_loglik, "_pwexp")
  logliksum <- 0
  for(i in 1:n_int){
    loglik_i <- do.call(fn_name, list(par = c(lambda_mat[i,], beta, beta_uniques), breakpoints = breakpoints[i,],  data = data,
                                      NomCovariable = NomCovariable))
    logliksum <- logliksum + loglik_i
  }
  return(logliksum)
}
#' @noRd
phildII1 <- function(Lambda= NULL, betax = NULL, beta=NULL,distribution= "pla",
                     nom_loglik = "logcomplet", dataset = NULL,
                     NomCovariable = NULL, breakpoints = NULL){
  beta_names <- NULL
  loi <- tolower(distribution)
  if (nom_loglik == "logincomplet") {
    indices  <- 1:2
    suffixes <- c("12", "13")
  } else {
    indices  <- 1:ncol(Lambda)
    suffixes <- c("12", "13", "23")
  }
  Lambda_red <- Lambda[, indices, drop = FALSE]
  n_rows = nrow(betax)
  if (is.null(betax)) {
    if (nom_loglik == "logincomplet"){
      par_init <- c(as.vector(Lambda_red),beta)
      beta_names <- "beta"
    } else if (nom_loglik == "logcomplet"){
      par_init <- as.vector(Lambda_red)
      beta_names <- NULL
    }
  } else {
    betax_red <- betax[, indices, drop = FALSE]
    if (nom_loglik == "logincomplet"){
      par_init <- c(as.vector(Lambda_red),beta, as.vector(betax_red))
      beta_names <- c("beta",paste0("beta", rep(1:n_rows, length(suffixes)),
                                    rep(suffixes, each = n_rows)))
    } else if (nom_loglik == "logcomplet"){
      par_init <- c(as.vector(Lambda_red), as.vector(betax_red))

      beta_names <- paste0("beta", rep(1:n_rows, length(suffixes)),
                           rep(suffixes, each = n_rows))
    }
  }

  stopifnot(all(Lambda > 0))
  n_intervalles <- nrow(Lambda_red)
  n_transitions <- ncol(Lambda_red)
  n_covariables <- if (is.null(betax)) 0 else ncol(betax[, indices, drop=FALSE])*nrow(betax)
  n_trans <- length(suffixes)
  lambda_names <- paste0("lambda",
                         rep(suffixes, times = n_intervalles),
                         "_intervalle_",
                         rep(1:n_intervalles, each = n_trans))
  parametres_names <- c(lambda_names, beta_names)
  low_lambda <- rep(0.00001, n_intervalles * n_transitions)
  if (nom_loglik == "logincomplet"){
    low_beta   <- rep(-Inf,  n_covariables+1)
  } else if (nom_loglik == "logcomplet"){
    low_beta   <- rep(-Inf,  n_covariables)
  }
  Low <- c(low_lambda, low_beta)
  Upe <- rep(Inf, length(Low))
  fn_name <- switch (loi,
                     "pwexp" = "totalLoglik",
                     "pla" = paste0(nom_loglik, "_", loi),
                     stop("Loi non supportée")
  )
  if (!exists(fn_name, mode = "function")){
    stop(sprintf("La loi '%s' n'existe pas. Fonction attendue : %s()",loi,
                 fn_name))
  }
  fn_objective <- get(fn_name, mode = "function")
  nbeta = n_covariables
  if(nom_loglik == "logincomplet"){
    nbeta = nbeta + 1
  }

  arguments_specifiques <- switch (loi,
                                   "pwexp" = list(nom_loglik = nom_loglik, nbeta=nbeta),
                                   "pla" = list()
  )
  args_optim <- list(
    par     = par_init,
    fn      = fn_objective,
    method  = "L-BFGS-B",
    breakpoints = breakpoints,
    data = dataset,
    NomCovariable = NomCovariable,
    lower   = Low,
    upper   = Upe,
    control = list(maxit = 1000, ndeps = rep(1e-4, length(par_init))),
    hessian = TRUE
  )
  optim_log <- do.call(optim, c(args_optim, arguments_specifiques))
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
  return(OptimResult)
}



