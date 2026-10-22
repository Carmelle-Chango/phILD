#' phild
#'
#' Estimates the parameters of constrained and unconstrained illness--death models
#' using the maximum likelihood method, considering Weibull, Gompertz, and
#' generalized log-logistic distributions.
#'
#' @param lambda_int Vector of initial values for the scale parameters of the
#' baseline intensity functions for the considered distributions (Weibull,
#' Gompertz, and generalized log-logistic).
#'
#' @param alpha_int Vector of initial values for the shape parameters of the
#' considered distributions (Weibull, Gompertz, and generalized log-logistic).
#'
#' @param rho_int Vector of initial values for the location parameters of the
#' generalized log-logistic distribution.
#'
#' @param beta_int Initial value of the proportionality parameter.
#' This argument is used only for estimating the constrained
#' illness--death model.
#'
#' @param betax_int Matrix of initial values for covariate coefficients.
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
#' @param cov_effect Character string specifying the treatment of fixed
#' covariate effects in the constrained illness-death model. Use "common" (default)
#'  to impose equal covariate effects on transitions 1->3 and 2->3, or "specific"
#'  to estimate separate covariate effects for these transitions.
#'  This argument has no effect on the unconstrained model.
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
#'   lambda_int = rep(0.001, 3),
#'   alpha_int = rep(0.001, 3),
#'   betax_int = matrix(
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
#' @seealso phildFlex
#'
#' @export
phild <- function(lambda_int = NULL, alpha_int = NULL,
                  beta_int = NULL, rho_int = NULL,
                  betax_int = NULL,
                  distribution = "weibull",
                  nom_loglik = "logcomplet",
                  dataset = NULL,
                  NomCovariable = NULL,
                  cov_effect = "common") {

  cov_effect <- match.arg(cov_effect, c("common", "specific"))
  nom_loglik <- match.arg(nom_loglik,
                          c("logcomplet", "logincomplet"))
  loi <- match.arg(tolower(distribution),
                   c("weibull", "gompertz", "gll"))

  if (nom_loglik == "logcomplet" || cov_effect == "specific") {
    suffixes <- c("12", "13", "23")
  } else {
    suffixes <- c("12", "13")
  }

  beta_names <- character(0)
  betavector <- numeric(0)

  if (!is.null(NomCovariable)) {

    if (is.null(betax_int)) {
      stop("betax_int must be provided when NomCovariable is not NULL.")
    }

    betax_int <- as.matrix(betax_int)

    if (nrow(betax_int) != length(NomCovariable)) {
      stop(
        "The number of rows in betax_int must match ",
        "the number of covariates."
      )
    }

    if (nom_loglik == "logcomplet" || cov_effect == "specific") {

      if (ncol(betax_int) != 3) {
        stop("betax_int must contain three columns.")
      }

    } else {

      if (ncol(betax_int) < 2) {
        stop("betax_int must contain at least two columns.")
      }

      betax_int <- betax_int[, 1:2, drop = FALSE]
    }

    betavector <- as.vector(betax_int)

    beta_names <- paste0(
      "beta",
      rep(seq_len(nrow(betax_int)), length(suffixes)),
      rep(suffixes, each = nrow(betax_int))
    )

  } else if (!is.null(betax_int)) {

    stop("NomCovariable must be provided when betax_int is not NULL.")
  }

  ntrans <- if (nom_loglik == "logcomplet") 3 else 2

  if (loi %in% c("weibull", "gll")) {

    if (length(lambda_int) != ntrans ||
        length(alpha_int) != ntrans ||
        any(!is.finite(lambda_int)) ||
        any(!is.finite(alpha_int)) ||
        any(lambda_int <= 0) ||
        any(alpha_int <= 0)) {

      stop(
        "lambda_int and alpha_int must contain positive, ",
        "finite values and have the appropriate lengths."
      )
    }

  } else {

    # Gompertz: alpha_int = scale parameters
    #           rho_int   = shape parameters

    if (length(alpha_int) != ntrans ||
        length(rho_int) != ntrans ||
        any(!is.finite(alpha_int)) ||
        any(!is.finite(rho_int)) ||
        any(alpha_int <= 0)) {

      stop(
        "For the Gompertz distribution, alpha_int must be ",
        "strictly positive, and alpha_int and rho_int must ",
        "have the appropriate lengths and contain finite values."
      )
    }
  }

  if (loi == "gll") {

    if (length(rho_int) != ntrans ||
        any(!is.finite(rho_int)) ||
        any(rho_int <= 0)) {

      stop(
        "For the generalized log-logistic distribution, ",
        "rho_int must contain strictly positive, finite values ",
        "and have the appropriate length."
      )
    }
  }

  if (nom_loglik == "logincomplet") {

    if (length(beta_int) != 1 || !is.finite(beta_int)) {
      stop("beta_int must contain exactly one finite value.")
    }
  }

  suffixes_base <- if (nom_loglik == "logcomplet") {
    c("12", "13", "23")
  } else {
    c("12", "13")
  }

  if (loi == "weibull") {

    parinit <- c(lambda_int, alpha_int)

    parametres_names <- c(
      paste0("lambda", suffixes_base),
      paste0("alpha", suffixes_base)
    )

    Low_base <- rep(1e-5, length(parinit))

  } else if (loi == "gompertz") {

    parinit <- c(alpha_int, rho_int)

    parametres_names <- c(
      paste0("alpha", suffixes_base),
      paste0("rho", suffixes_base)
    )

    Low_base <- c(
      rep(1e-5, length(alpha_int)),
      rep(-Inf, length(rho_int))
    )

  } else {

    # Generalized log-logistic: lambda, alpha, rho

    parinit <- c(lambda_int, alpha_int, rho_int)

    parametres_names <- c(
      paste0("lambda", suffixes_base),
      paste0("alpha", suffixes_base),
      paste0("rho", suffixes_base)
    )

    Low_base <- rep(1e-5, length(parinit))
  }

  if (nom_loglik == "logincomplet") {

    parinit <- c(parinit, beta_int)
    parametres_names <- c(parametres_names, "prop")
    Low_base <- c(Low_base, -Inf)
  }

  parinit <- c(parinit, betavector)
  parametres_names <- c(parametres_names, beta_names)

  npar <- length(parinit)

  Low <- c(Low_base, rep(-Inf, length(betavector)))
  Upe <- rep(Inf, npar)

  # Select the appropriate likelihood function

  if (nom_loglik == "logincomplet" &&
      cov_effect == "specific") {

    fn_name <- paste0("logincomplet_spec_", loi)

  } else {

    fn_name <- paste0(nom_loglik, "_", loi)
  }

  if (!exists(fn_name, mode = "function", inherits = TRUE)) {

    stop(sprintf(
      "The likelihood function '%s' was not found.",
      fn_name
    ))
  }

  fn_objective <- get(
    fn_name,
    mode = "function",
    inherits = TRUE
  )

  # Maximum likelihood estimation

  optim_log <- optim(
    par = parinit,
    fn = fn_objective,
    method = "L-BFGS-B",
    data = dataset,
    NomCovariable = NomCovariable,
    lower = Low,
    upper = Upe,
    control = list(maxit = 1000),
    hessian = TRUE
  )

  cat("\n>>> Optimization diagnostics:\n")
  cat(
    "Convergence status (0 = success):",
    optim_log$convergence, "\n"
  )
  cat(
    "Maximized log-likelihood:",
    -optim_log$value, "\n"
  )
  cat(
    "Number of parameters:",
    npar, "\n\n"
  )

  parametres_estimes <- optim_log$par

  # Compute standard errors from the observed Hessian

  variance_estime <- tryCatch({

    diag(solve(optim_log$hessian))

  }, error = function(e) {

    NULL
  })

  if (!is.null(variance_estime) &&
      all(is.finite(variance_estime)) &&
      all(variance_estime > 0)) {

    varetst <- sqrt(variance_estime)
    tetst <- parametres_estimes / varetst
    zts <- qnorm(1 - 0.05 / 2)

    OptimResult <- data.frame(
      Estimate = parametres_estimes,
      Se = varetst,
      Z_stat = tetst,
      ICL = parametres_estimes - varetst * zts,
      ICU = parametres_estimes + varetst * zts,
      row.names = parametres_names
    )

  } else {

    warning("Unable to compute valid standard errors.")

    OptimResult <- data.frame(
      Estimate = parametres_estimes,
      Se = NA_real_,
      Z_stat = NA_real_,
      ICL = NA_real_,
      ICU = NA_real_,
      row.names = parametres_names
    )
  }

  # Organize parameter estimates by transition

  tableaux_separes <- lapply(suffixes, function(s) {

    indices <- grep(
      paste0(s, "$"),
      rownames(OptimResult)
    )

    sub_tab <- OptimResult[indices, , drop = FALSE]

    current_names <- rownames(sub_tab)
    is_beta <- grepl("^beta[0-9]+", current_names)

    if (any(is_beta)) {
      current_names[is_beta] <- paste0(NomCovariable, s)
    }

    clean_names <- gsub(
      paste0(s, "$"),
      "",
      current_names
    )

    rownames(sub_tab) <- clean_names

    sub_tab
  })

  names(tableaux_separes) <- paste0(
    "Transition ", suffixes
  )

  if ("prop" %in% rownames(OptimResult)) {

    tab_prop <- OptimResult["prop", , drop = FALSE]

    tableaux_separes[["Proportionality coefficient"]] <- tab_prop
  }

  return(tableaux_separes)
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
logincomplet_spec_weibull <- function(par, data, NomCovariable=NULL){
  l12 <- par[1]
  l13 <- par[2]
  a12 <- par[3]
  a13 <- par[4]
  beta<- par[5]
  data <- data.frame(data)
  if(!is.null(NomCovariable)){
    X <- as.matrix(data[,NomCovariable])
    ncov <- ncol(X)
    data$lin12 <- as.vector(X %*% par[5 + (1:ncov)])
    data$lin13 <- as.vector(X %*% par[5 + ((ncov+1):(2*ncov))])
    data$lin23 <- as.vector(X %*% par[5 + ((2*ncov+1):(3*ncov))])
  }else{
    data$lin12 <- data$lin13 <- data$lin23 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*(log(a13*l13*(stop)^(a13 - 1)) +  lin13)
                                                           - l13*((stop)^a13)*exp(lin13) ))
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(status*(log(a12*l12*(stop)^(a12-1)  ) + lin12)
                                                           - l12*((stop)^a12)*exp(lin12)))
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(status*(log(a13*l13*exp(beta)*(stop)^(a13 - 1) ) + lin23)
                                                           -l13*exp(beta)*(stop^a13 - start^a13 )*exp(lin23)))

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
logincomplet_spec_gompertz <- function(par,data, NomCovariable=NULL){
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
    data$lin23 <- as.vector(X %*% par[5 + ((2*ncov+1):(3*ncov))])
  }else{
    data$lin12 <- data$lin13 <- data$lin23 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*(log(a13*exp(r13*stop) ) + lin13 )
                                                           -(a13/r13)*(exp(r13*stop) - 1)*exp(lin13)) )
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(status*(log(a12*exp(r12*stop)) + lin12 )
                                                           -(a12/r12)*(exp(r12*stop) - 1)*exp(lin12)) )
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(status*(log(a13*exp(beta)*exp(r13*stop)) + lin23 )
                                                           - ((a13*exp(beta))/r13)*( exp(r13*stop) - exp(r13*start) )*exp(lin23) ))

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
logincomplet_spec_gll <- function(par,data, NomCovariable=NULL){
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
    data$lin23 <- as.vector(X %*% par[7 + ((2*ncov+1):(3*ncov))])
  }else{
    data$lin12 <- data$lin13 <- data$lin23 <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(status*( log(a13*l13*(stop)^(a13 -1))
                                                                    - log(1 + r13*(stop)^a13) + lin13)
                                                           - (l13/r13)*log(1 + r13*(stop)^a13)*exp(lin13)  ))
  L12 <- with(data[(data$from == 1)&(data$to ==2), ], sum(status*( log(a12*l12*(stop)^(a12 -1))
                                                                   - log(1 + r12*(stop)^a12) + lin12)
                                                          - (l12/r12)*log(1 + r12*(stop)^a12)*exp(lin12)  ))
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum( status*( log(a13*l13*exp(beta)*(stop)^(a13 -1))
                                                                     - log(1 + r13*(stop)^a13) + lin23)
                                                            - (l13*exp(beta)/r13)* ( log(1 + r13*(stop)^a13)
                                                                                     - log(1 + r13*(start)^a13))*exp(lin23)  ))
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
phild1 <- function(lambda_int = NULL, alpha_int = NULL,
                  beta_int = NULL, rho_int = NULL,
                  betax_int = NULL,
                  distribution = "weibull",
                  nom_loglik = "logcomplet",
                  dataset = NULL,
                  NomCovariable = NULL,
                  cov_effect = "common") {

  cov_effect <- match.arg(cov_effect, c("common", "specific"))
  nom_loglik <- match.arg(nom_loglik,
                          c("logcomplet", "logincomplet"))
  loi <- match.arg(tolower(distribution),
                   c("weibull", "gompertz", "gll"))

  if (nom_loglik == "logcomplet" || cov_effect == "specific") {
    suffixes <- c("12", "13", "23")
  } else {
    suffixes <- c("12", "13")
  }

  beta_names <- character(0)
  betavector <- numeric(0)

  if (!is.null(NomCovariable)) {

    if (is.null(betax_int)) {
      stop("betax_int must be provided when NomCovariable is not NULL.")
    }

    betax_int <- as.matrix(betax_int)

    if (nrow(betax_int) != length(NomCovariable)) {
      stop(
        "The number of rows in betax_int must match ",
        "the number of covariates."
      )
    }

    if (nom_loglik == "logcomplet" || cov_effect == "specific") {

      if (ncol(betax_int) != 3) {
        stop("betax_int must contain three columns.")
      }

    } else {

      if (ncol(betax_int) < 2) {
        stop("betax_int must contain at least two columns.")
      }

      betax_int <- betax_int[, 1:2, drop = FALSE]
    }

    betavector <- as.vector(betax_int)

    beta_names <- paste0(
      "beta",
      rep(seq_len(nrow(betax_int)), length(suffixes)),
      rep(suffixes, each = nrow(betax_int))
    )

  } else if (!is.null(betax_int)) {

    stop("NomCovariable must be provided when betax_int is not NULL.")
  }

  ntrans <- if (nom_loglik == "logcomplet") 3 else 2

  if (loi %in% c("weibull", "gll")) {

    if (length(lambda_int) != ntrans ||
        length(alpha_int) != ntrans ||
        any(!is.finite(lambda_int)) ||
        any(!is.finite(alpha_int)) ||
        any(lambda_int <= 0) ||
        any(alpha_int <= 0)) {

      stop(
        "lambda_int and alpha_int must contain positive, ",
        "finite values and have the appropriate lengths."
      )
    }

  } else {

    # Gompertz: alpha_int = scale parameters
    #           rho_int   = shape parameters

    if (length(alpha_int) != ntrans ||
        length(rho_int) != ntrans ||
        any(!is.finite(alpha_int)) ||
        any(!is.finite(rho_int)) ||
        any(alpha_int <= 0)) {

      stop(
        "For the Gompertz distribution, alpha_int must be ",
        "strictly positive, and alpha_int and rho_int must ",
        "have the appropriate lengths and contain finite values."
      )
    }
  }

  if (loi == "gll") {

    if (length(rho_int) != ntrans ||
        any(!is.finite(rho_int)) ||
        any(rho_int <= 0)) {

      stop(
        "For the generalized log-logistic distribution, ",
        "rho_int must contain strictly positive, finite values ",
        "and have the appropriate length."
      )
    }
  }

  if (nom_loglik == "logincomplet") {

    if (length(beta_int) != 1 || !is.finite(beta_int)) {
      stop("beta_int must contain exactly one finite value.")
    }
  }

  suffixes_base <- if (nom_loglik == "logcomplet") {
    c("12", "13", "23")
  } else {
    c("12", "13")
  }

  if (loi == "weibull") {

    parinit <- c(lambda_int, alpha_int)

    parametres_names <- c(
      paste0("lambda", suffixes_base),
      paste0("alpha", suffixes_base)
    )

    Low_base <- rep(1e-5, length(parinit))

  } else if (loi == "gompertz") {

    parinit <- c(alpha_int, rho_int)

    parametres_names <- c(
      paste0("alpha", suffixes_base),
      paste0("rho", suffixes_base)
    )

    Low_base <- c(
      rep(1e-5, length(alpha_int)),
      rep(-Inf, length(rho_int))
    )

  } else {

    # Generalized log-logistic: lambda, alpha, rho

    parinit <- c(lambda_int, alpha_int, rho_int)

    parametres_names <- c(
      paste0("lambda", suffixes_base),
      paste0("alpha", suffixes_base),
      paste0("rho", suffixes_base)
    )

    Low_base <- rep(1e-5, length(parinit))
  }

  if (nom_loglik == "logincomplet") {

    parinit <- c(parinit, beta_int)
    parametres_names <- c(parametres_names, "prop")
    Low_base <- c(Low_base, -Inf)
  }

  parinit <- c(parinit, betavector)
  parametres_names <- c(parametres_names, beta_names)

  npar <- length(parinit)

  Low <- c(Low_base, rep(-Inf, length(betavector)))
  Upe <- rep(Inf, npar)

  # Select the appropriate likelihood function

  if (nom_loglik == "logincomplet" &&
      cov_effect == "specific") {

    fn_name <- paste0("logincomplet_spec_", loi)

  } else {

    fn_name <- paste0(nom_loglik, "_", loi)
  }

  if (!exists(fn_name, mode = "function", inherits = TRUE)) {

    stop(sprintf(
      "The likelihood function '%s' was not found.",
      fn_name
    ))
  }

  fn_objective <- get(
    fn_name,
    mode = "function",
    inherits = TRUE
  )

  # Maximum likelihood estimation

  optim_log <- optim(
    par = parinit,
    fn = fn_objective,
    method = "L-BFGS-B",
    data = dataset,
    NomCovariable = NomCovariable,
    lower = Low,
    upper = Upe,
    control = list(maxit = 1000),
    hessian = TRUE
  )

  cat("\n>>> Optimization diagnostics:\n")
  cat(
    "Convergence status (0 = success):",
    optim_log$convergence, "\n"
  )
  cat(
    "Maximized log-likelihood:",
    -optim_log$value, "\n"
  )
  cat(
    "Number of parameters:",
    npar, "\n\n"
  )

  parametres_estimes <- optim_log$par

  # Compute standard errors from the observed Hessian

  variance_estime <- tryCatch({

    diag(solve(optim_log$hessian))

  }, error = function(e) {

    NULL
  })

  if (!is.null(variance_estime) &&
      all(is.finite(variance_estime)) &&
      all(variance_estime > 0)) {

    varetst <- sqrt(variance_estime)
    tetst <- parametres_estimes / varetst
    zts <- qnorm(1 - 0.05 / 2)

    OptimResult <- data.frame(
      Estimate = parametres_estimes,
      Se = varetst,
      Z_stat = tetst,
      ICL = parametres_estimes - varetst * zts,
      ICU = parametres_estimes + varetst * zts,
      row.names = parametres_names
    )

  } else {

    warning("Unable to compute valid standard errors.")

    OptimResult <- data.frame(
      Estimate = parametres_estimes,
      Se = NA_real_,
      Z_stat = NA_real_,
      ICL = NA_real_,
      ICU = NA_real_,
      row.names = parametres_names
    )
  }

  return(OptimResult)
}
