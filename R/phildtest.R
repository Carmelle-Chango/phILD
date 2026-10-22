
#' phildtest
#'
#' Performs a likelihood ratio test between the constrained and unconstrained
#' illness--death models, considering Weibull, Gompertz, and generalized
#' log-logistic distributions.
#'
#' @param lambda_int Vector of initial values for the scale parameters of the
#' baseline intensity functions for the considered distributions (Weibull
#' and generalized log-logistic).
#'
#' @param alpha_int Vector of initial values for the shape parameters of the
#' considered distributions.
#'
#' @param rho_int Vector of initial values for the additional parameters
#' of the Gompertz and generalized log-logistic distributions.
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
#' @param dataset Dataset used for estimation.
#'
#' @param NomCovariable Vector containing the names of the covariates
#' included in the dataset.
#'
#' @param cov_effect Character string specifying the treatment of fixed
#' covariate effects in the constrained illness-death model. Use
#' \code{"common"} (default) to impose equal covariate effects on
#' transitions 1->3 and 2->3, or \code{"specific"} to estimate separate
#' covariate effects for these transitions. This argument has no effect
#' on the unconstrained model.
#'
#' @return An object of class \code{"htest"} containing the likelihood
#' ratio test statistic, degrees of freedom, and associated p-value.
#'
#' @seealso phildtestFlex
#' @export

phildtest <- function(lambda_int = NULL,
                      alpha_int = NULL,
                      rho_int = NULL,
                      beta_int = NULL,
                      betax_int = NULL,
                      distribution = "weibull",
                      dataset = NULL,
                      NomCovariable = NULL,
                      cov_effect = "common") {

  cov_effect <- match.arg(
    cov_effect, c("common", "specific")
  )

  loi <- match.arg(
    tolower(distribution),
    c("weibull", "gompertz", "gll")
  )

  cat("\n", rep("=", 100), "\n", sep = "")

  # Step 1: Fit the constrained model (H0)
  cat(">> Step 1: Optimization of the restricted model (H0)...\n")

  incomplet <- phild1(
    lambda_int = if (is.null(lambda_int)) NULL else lambda_int[1:2],
    alpha_int = if (is.null(alpha_int)) NULL else alpha_int[1:2],
    rho_int = if (is.null(rho_int)) NULL else rho_int[1:2],
    beta_int = beta_int,
    betax_int = betax_int,
    distribution = loi,
    nom_loglik = "logincomplet",
    dataset = dataset,
    NomCovariable = NomCovariable,
    cov_effect = cov_effect
  )

  parInComplet <- incomplet[, "Estimate"]

  if (cov_effect == "specific") {

    fn_nameincomplet <- paste0(
      "logincomplet_spec_", loi
    )

  } else {

    fn_nameincomplet <- paste0(
      "logincomplet_", loi
    )
  }

  fn_incomplet <- get(
    fn_nameincomplet,
    mode = "function"
  )

  H0G <- do.call(
    fn_incomplet,
    list(
      par = parInComplet,
      data = dataset,
      NomCovariable = NomCovariable
    )
  )

  cat("Negative log-likelihood (H0):", H0G, "\n")

  # Step 2: Fit the unconstrained model (H1)
  cat(">> Step 2: Optimization of the full model (H1)...\n")

  complet <- phild1(
    lambda_int = lambda_int,
    alpha_int = alpha_int,
    rho_int = rho_int,
    beta_int = beta_int,
    betax_int = betax_int,
    distribution = loi,
    nom_loglik = "logcomplet",
    dataset = dataset,
    NomCovariable = NomCovariable,
    cov_effect = cov_effect
  )

  parComplet <- complet[, "Estimate"]

  fn_namecomplet <- paste0(
    "logcomplet_", loi
  )

  fn_complet <- get(
    fn_namecomplet,
    mode = "function"
  )

  H1G <- do.call(
    fn_complet,
    list(
      par = parComplet,
      data = dataset,
      NomCovariable = NomCovariable
    )
  )

  cat("Negative log-likelihood (H1):", H1G, "\n")

  # Step 3: Compute the likelihood ratio test
  cat(">> Step 3: Likelihood Ratio Test calculation...\n")

  ddl <- nrow(complet) - nrow(incomplet)

  DG <- 2 * (H0G - H1G)

  if (!is.finite(DG)) {
    stop("The likelihood ratio test statistic is not finite.")
  }

  if (ddl <= 0) {
    stop("The degrees of freedom must be positive.")
  }

  if (DG < 0) {
    warning(
      "Negative likelihood ratio test statistic: ",
      "check optimization convergence and the maximized likelihoods."
    )
  }

  pvalG <- if (DG >= 0) {
    pchisq(DG, df = ddl, lower.tail = FALSE)
  } else {
    NA_real_
  }

  cat("LRT statistic:", DG, "\n")
  cat("Degrees of freedom:", ddl, "\n")
  cat("P-value:", pvalG, "\n")

  cat("\n>>> Test completed.\n")
  cat("\n", rep("=", 100), "\n", sep = "")

  structure(
    list(
      statistic = c(LRT = DG),
      p.value = as.numeric(pvalG),
      parameter = c(df = ddl),
      method = paste0(
        "Likelihood Ratio Test for Proportional Hazards assumption (",
        loi, "; cov_effect = ", cov_effect, ")"
      ),
      data.name = deparse(substitute(dataset))
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
