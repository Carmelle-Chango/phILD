#'phildtestFlex
#'
#'  Likelihood Ratio Test for Flexible Illness-Death Models
#'
#' Compares constrained and unconstrained illness-death models
#' with PWEXP or PLA baseline hazards.
#'
#' @param lambda_int Initial baseline matrix (intervals x 3 transitions).
#' @param betax_int Initial covariate coefficient matrix (covariates x 3).
#' @param beta_int Initial proportionality coefficient for H0.
#' @param distribution Either "pwexp" or "pla".
#' @param dataset Long-format dataset.
#' @param breakpoints PWEXP interval matrix or PLA cutpoint vector.
#' @param NomCovariable Character vector of covariate names, or NULL.
#' @param cov_effect "common" or "specific" under H0.
#' @param weights NULL, numeric observation weights, or a weight-column name.
#' @return An object of class "htest".
#' @seealso phildFlex
#' @export
phildtestFlex <- function(lambda_int = NULL,
                          betax_int = NULL,
                          beta_int = NULL,
                          distribution = "pwexp",
                          dataset = NULL,
                          breakpoints = NULL,
                          NomCovariable = NULL,
                          cov_effect = "common",
                          weights = NULL) {

  loi <- match.arg(tolower(distribution), c("pwexp", "pla"))
  cov_effect <- match.arg(cov_effect, c("common", "specific"))

  if (is.null(dataset)) stop("dataset must be provided.")
  data_name <- deparse(substitute(dataset))
  dataset <- as.data.frame(dataset)
  weights <- get_weights(dataset, weights)

  p <- if (is.null(NomCovariable)) 0L else length(NomCovariable)
  if (p > 0L && is.null(betax_int)) {
    stop("betax_int must be provided when covariates are specified.")
  }
  if (p == 0L && !is.null(betax_int)) {
    stop("NomCovariable must be provided when betax_int is specified.")
  }

  cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
  cat(">> Step 1: Optimization of the restricted model (H0)...\n")

  incomplet <- phildII1(
    lambda_int = lambda_int,
    betax_int = betax_int,
    beta_int = beta_int,
    distribution = loi,
    nom_loglik = "logincomplet",
    dataset = dataset,
    NomCovariable = NomCovariable,
    breakpoints = breakpoints,
    cov_effect = cov_effect,
    weights = weights
  )

  parInComplet <- incomplet[, "Estimate"]

  cat(">> Step 2: Optimization of the full model (H1)...\n")

  complet <- phildII1(
    lambda_int = lambda_int,
    betax_int = betax_int,
    distribution = loi,
    nom_loglik = "logcomplet",
    dataset = dataset,
    NomCovariable = NomCovariable,
    breakpoints = breakpoints,
    cov_effect = cov_effect,
    weights = weights
  )

  parComplet <- complet[, "Estimate"]

  if (any(!is.finite(parInComplet)) || any(!is.finite(parComplet))) {
    stop("Non-finite parameter estimates: check optimization.")
  }

  cat(">> Step 3: Likelihood Ratio Test calculation...\n")

  if (loi == "pwexp") {
    nbeta_H0 <- 1L + if (cov_effect == "common") 2L * p else 3L * p
    nbeta_H1 <- 3L * p

    H0G <- totalLoglik(
      par = parInComplet,
      breakpoints = breakpoints,
      data = dataset,
      nom_loglik = "logincomplet",
      NomCovariable = NomCovariable,
      nbeta = nbeta_H0,
      cov_effect = cov_effect,
      weights = weights
    )

    H1G <- totalLoglik(
      par = parComplet,
      breakpoints = breakpoints,
      data = dataset,
      nom_loglik = "logcomplet",
      NomCovariable = NomCovariable,
      nbeta = nbeta_H1,
      cov_effect = cov_effect,
      weights = weights
    )
  } else {
    fn_H0 <- if (cov_effect == "specific") {
      logincomplet_spec_pla
    } else {
      logincomplet_pla
    }

    H0G <- fn_H0(
      par = parInComplet,
      breakpoints = breakpoints,
      data = dataset,
      NomCovariable = NomCovariable,
      weights = weights
    )

    H1G <- logcomplet_pla(
      par = parComplet,
      breakpoints = breakpoints,
      data = dataset,
      NomCovariable = NomCovariable,
      weights = weights
    )
  }

  if (length(H0G) != 1L || length(H1G) != 1L ||
      !all(is.finite(c(H0G, H1G)))) {
    stop("The optimized negative log-likelihoods must be finite scalars.")
  }

  ddl <- length(parComplet) - length(parInComplet)
  ddl_attendu <- if (cov_effect == "common") 1L + p else 1L
  if (ddl <= 0L || ddl != ddl_attendu) {
    stop(sprintf("Unexpected degrees of freedom: got %d; expected %d.",
                 ddl, ddl_attendu))
  }

  DG <- 2 * (H0G - H1G)
  if (DG < -1e-6) {
    warning("H0 has a higher optimized log-likelihood than H1. ",
            "Check optimization and parameter ordering.")
  }
  # Do not silently transform a negative LRT into zero.
  pvalG <- if (DG >= 0) {
    pchisq(DG, df = ddl, lower.tail = FALSE)
  } else {
    NA_real_
  }

  cat("Negative log-likelihood (H0): ", H0G, "\n", sep = "")
  cat("Negative log-likelihood (H1): ", H1G, "\n", sep = "")
  cat("Number of parameters (H0): ", length(parInComplet), "\n", sep = "")
  cat("Number of parameters (H1): ", length(parComplet), "\n", sep = "")
  cat("LRT statistic: ", DG, "\n", sep = "")
  cat("Degrees of freedom: ", ddl, "\n", sep = "")
  cat("P-value: ", pvalG, "\n", sep = "")
  cat("\n>>> Test completed.\n")
  cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")

  structure(
    list(
      statistic = c(LRT = DG),
      p.value = as.numeric(pvalG),
      parameter = c(df = ddl),
      method = paste0("Likelihood Ratio Test for Proportional Hazards ",
                      "assumption (", toupper(loi),
                      "; cov_effect = ", cov_effect, ")"),
      data.name = data_name
    ),
    class = "htest"
  )
}


#' @noRd
get_weights <- function(data, weights = NULL) {

  if (is.null(weights)) {
    return(rep(1, nrow(data)))
  }

  if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(data)) {
      stop("Weight column not found in dataset.")
    }
    weights <- data[[weights]]
  }

  if (!is.numeric(weights) ||
      length(weights) != nrow(data) ||
      anyNA(weights) ||
      any(!is.finite(weights)) ||
      any(weights < 0)) {
    stop("weights must be numeric, finite, non-negative, ",
         "and have one value per row of dataset.")
  }

  as.numeric(weights)
}

#' @noRd
logincomplet_pwexp <- function(par, data, breakpoints, NomCovariable=NULL ,weights = NULL){
  q12 <- par[1]
  q13 <- par[2]
  beta <- par[3]

  data <- data.frame(data)
  data$weig <- get_weights(data, weights)

  if(!is.null(NomCovariable)){
    X <- as.matrix(data[, NomCovariable])
    ncov <- ncol(X)
    sequence <- seq(from = 4, length.out = 2*ncov)
    data$lin12 <- as.vector(X %*% par[sequence[1:ncov]])
    data$lin13 <- as.vector(X %*% par[sequence[(ncov+1):(2*ncov)]])
  }else{
    data$lin12 <- data$lin13  <- 0
  }
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(weig*(status*indi_vec(stop, breakpoints)*( log(q13) + lin13)
                                                                 - q13*int2_vec(start, stop, breakpoints)*exp(lin13)  )) )
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(weig*(status*indi_vec(stop, breakpoints)*(log(q12) + lin12 )
                                                                 - q12*int2_vec(start, stop, breakpoints)*exp(lin12) )) )
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(weig*(status*indi_vec(stop, breakpoints)*( log(q13*exp(beta))  + lin13 )
                                                                 - q13*exp(beta)*int2_vec(start, stop, breakpoints)*exp(lin13) )) )
  L <- L12 + L23 + L13
  return(-L)
}
#' @noRd

logincomplet_spec_pwexp <- function(par, data, breakpoints,NomCovariable = NULL,weights = NULL) {

  q12 <- par[1]
  q13 <- par[2]
  beta <- par[3]

  data <- data.frame(data)
  data$weig <- get_weights(data, weights)

  if (!is.null(NomCovariable)) {
    X <- as.matrix(data[, NomCovariable, drop = FALSE])
    ncov <- ncol(X)
    sequence <- seq(from = 4, length.out = 3 * ncov)

    data$lin12 <- as.vector(X %*% par[sequence[1:ncov]])
    data$lin13 <- as.vector(X %*% par[sequence[(ncov + 1):(2 * ncov)]])
    data$lin23 <- as.vector(X %*% par[sequence[(2 * ncov + 1):(3 * ncov)]])

  } else {
    data$lin12 <- data$lin13 <- data$lin23 <- 0
  }

  L13 <- with(data[(data$from == 1) & (data$to == 3), ], sum(weig*( status * indi_vec(stop, breakpoints) *(log(q13) + lin13)
                                                                    - q13 * int2_vec(start, stop, breakpoints) * exp(lin13) )) )

  L12 <- with( data[(data$from == 1) & (data$to == 2), ], sum(weig*(status * indi_vec(stop, breakpoints) * (log(q12) + lin12) -
                                                                      q12 * int2_vec(start, stop, breakpoints) * exp(lin12))) )

  L23 <- with( data[(data$from == 2) & (data$to == 3), ],
               sum(weig*( status * indi_vec(stop, breakpoints) * (log(q13 * exp(beta)) + lin23) -
                            q13 * exp(beta) * int2_vec(start, stop, breakpoints) * exp(lin23) ) ) )

  L <- L12 + L23 + L13

  return(-L)
}
#' @noRd
logincomplet_pla <- function(par, breakpoints, data,
                             NomCovariable = NULL ,weights = NULL) {

  ind <- length(breakpoints) - 1
  lenpar <- length(par)

  par12 <- par[1:ind]
  par13 <- par[(ind + 1):(2 * ind)]
  beta <- par[2 * ind + 1]
  par23 <- par13 * exp(beta)

  # Pentes sur les intervalles finis
  b12 <- diff(par12) / diff(breakpoints[1:ind])
  b13 <- diff(par13) / diff(breakpoints[1:ind])
  b23 <- diff(par23) / diff(breakpoints[1:ind])

  # Ordonnées à l'origine sur les intervalles finis
  a12 <- par12[-1] - b12 * breakpoints[-c(1, ind + 1)]
  a13 <- par13[-1] - b13 * breakpoints[-c(1, ind + 1)]
  a23 <- par23[-1] - b23 * breakpoints[-c(1, ind + 1)]

  # Dernier intervalle : risque constant égal au dernier psi
  rates12 <- rbind(cbind(a12, b12), c(par12[ind], 0))
  rates13 <- rbind(cbind(a13, b13), c(par13[ind], 0))
  rates23 <- rbind(cbind(a23, b23), c(par23[ind], 0))


  data <- data.frame(data)
  data$weig <- get_weights(data, weights)


  if (!is.null(NomCovariable)) {

    X <- as.matrix( data[, NomCovariable, drop = FALSE] )
    ncov <- ncol(X)

    c <- par[(2 * ind + 2):lenpar]

    if (length(c) != 2 * ncov) {
      stop( "The common-effects constrained PLA model requires ",
            "two coefficients per covariate." )
    }

    data$lin12 <- as.vector(X %*% c[seq_len(ncov)])
    data$lin13 <- as.vector( X %*% c[ncov + seq_len(ncov)])

  } else {
    data$lin12 <- data$lin13 <- 0
  }

  L13 <- with(data[(data$from == 1) & (data$to == 3), ],sum(weig*(status * ( log(hpla(stop, breakpoints, rates13)) + lin13 ) -
                                                                    ( Hpla(stop, breakpoints, rates13) - Hpla(start, breakpoints, rates13) ) * exp(lin13))) )

  L12 <- with( data[(data$from == 1) & (data$to == 2), ], sum(weig*( status * (log(hpla(stop, breakpoints, rates12)) + lin12 ) -
                                                                       (Hpla(stop, breakpoints, rates12) - Hpla(start, breakpoints, rates12)) * exp(lin12))) )

  L23 <- with(data[(data$from == 2) & (data$to == 3), ], sum(weig*( status * ( log(hpla(stop, breakpoints, rates23)) + lin13 ) -
                                                                      ( Hpla(stop, breakpoints, rates23) -  Hpla(start, breakpoints, rates23) ) * exp(lin13))) )

  L <- L12 + L23 + L13

  return(-L)
}
#' @noRd
logincomplet_spec_pla <- function(par, breakpoints, data,
                                  NomCovariable = NULL ,weights = NULL) {

  ind <- length(breakpoints) - 1
  lenpar <- length(par)

  par12 <- par[1:ind]
  par13 <- par[(ind + 1):(2 * ind)]
  beta <- par[2 * ind + 1]
  par23 <- par13 * exp(beta)

  # Pentes sur les intervalles finis
  b12 <- diff(par12) / diff(breakpoints[1:ind])
  b13 <- diff(par13) / diff(breakpoints[1:ind])
  b23 <- diff(par23) / diff(breakpoints[1:ind])

  # Ordonnées à l'origine sur les intervalles finis
  a12 <- par12[-1] - b12 * breakpoints[-c(1, ind + 1)]
  a13 <- par13[-1] - b13 * breakpoints[-c(1, ind + 1)]
  a23 <- par23[-1] - b23 * breakpoints[-c(1, ind + 1)]

  # Dernier intervalle : risque constant égal au dernier psi
  rates12 <- rbind(cbind(a12, b12), c(par12[ind], 0))
  rates13 <- rbind(cbind(a13, b13), c(par13[ind], 0))
  rates23 <- rbind(cbind(a23, b23), c(par23[ind], 0))

  data <- data.frame(data)
  data$weig <- get_weights(data, weights)

  if (!is.null(NomCovariable)) {

    X <- as.matrix( data[, NomCovariable, drop = FALSE] )
    ncov <- ncol(X)

    c <- par[(2 * ind + 2):lenpar]

    if (length(c) != 3 * ncov) {
      stop( "The specific-effects constrained PLA model requires ",
            "three coefficients per covariate." )
    }

    data$lin12 <- as.vector(X %*% c[seq_len(ncov)])
    data$lin13 <- as.vector( X %*% c[ncov + seq_len(ncov)])
    data$lin23 <- as.vector(X %*% c[2 * ncov + seq_len(ncov)])

  } else {
    data$lin12 <- data$lin13 <-  data$lin23 <- 0
  }

  L13 <- with(data[(data$from == 1) & (data$to == 3), ],sum(weig*(status * ( log(hpla(stop, breakpoints, rates13)) + lin13 ) -
                                                                    ( Hpla(stop, breakpoints, rates13) - Hpla(start, breakpoints, rates13) ) * exp(lin13))) )

  L12 <- with( data[(data$from == 1) & (data$to == 2), ], sum(weig*( status * (log(hpla(stop, breakpoints, rates12)) + lin12 ) -
                                                                       (Hpla(stop, breakpoints, rates12) - Hpla(start, breakpoints, rates12)) * exp(lin12))) )

  L23 <- with(data[(data$from == 2) & (data$to == 3), ], sum(weig*( status * ( log(hpla(stop, breakpoints, rates23)) + lin23 ) -
                                                                      ( Hpla(stop, breakpoints, rates23) -  Hpla(start, breakpoints, rates23) ) * exp(lin23))) )

  L <- L12 + L23 + L13

  return(-L)
}
#' @noRd
logcomplet_pwexp <- function(par, data, breakpoints, NomCovariable=NULL ,weights = NULL){
  q12 <- par[1]
  q13 <- par[2]
  q23 <- par[3]
  data <- data.frame(data)
  data$weig <- get_weights(data, weights)
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
  L13 <-  with(data[(data$from == 1)&(data$to ==3), ], sum(weig*(status*indi_vec(stop, breakpoints)*( log(q13) + lin13)
                                                                 - q13*int2_vec(start, stop, breakpoints)*exp(lin13)  )) )
  L12 <-  with(data[(data$from == 1)&(data$to ==2), ], sum(weig*(status*indi_vec(stop, breakpoints)*(log(q12) + lin12 )
                                                                 - q12*int2_vec(start, stop, breakpoints)*exp(lin12) )) )
  L23 <-  with(data[(data$from == 2)&(data$to ==3), ], sum(weig*(status*indi_vec(stop, breakpoints)*( log(q23)  + lin23 )
                                                                 - q23*int2_vec(start, stop, breakpoints)*exp(lin23) )) )
  L <- L12 + L23 + L13

  return(-L)
}
#' @noRd

logcomplet_pla <- function(par, breakpoints, data = NULL,
                           NomCovariable = NULL ,weights = NULL) {

  ind <- length(breakpoints) - 1

  lenpar <- length(par)

  par12 <- par[1:ind]
  par13 <- par[(ind + 1):(2 * ind)]
  par23 <- par[(2 * ind + 1):(3 * ind)]

  # Pentes sur les intervalles finis
  b12 <- diff(par12) / diff(breakpoints[1:ind])
  b13 <- diff(par13) / diff(breakpoints[1:ind])
  b23 <- diff(par23) / diff(breakpoints[1:ind])

  # Ordonnées à l'origine sur les intervalles finis
  a12 <- par12[-1] - b12 * breakpoints[-c(1, ind + 1)]
  a13 <- par13[-1] - b13 * breakpoints[-c(1, ind + 1)]
  a23 <- par23[-1] - b23 * breakpoints[-c(1, ind + 1)]

  # Dernier intervalle : risque constant égal au dernier psi
  rates12 <- rbind(cbind(a12, b12), c(par12[ind], 0))
  rates13 <- rbind(cbind(a13, b13), c(par13[ind], 0))
  rates23 <- rbind(cbind(a23, b23), c(par23[ind], 0))

  data <- data.frame(data)
  data$weig <- get_weights(data, weights)

  if (!is.null(NomCovariable)) {

    X <- as.matrix(data[, NomCovariable, drop = FALSE])
    ncov <- ncol(X)

    c <- par[(3 * ind + 1):lenpar]

    if (length(c) != 3 * ncov) {
      stop(
        "The unconstrained PLA model requires ",
        "three coefficients per covariate."
      )
    }

    data$lin12 <- as.vector( X %*% c[seq_len(ncov)] )
    data$lin13 <- as.vector( X %*% c[ncov + seq_len(ncov)] )
    data$lin23 <- as.vector( X %*% c[2 * ncov + seq_len(ncov)])

  } else {

    data$lin12 <- data$lin13 <- data$lin23 <- 0
  }


  L13 <- with(data[(data$from == 1) & (data$to == 3), ],sum(weig*(status * ( log(hpla(stop, breakpoints, rates13)) + lin13 ) -
                                                                    ( Hpla(stop, breakpoints, rates13) - Hpla(start, breakpoints, rates13) ) * exp(lin13))) )

  L12 <- with( data[(data$from == 1) & (data$to == 2), ], sum(weig*( status * (log(hpla(stop, breakpoints, rates12)) + lin12 ) -
                                                                       (Hpla(stop, breakpoints, rates12) - Hpla(start, breakpoints, rates12)) * exp(lin12))) )

  L23 <- with(data[(data$from == 2) & (data$to == 3), ], sum(weig*( status * ( log(hpla(stop, breakpoints, rates23)) + lin23 ) -
                                                                      ( Hpla(stop, breakpoints, rates23) -  Hpla(start, breakpoints, rates23) ) * exp(lin23))) )


  L <- L12 + L23 + L13

  return(-L)
}

#' @noRd

totalLoglik <- function(par, breakpoints, data = NULL,
                        nom_loglik = "logcomplet",
                        NomCovariable = NULL,
                        nbeta = NULL,
                        cov_effect = "common" ,weights = NULL) {


  nom_loglik <- match.arg( nom_loglik, c("logcomplet", "logincomplet"))

  cov_effect <- match.arg(cov_effect, c("common", "specific"))

  n_int <- nrow(breakpoints)
  n_L <- length(par) - nbeta
  lambda_vec <- par[seq_len(n_L)]
  n_trans <- if (nom_loglik == "logcomplet") 3 else 2

  if (n_L != n_int * n_trans) {
    stop(
      "The number of baseline intensity parameters ",
      "does not match the number of intervals and transitions."
    )
  }
  lambda_mat <- matrix(
    lambda_vec,
    nrow = n_int,
    byrow = FALSE
  )

  if (nom_loglik == "logincomplet") {

    beta <- par[n_L + 1]
    beta_uniques <- if (length(par) > n_L + 1) {
      par[(n_L + 2):length(par)]
    } else {
      numeric(0)
    }

  } else {

    beta <- NULL
    beta_uniques <- if (length(par) > n_L) {
      par[(n_L + 1):length(par)]
    } else {
      numeric(0)
    }
  }
  if (nom_loglik == "logincomplet" &&
      cov_effect == "specific") {

    fn_name <- "logincomplet_spec_pwexp"

  } else {

    fn_name <- paste0(nom_loglik, "_pwexp")
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

  logliksum <- 0

  for (i in seq_len(n_int)) {

    loglik_i <- do.call(
      fn_objective,
      list(
        par = c(lambda_mat[i, ], beta, beta_uniques),
        breakpoints = breakpoints[i, ],
        data = data,
        NomCovariable = NomCovariable,
        weights = weights
      )
    )

    logliksum <- logliksum + loglik_i
  }

  return(logliksum)
}


phildII1 <- function(lambda_int = NULL,
                     betax_int = NULL,
                     beta_int = NULL,
                     distribution = "pwexp",
                     nom_loglik = "logcomplet",
                     dataset = NULL,
                     NomCovariable = NULL,
                     breakpoints = NULL,
                     cov_effect = "common",
                     weights = NULL) {

  loi <- tolower(distribution)

  cov_effect <- match.arg(
    cov_effect,
    choices = c("common", "specific")
  )

  if (!nom_loglik %in% c("logincomplet", "logcomplet")) {
    stop("nom_loglik must be 'logincomplet' or 'logcomplet'.")
  }

  if (!loi %in% c("pwexp", "pla")) {
    stop("Unsupported distribution. Use 'pwexp' or 'pla'.")
  }

  if (is.null(dataset)) {
    stop("dataset must be provided.")
  }

  dataset <- as.data.frame(dataset)
  weights <- get_weights(dataset, weights)

  if (is.null(lambda_int) || !is.matrix(lambda_int)) {
    stop("lambda_int must be a matrix.")
  }

  if (ncol(lambda_int) != 3) {
    stop(
      "lambda_int must contain three columns ",
      "corresponding to transitions 12, 13 and 23."
    )
  }

  if (any(!is.finite(lambda_int)) ||
      any(lambda_int <= 0)) {
    stop(
      "All initial baseline parameters must be ",
      "positive and finite."
    )
  }

  if (is.null(breakpoints)) {
    stop("breakpoints must be provided.")
  }

  # Validate breakpoints
  if (loi == "pwexp") {

    if (!is.matrix(breakpoints) ||
        ncol(breakpoints) != 2) {
      stop(
        "For PWEXP, breakpoints must be a two-column matrix."
      )
    }

    if (nrow(lambda_int) != nrow(breakpoints)) {
      stop(
        "For PWEXP, the number of rows in lambda_int ",
        "must match the number of intervals."
      )
    }

  } else {

    if (!is.numeric(breakpoints) ||
        !is.null(dim(breakpoints))) {
      stop(
        "For PLA, breakpoints must be a numeric vector."
      )
    }

    if (length(breakpoints) < 3 ||
        anyNA(breakpoints) ||
        any(!is.finite(utils::head(breakpoints, -1))) ||
        !is.infinite(utils::tail(breakpoints, 1)) ||
        tail(breakpoints, 1) < 0 ||
        any(diff(breakpoints) <= 0)) {

      stop(
        "For PLA, breakpoints must contain at least ",
        "two finite cut points followed by Inf, ",
        "in strictly increasing order."
      )
    }

    if (nrow(lambda_int) != length(breakpoints) - 1) {
      stop(
        "For PLA, nrow(lambda_int) must equal ",
        "length(breakpoints) - 1."
      )
    }
  }

  # Select baseline transitions
  if (nom_loglik == "logincomplet") {

    if (is.null(beta_int) ||
        length(beta_int) != 1 ||
        !is.finite(beta_int)) {
      stop(
        "A single finite beta_int value is required ",
        "for the constrained model."
      )
    }

    indices <- 1:2
    suffixes <- c("12", "13")

  } else {

    indices <- 1:3
    suffixes <- c("12", "13", "23")
  }

  # Initial baseline parameters
  lambda_red <- lambda_int[, indices, drop = FALSE]

  n_intervalles <- nrow(lambda_red)
  n_transitions <- ncol(lambda_red)

  par_init <- as.vector(lambda_red)

  lambda_names <- paste0(
    "lambda",
    rep(suffixes, each = n_intervalles),
    "_intervalle_",
    rep(seq_len(n_intervalles),
        times = n_transitions)
  )

  parametres_names <- lambda_names

  # Initial proportionality coefficient
  if (nom_loglik == "logincomplet") {

    par_init <- c(par_init, beta_int)
    parametres_names <- c(parametres_names, "beta")
  }

  # Initial covariate coefficients
  n_covariables <- 0
  suffixes_result <- suffixes

  if (!is.null(betax_int)) {

    if (is.null(NomCovariable)) {
      stop(
        "NomCovariable must be provided ",
        "when betax_int is not NULL."
      )
    }

    betax_int <- as.matrix(betax_int)
    n_rows <- nrow(betax_int)

    if (length(NomCovariable) != n_rows) {
      stop(
        "The number of rows in betax_int ",
        "must match the number of covariates."
      )
    }

    if (ncol(betax_int) != 3) {
      stop(
        "betax_int must contain three columns ",
        "corresponding to transitions 12, 13 and 23."
      )
    }

    if (any(!is.finite(betax_int))) {
      stop(
        "All initial covariate coefficients must be finite."
      )
    }

    # Select covariate effects
    if (nom_loglik == "logincomplet" &&
        cov_effect == "common") {

      beta_indices <- 1:2

    } else {

      beta_indices <- 1:3
    }

    betax_red <- betax_int[
      , beta_indices, drop = FALSE
    ]

    beta_suffixes <- c("12", "13", "23")[beta_indices]

    par_init <- c(
      par_init,
      as.vector(betax_red)
    )

    beta_names <- paste0(
      "beta",
      rep(seq_len(n_rows),
          times = length(beta_suffixes)),
      rep(beta_suffixes, each = n_rows)
    )

    parametres_names <- c(
      parametres_names,
      beta_names
    )

    n_covariables <- length(betax_red)

    # Display transition 23 when its covariate effects
    # are estimated separately
    if (nom_loglik == "logincomplet" &&
        cov_effect == "specific") {

      suffixes_result <- c("12", "13", "23")
    }

  } else if (!is.null(NomCovariable)) {

    stop(
      "betax_int must be provided when ",
      "NomCovariable is specified."
    )
  }

  # Parameter bounds
  low_lambda <- rep(
    1e-5,
    n_intervalles * n_transitions
  )

  n_other <- length(par_init) - length(low_lambda)

  low_beta <- rep(-Inf, n_other)

  Low <- c(low_lambda, low_beta)
  Upe <- rep(Inf, length(par_init))

  # Select likelihood function
  fn_name <- switch(
    loi,

    "pwexp" = "totalLoglik",

    "pla" = {
      if (nom_loglik == "logincomplet" &&
          cov_effect == "specific") {

        "logincomplet_spec_pla"

      } else {

        paste0(nom_loglik, "_pla")
      }
    }
  )

  if (!exists(fn_name, mode = "function")) {
    stop(
      sprintf(
        "Expected likelihood function '%s' was not found.",
        fn_name
      )
    )
  }

  fn_objective <- get(fn_name, mode = "function")

  nbeta <- n_covariables

  if (nom_loglik == "logincomplet") {
    nbeta <- nbeta + 1
  }

  arguments_specifiques <- switch(
    loi,

    "pwexp" = list(
      nom_loglik = nom_loglik,
      nbeta = nbeta,
      cov_effect = cov_effect
    ),

    "pla" = list()
  )

  # Maximum likelihood optimization
  args_optim <- list(
    par = par_init,
    fn = fn_objective,
    method = "L-BFGS-B",
    breakpoints = breakpoints,
    data = dataset,
    NomCovariable = NomCovariable,
    weights = weights,
    lower = Low,
    upper = Upe,
    control = list(maxit = 1000),
    hessian = TRUE
  )

  optim_log <- do.call(
    optim,
    c(args_optim, arguments_specifiques)
  )

  cat("\n>>> Optimization diagnostics:\n")
  cat(
    "Convergence status (0 = success): ",
    optim_log$convergence,
    "\n\n"
  )

  parametres_estimes <- optim_log$par

  # Variance estimation
  variance_estime <- tryCatch(
    diag(solve(optim_log$hessian)),
    error = function(e) NULL
  )

  if (!is.null(variance_estime) &&
      all(is.finite(variance_estime)) &&
      all(variance_estime >= 0)) {

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

    message(
      "Cannot calculate variance ",
      "(singular or invalid Hessian)."
    )

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

