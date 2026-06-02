#' transILD
#'
#' Transforms a dataset into a long format compatible with the package functions.
#' The likelihood functions implemented in the package are constructed under the
#' assumption that the data are organized in long format. This function also allows
#' any dataset, whether initially in short or long format, to be reformatted
#' according to the structure required by the package functions. We therefore
#' recommend that users first apply the \code{transILD} function before passing
#' a dataset to the functions \code{phild}, \code{phildII},
#' \code{phildtest}, and \code{phildtestII}.
#'
#' @param initial_data Dataset in short or long format.
#'
#' @param identifiant Name of the column corresponding to the individual
#' identifier.
#'
#' @param observation Name of the column corresponding to the event
#' indicator.
#'
#' @param indicateur Name of the column corresponding to the illness
#' indicator.
#'
#' @param tempsUn Name of the column corresponding to the entry time
#' into state $k$.
#'
#' @param tempsDeux Name of the column corresponding to the exit time
#' from state $k$.
#'
#' @param format_long Indicates whether the dataset is already in
#' long format (\code{TRUE}) or in short format (\code{FALSE}).
#'
#' @return The function returns a dataset containing the variables:
#' \code{id} (individual identifier), \code{start} (entry time into state $k$),
#' \code{stop} (exit time from state $k$), \code{status} (indicates whether the
#' transition $k -> l$ is observed), \code{from} (starting state), \code{to}
#' (destination state), as well as the covariates present in the initial dataset.
#'
#' @examples
#' library(survival)
#' data(heart)
#'
#' transILD(
#'   heart,
#'   identifiant = "id",
#'   observation = "event",
#'   indicateur = "transplant",
#'   tempsUn = "start",
#'   tempsDeux = "stop",
#'   format_long = FALSE
#' )
#'
#' @export
#'
#' @import dplyr
transILD <- function(initial_data, identifiant = NULL, observation = NULL,
                     indicateur = NULL, tempsUn = NULL,
                     tempsDeux = NULL, format_long){

  initial_data <- data.frame(initial_data)
  if(isTRUE(format_long)){
    nouvelle_base <- initial_data
    nom_cible <- c("id", "from", "to","status","start","stop")
    for(nom_final in nom_cible){
      cat("Enter the name of the column in your dataset corresponding to ", nom_final," : ", sep = " ")
      nouveau_nom <- scan(what = character(),nlines = 1, quiet = TRUE)
      teste <- nouveau_nom %in% colnames(initial_data)
      conteur <- 0
      while(!teste){
        conteur <- conteur + 1
        if(conteur > 3){
          stop("Please restart after more than 3 errors")
        }
        cat(nouveau_nom , "was not found in your dataset")
        cat("Please re-enter the name of the column in your dataset corresponding to ", nom_final," : ", sep = " ")
        nouveau_nom <- scan(what = character(), nlines = 1, quiet = TRUE)
        teste <- nouveau_nom %in% colnames(initial_data)
      }
      colnames( nouvelle_base)[colnames( nouvelle_base) == nouveau_nom] <- nom_final
    }
    transformed_data <- nouvelle_base
  }else{
    transformed_data <- initial_data %>%
      group_by(id) %>%
      do({
        id_data <- .
        autres_colonnes <- id_data[1, !names(id_data) %in% c(identifiant,observation,indicateur,tempsUn,tempsDeux)]
        if(nrow(id_data) > 1) {
          data1 =  data.frame(
            id = id_data[[identifiant]][1],
            from = c(1, 1, 2),
            to = c(3, 2, 3),
            status = c(0, 1, id_data[[observation]][2]),
            start = c(0,0,id_data[[tempsUn]][2]),
            stop = c(id_data[[tempsDeux]][1],id_data[[tempsDeux]][1],id_data[[tempsDeux]][2])
          )
        } else {
          data1 =   data.frame(
            id = id_data[[identifiant]][1],
            from = c(1, 1),
            to = c(2, 3),
            status = c(0, id_data[[observation]][1]),
            start = c(0,id_data[[tempsUn]][1]),
            stop = c(id_data[[tempsDeux]][1],id_data[[tempsDeux]][1])
          )
        }
        cbind(data1,autres_colonnes)
      }) %>%
      ungroup()
  }
  transformed_data = data.frame(transformed_data)
  return( transformed_data )
}
