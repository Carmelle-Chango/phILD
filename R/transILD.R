#' transILD
#'
#' Transforme une base de données sous un format long compatible avec
#' les fonctions du package. Les fonctions de vraisemblance du package
#' sont construites en supposant que les données sont organisées sous
#' un format long. Cette fonction permet également de reformater toute base de données,
#' qu'elle soit initialement sous format court ou long, selon la structure
#' requise par les fonctions du package. Nous recommandons donc aux utilisateurs
#' d'utiliser préalablement la fonction \code{transILD} avant de passer un jeu
#' de données aux fonctions \code{phild}, \code{phildII},
#' \code{phildtest} et \code{phildtestII}.
#'
#' @param initial_data Base de données sous format court ou long.
#'
#' @param identifiant Nom de la colonne correspondant à l'identifiant
#' de chaque individu.
#'
#' @param observation Nom de la colonne correspondant à l'indicateur
#' d'événement.
#'
#' @param indicateur Nom de la colonne correspondant à l'indicateur
#' de maladie.
#'
#' @param tempsUn Nom de la colonne correspondant au temps d'entrée
#' dans un état $k$.
#'
#' @param tempsDeux Nom de la colonne correspondant au temps de sortie
#' d'un état $k$.
#'
#' @param format_long Indique si la base de données est déjà sous
#' format long (\code{TRUE}) ou sous format court (\code{FALSE}).
#'
#' @return La fonction retourne une base de données contenant les variables :
#' \code{id} (identifiant de l'individu), \code{start} (temps d'entrée dans l'état $k$),
#' \code{stop} (temps de sortie de l'état $k$), \code{status} (indique si la transition
#' $k -> l$ est observée), \code{from} (état de départ), \code{to} (état d'arrivée),
#' ainsi que les covariables présentes dans la base de données initiale.
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
      cat("Saisir le nom de la colonne de votre base correspondant a ", nom_final," : ", sep = " ")
      nouveau_nom <- scan(what = character(),nlines = 1, quiet = TRUE)
      teste <- nouveau_nom %in% colnames(initial_data)
      conteur <- 0
      while(!teste){
        conteur <- conteur + 1
        if(conteur > 3){
          stop("Veuillez reprendre à plus de 3 erreurs")
        }
        cat(nouveau_nom , "est introuvable dans votre base de donnée")
        cat("Veuillez saisir de nouveau  le nom de la colonne de votre base correspondant a ", nom_final," : ", sep = " ")
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


