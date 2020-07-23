### =========================================================================
### Set Class
### =========================================================================
#' An S4 class to store pseudoabsence data
#'
#' @slot meta A list with meta information
#' @slot type The type of pseudoabsence sampling
#' @slot pa a vector with presences and pseudo absences
#' @slot env_vars data.frame with environmental predictors
#' @slot xy a matrix with coordinates of the points
#' @author Philipp
#' @export
wsl.pseudoabsences<-setClass("wsl.pseudoabsences",slots=c(meta="list", # Meta information
                                                          pa="numeric", # store presence/pseudoabsence information
                                                          env_vars="data.frame", # store extracted env variables
                                                          xy="matrix", # store coordiantes
                                                          call="call")) # conserve function call