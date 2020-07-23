### =========================================================================
### Set Class
### =========================================================================
#' An S4 class to store prediction data
#'
#' @slot meta A list with meta information
#' @slot thres A vector with externally supplied thresholds
#' @slot predictions A list with model predictions
#' @author Philipp
#'
wsl.prediction<-setClass("wsl.prediction",slots=c(meta="list", # Meta information
                                                  thres="numeric", # supply external threshold
                                                  predictions="list")) # conserve function call
