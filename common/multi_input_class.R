### =========================================================================
### define multi.input class
### =========================================================================
#' An S4 class to store evaluation data
#'
#' @slot mod Name of model algorithm to be called (character) e.g., 'glm'
#' @slot args A list of arguments to be supplied to the model fitting alrogrithm
#' @slot tag The name of the model setup (character)
#' @slot step A logical indicating whether step function should be ran (for glm
#' and gam)
#' @slot weight Should observations be weighted based on their prevalence?
#' @author Philipp
#' @export
multi.input<-setClass("multi.input",slots=c(mod="character", # Model function
                                            args="list", # Model function arguments
                                            tag="character", # Model set-up name
                                            step="logical", # Should step function be added?
                                            weight="logical")) # Should observations be weighted?

