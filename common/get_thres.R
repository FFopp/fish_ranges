### =========================================================================
### pull out mean thresholds from evaluation
### =========================================================================

#' Get threshold
#'
#' Extracts thresholds from wsl.evaluation objects and names them so they can be fed
#' to the wsl.predict function. At the moment only averages over replicates can be obtained.
#'
#' @param x a object of class wsl.evaluation
#' @return a vector with threshols
#' @author Philipp
#' #'
#' @export
get_thres=function(x){

  thrs=lapply(x@performance,function(y){
    a=sapply(y,function(z){
      return(z["threshold"])
    })
    return(a)
  })

  mat=do.call("rbind",thrs)

  out=colMeans(mat)

  names(out)=gsub(".threshold","",names(out))

  return(out)
}
