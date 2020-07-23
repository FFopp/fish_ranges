### =========================================================================
### define summary function for wsl.evaluation objects
### =========================================================================

#' @export
summary.wsl.evaluation=function(object){

  cat("\nMeta information: \n")
  df=data.frame(object@meta[c("author","date")],object@meta$wsl.fit[c("project","replicatetype","replicates")])

  rownames(df)=""
  print(df)

  cat("\nThreshold: \n")
  df=as.data.frame(object@meta[c("cutoff")])

  rownames(df)=""
  print(df)

  cat("\nMean skill: \n")

  mats=list()
  for(i in 1:length(object@performance)){
    mats[[i]]=do.call("cbind",object@performance[[i]])
  }
  mn=Reduce("+", mats) / length(mats)

  print(mn)

}
