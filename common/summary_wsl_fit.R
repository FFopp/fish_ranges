### =========================================================================
### define summary function for wsl.fit objects
### =========================================================================

#' @export
summary.wsl.fit=function(object){

  cat("Call: \n")
  print(object@call)

  cat("\nMeta information: \n")
  df=as.data.frame(object@meta[c("author","date","project")])

  if(as.character(object@call)[1]=="wsl.flex"){
    df$model_tags=paste(names(object@fits[[1]]),collapse=", ")
  } else {
    df$model_tags=object@meta["model_tag"]
  }

  rownames(df)=""
  print(df)

  cat("\nVariables used: \n")

  df=as.data.frame(object@meta[c("taxon","env_vars")])
  rownames(df)=""
  print(df)

  cat("\nResampling: \n")

  df=as.data.frame(object@meta[c("replicatetype","replicates")])
  rownames(df)=""
  print(df)

  cat("\nOther: \n")

  if(as.character(object@call)[1]=="wsl.glm"){
    df=as.data.frame(object@meta[c("step")])
    rownames(df)=""
    print(df)
  }


}
