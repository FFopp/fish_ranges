### =========================================================================
### prediction-evaluation prediction function
### =========================================================================
#' Correctly feed the predict functions depending on model type (glm, gbm, maxent...)
#'
#' Not to be called directly by the user
#' @author Philipp
#' @export
prd=function(mod,tst){

  # Generate probabilistic precitions
  if("MaxEnt"%in%class(mod)){

    pred<-df_or_rast(mod,tst)

  } else if("glm"%in%class(mod)){

    pred<-df_or_rast(mod=mod,
                     nwdat=tst,
                     type="response")

  } else if("gbm"%in%class(mod)){

    pred<-df_or_rast(mod,
                     nwdat=tst,
                     n.trees=mod$n.trees,
                     type="response")

  } else if("randomForest"%in%class(mod)){

    pred<-df_or_rast(mod,
                     nwdat=tst,
                     type="prob")
  }

  # Convert to numeric
  if(class(tst)=="data.frame"){
    pred<-as.numeric(pred)
  }

  return(pred)

}
