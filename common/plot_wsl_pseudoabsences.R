### =========================================================================
### define summary function for wsl.evaluation objects
### =========================================================================

#' @export
plot.wsl.pseudoabsences=function(object){

  xy_pres=object@xy[which(object@pa==1),]

  xy_abs=object@xy[which(object@pa==0),]
  if(nrow(xy_abs)>10000){
    xy_abs=xy_abs[sample(1:nrow(xy_abs),10000),]
  }

  if(!is.na(object@meta$template_file)){
    rst=raster(object@meta$template_file)
    plot(rst,col=c("#f0f0f0","#99d8c9"),
         main=object@meta$type,legend=F)

  } else {
    plot(xy_abs,
         pch=16,
         xlab="",
         ylab="",
         main=object@meta$type,
         type="n")
  }

  points(xy_abs,pch=16,cex=.2,col="#00000020")
  points(xy_pres,pch=3,cex=.5,col="darkred")

  legend("bottomleft",
         pch=c(16,3),
         col=c("#00000020","darkred"),
         c(paste0("absences (n=",length(which(object@pa==0)),")"),
           paste0("presences (n=",length(which(object@pa==1)),")")))

}
