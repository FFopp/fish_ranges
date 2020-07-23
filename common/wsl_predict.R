### =========================================================================
### define wsl.predict function
### =========================================================================
#' Make predictions
#'
#' Make predictions with all models from a wsl.fit object. If thresholds are supplied
#' binary predictions are made, otherwise continuous predictions are returned.
#'
#' @param x an object of class wsl.fit ...
#' @param predat data.frame or raster for which predictions should be made
#' @param thres vector of the same length as number of models chosen with custom
#' thresholds for model evaluation. for wsl.flex outputs the thresholds have to
#' be labelled with the same names provided to models
#' @return Object of class wsl.prediction with slots for meta info, and model predictions
#' @author Philipp
#' @export
#'
wsl.predict<-function(x,predat=data.frame(),thres=numeric(),write=FALSE,output_dir=""){

  ### ------------------------
  ### Check input and prepare
  ### ------------------------

  if(nrow(predat)==0){
    stop("Prediction data missing!")
  }

  # thres has to be a vector with named elements (same names
  # as in evaluation matrix)
  if(length(thres)>0){
    if(length(x@fits[[1]])!=length(thres)){
      stop("Wrong number of thresholds supplied! Should be one threshold per model type...")
    }
  }

  if(write){

    # Create folder with species name in output directory
    spnm=gsub(" ","_",x@meta$taxon)
    dr=paste0(output_dir,"/",spnm)
    if(!dir.exists(dr)){
      dir.create(dr, recursive=T)
    }

    # Delete any binary predictions that are already in the folder
    oldfl=list.files(dr,pattern="BinaryPredictions",full.names = T)
    if(length(oldfl)>0){
      unlink(oldfl)
    }

  }

  ### ------------------------
  ### generate wsl.evaluation and add meta info
  ### ------------------------

  out<-preva.meta(type="prediction")

  ### -------------------------------------------
  ### Evaluate models
  ### -------------------------------------------

  lis<-list()
  # loop over replicates
  for(i in 1:length(x@fits)){

    lisa<-list()
    # Loop over model types
    for(j in 1:length(x@fits[[1]])){

      pred=prd(x@fits[[i]][[j]],predat)

      # Convert to binary predictions if thresholds were supplied
      if(length(thres)>0){
        the.tre=thres[which(names(thres)==names(x@fits[[i]])[j])]
      }

      if(write){

        flnm=paste0(dr,"/BinaryPredictions_",spnm,
                    "_PA_",x@meta$pseudoabsence_typ,
                    "_Mod_",names(x@fits[[1]])[j],".tif")
        if(exists("the.tre")){
          reclassify(pred,matrix(c(0,the.tre,0,the.tre,1,1),ncol=3,byrow=T),
                          filename = flnm,datatype="INT1U")
        } else {
          writeRaster(pred,filename = flnm)
        }

        rm(pred)
        # Return file path
        lisa[[j]]<-flnm
      } else {
        if(exists("the.tre")){
          if(class(pred)=="numeric"){
            pred[pred<the.tre]=0
            pred[pred>=the.tre]=1
          } else {
            pred=reclassify(pred,matrix(c(0,the.tre,0,the.tre,1,1),ncol=3,byrow=T))
          }

        }
        lisa[[j]]<-pred
      }

    }
    names(lisa)=names(x@fits[[i]])
    lis[[i]]<-lisa
  }
  names(lis)<-names(x@fits)
  out@predictions<-lis

  return(out)
}

