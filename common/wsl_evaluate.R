#' Evaluate models
#'
#' Assess several model skill metrics for all models in a wsl.fit object. Currently
#' AUC, RMSE, TSS, PPV, Accuracy, and Cohen's Kappa are evaluated. Furthermore, the
#' threshold applied is returned.
#'
#' @param x a wsl.fit object
#' @param tester data.frame with testing data (only mandatory if replicatetype='none'
#' was chosen when models were fitted)
#' @param threshold vector of the same length as number of models chosen with custom
#' thresholds for model evaluation. for wsl.flex outputs the thresholds have to be labelled
#' with the same names provided to models
#' @param crit which threshold criterion should be considered? Currently 'pp=op'
#' (predicted prevalence = observed prevalence), 'maxTSS' (threshold yielding maximum TSS),
#' and 'external' (thresholds manually supplied) are possible
#' @param pseuAbsCor logical. Should a pseudoabsence correction be performed on the test
#' data. In a pseudoabsence correction the same amount of pseudoabsences as presences are
#' subsampled for model validation.
#' @return an obejct of class 'wsl.evaluation'
#' @author Philipp
#' @export
wsl.evaluate<-function(x,tester=data.frame(),thres=numeric(),crit="pp=op",pseuAbsCor=FALSE){

  ### ------------------------
  ### check tresholds
  ### ------------------------

  # thres has to be a vector with named elements (same names
  # as in evaluation matrix)
  if(length(thres)>0){
    if(length(x@fits[[1]])!=length(thres)){
      stop("Wrong number of thresholds supplied! Should be one threshold per model type...")
    }
    if(crit!="external"){
      warning("Assuming you want external tresholds to be used - setting crit='external'!")
      crit="external"
    }
  }

  if(!(crit%in%c("pp=op","maxTSS","external"))){
    stop("Invalid threshold criterion chosen!")
  }


  ### ------------------------
  ### Check testing data and prepare for evaluation
  ### ------------------------

  if(x@meta$replicatetype=="none" && nrow(tester)==0){
    stop("External testing data must be supplied for replicatetype 'none'")
  } else if(x@meta$replicatetype%in%c("cv","block-cv","splitsample")) {
    
    if(pseuAbsCor){
      x@tesdat=lapply(x@tesdat,function(y){
        tdpres=y[which(y$Presence==1),]
        tdabs=y[which(y$Presence==0),]
        tdabs=tdabs[sample(1:nrow(tdabs),nrow(tdpres)),]
        return(rbind(tdpres,tdabs))
      })
    }
    
    outerloop<-length(x@tesdat)
    testa<-lapply(x@tesdat,function(x){
      y<-x[,-which(colnames(x)=="Presence"),drop=FALSE]
    })
    papa<-lapply(x@tesdat,function(x){
      y<-x[,"Presence"]
    })

  } else if(x@meta$replicatetype=="none"){

    outerloop<-1
    testa<-list(tester[,-which(colnames(tester)=="Presence"),drop=FALSE])
    papa<-list(tester[,"Presence"])

  }

  ### ------------------------
  ### generate wsl.evaluation and add meta info
  ### ------------------------

  out<-preva.meta(type="evaluation")

  ### -------------------------------------------
  ### Evaluate models
  ### -------------------------------------------

  lis<-list()
  # loop over replicates
  for(i in 1:length(x@fits)){

    lisa<-list()
    # Loop over model types
    for(j in 1:length(x@fits[[1]])){

      # Make prediction
      pred=prd(x@fits[[i]][[j]],testa[[i]])

      #Feed with external threshold if available
      if(length(thres)==0){

        scores<-ceval(f=pred,
                      pa=papa[[i]],
                      tesdat=testa[[i]],
                      crit=crit)

      } else{

        scores<-ceval(f=pred,
                      pa=papa[[i]],
                      tesdat=testa[[i]],
                      tre=thres[which(names(thres)==names(x@fits[[i]])[j])],
                      crit=crit)

      }
      
      if(scores["threshold"]==Inf){
        scores["threshold"]=.5
      }
      lisa[[j]]<-scores

    }
    names(lisa)=names(x@fits[[i]])
    lis[[i]]<-lisa
  }
  names(lis)<-names(x@fits)
  out@performance<-lis

  return(out)
}
