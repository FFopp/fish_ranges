### =========================================================================
### define wsl.glm function
### =========================================================================
#' Fit models
#'
#' Flexibly fit various types of functions but let framework take care of resampling,
#' meta-info storage, and file saving. Allows supplying any possible
#' model, however, functionality has been tested onlyn for GLM, GAM, GBM, randomForest,
#' and MaxEnt, and there may be problems with prediction/evaluation for exotic functions.
#'
#' @param pa vector with presence/absence values
#' @param env_vars data.frame with environmental predictors
#' @param taxon  name of the taxon for which models are fitted
#' @param replicatetype (how) should replicates be generated? may be 'none', 'splitsample',
#' 'cv' 'block-cv'
#' @param reps number of replicates
#' @param strata a numeric vector of the same length as observations with integers separating
#' cross validation replicates (used when replicatetype='block-cv')
#' @param save  should the model be saved in a structured way? (not implemented yet)
#' @param project character indicating the name of the project within which the models are run
#' (later used to define saving directories)
#' @param path where to save? (not implemented yet)
#' @param step (for glms and gams only) should the models be updated with the step function?
#' @param mod_tag (not in wsl.flex) label for the current model
#' @param mod_args list with elements of class 'multi.input' which specify models to be fitted
#' in wsl.flex
#' @return Object of class wsl.fit including slots for meta info, testing data for
#' evaluation, and model objects
#' @author Philipp
#' @name fitdoc
#' @examples
#'
#' # Take anguilla data set from dismo package
#' data("Anguilla_train")
#' vrs=c("SegSumT","USRainDays","USSlope")
#' env=Anguilla_train[,vrs]
#'
#' ### Check out wsl.flex
#' form.glm=as.formula(paste("Presence~",paste(paste0("poly(",vrs,",2)"),collapse="+")))
#' form.gam=as.formula(paste("Presence~",paste(paste0("s(",vrs,")"),collapse="+")))
#' form.tree=as.formula(Presence ~ .)
#' form.glm.2=as.formula(paste("Presence~",paste(vrs,collapse="+")))
#'
#' modinp=list(multi("glm",list(formula=form.glm,family="binomial"),"glm-simple",step=TRUE),
#' multi("gbm",list(formula=form.tree,
#'                  distribution = "bernoulli",
#'                  interaction.depth = 1,
#'                  shrinkage=.01,
#'                  n.trees = 3500),"gbm-simple"),
#' multi("gam",list(formula=form.gam,family="binomial"),"gam-simple",step=FALSE),
#' multi("randomForest",list(formula=form.tree,ntree=500,maxnodes=NULL),"waud1"),
#' multi("glm",list(formula=form.glm.2,family="binomial"),"glm-lin",step=TRUE))
#'
#' # Fit models
#' library(gam)
#' library(gbm)
#' library(randomForest)
#' library(MASS)
#' modi5=wsl.flex(pa=Anguilla_train$Angaus,
#'                env_vars = env,
#'                taxon="Angaus",
#'                replicatetype="block-cv",
#'                reps=3,
#'                strata=sample(1:3,nrow(env),replace=TRUE),
#'                project="multitest",
#'                mod_args=modinp)
#'
#' # Try out custom summary function
#' summary(modi5)
#'
#' # Access glm object of first replicate
#' summary(modi5@fits$replicate_01$`glm-simple`)
#'
#' # Evaluate the model
#' eval5<-wsl.evaluate(modi5,crit="pp=op")
#'
#' # Get evaluation summary
#' summary(eval5)
#'
#' # Get thresholds
#' thr.5=get_thres(eval5)
#'
#' ### Make some predictions
#' pred5a=wsl.predict(modi5,predat=env)
#' pred5b=wsl.predict(modi5,predat=env,thres=thr.5)
#'
NULL

#' @rdname fitdoc
#' @export
wsl.flex<-function(x=numeric(),
                   pa=numeric(),
                   env_vars=data.frame(),
                   taxon=character(),
                   replicatetype=character(),
                   reps,
                   strata=NA,
                   save=FALSE,
                   project=NA,
                   path=NA,
                   mod_args=list()){

  # Check supplied model types
  for(i in 1:length(mod_args)){
    if(!(mod_args[[i]]@mod%in%c("glm","gam","gbm","maxent","randomForest"))){
      warning(paste(mod_args[[i]]@mod,"not in focal model functions. You might run in to problems when evaluating/predicting..."))
    }
  }

  # Check if pseudo absence object is supplied
  if(class(x)=="wsl.pseudoabsences"){
    pa=x@pa
    env_vars=x@env_vars
    taxon=x@meta$taxon
  }

  # check and prepare data and output
  lis=preps(call=match.call())

  # loop over replicates
  fits=list()
  for(i in 1:reps){

    modi=list()
    # loop over models
    for(j in 1:length(mod_args)){

      if(mod_args[[j]]@mod=="maxent"){

        # Create directory for temporary MaxEnt data
        hde(mod_args[[j]]@args$me.temp.dir<-paste("tmp_Maxent",
                                                  print(as.numeric(Sys.time())*1000,digits=15),
                                                  sep="_"))
        dir.create(mod_args[[j]]@args$me.temp.dir)

        mod_args[[j]]@args$x<-lis$train[[i]][,-which(colnames(lis$train[[i]])=="Presence")]
        mod_args[[j]]@args$p<-lis$train[[i]][,"Presence"]

        hde(modi[[j]]<-do.call(mod_args[[j]]@mod,mod_args[[j]]@args))

        #Remove Temporary folder for Maxent
        unlink(mod_args[[j]]@args$me.temp.dir,recursive=T)

      } else {

        mod_args[[j]]@args$data=lis$train[[i]]

        # Make weight vector
        wi=which(mod_args[[j]]@args$data$Presence==1)
        wt=rep(1,nrow(mod_args[[j]]@args$data))
        wt[wi]<-round((nrow(mod_args[[j]]@args$data)-length(wi))/length(wi))

        if(mod_args[[j]]@weight){
          mod_args[[j]]@args$weights=wt
        }

        if(mod_args[[j]]@mod=="randomForest"){
          mod_args[[j]]@args$data$Presence=as.factor(mod_args[[j]]@args$data$Presence)
        }

        modi[[j]]=do.call(mod_args[[j]]@mod,mod_args[[j]]@args)
      }

      names(modi)[j]=ifelse(mod_args[[j]]@tag=="",paste0("model_",j),mod_args[[j]]@tag)

      if(mod_args[[j]]@step){
        modi[[j]]=stepAIC(modi[[j]],direction="both",trace=FALSE)
      }
    }

    fits[[i]]=modi

  }

  names(fits)=paste0("replicate_",sprintf("%02d",1:reps))

  # supply fitted objects
  lis$wslfi@fits=fits

  # Save
  #...

  return(lis$wslfi)

}
