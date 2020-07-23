### =========================================================================
### evaluate core function (peval)
### =========================================================================
#' Do the actual model evaluations
#'
#' Not to be called directly by the user
#'
#' @author Philipp
#' @export
ceval<-function(f,pa,tesdat,crit,tre=numeric()){

  # If there are any presences in the evaluation data
  if(any(pa==1)){

    z<-prediction(f,pa)
    # AUC
    auc<-performance(z,measure="auc")@y.values[[1]]
    rmse=performance(z,measure="rmse")@y.values[[1]]

    # optimum Threshold for conversion into binary outputs
    prbos<-seq(from=0,to=1,length.out=length(pa))
    zz<-performance(z,measure="sens",x.measure="spec",prbe=100000)
    zzz<-performance(z,measure="fpr",x.measure="fnr",prbe=100000)
    all.tss=zz@x.values[[1]]+zz@y.values[[1]]-1
    all.ppv<-performance(z,measure="ppv",prbe=100000)
    all.acc<-performance(z,measure="acc",prbe=100000)
    pn=zz@x.values[[1]]*zzz@x.values[[1]]
    py=zz@y.values[[1]]*zzz@y.values[[1]]

    all.kappa=(all.acc@y.values[[1]]-py*pn)/(1-py*pn)


    if(crit=="maxTSS"){
      thr=z@cutoffs[[1]][which.max(all.tss)]

    } else if(crit=="pp=op"){
      thr=quantile(f,probs=1-mean(pa))

    }else if(length(tre)!=0){
      thr=tre

    }

    wi=which.min(abs(thr-z@cutoffs[[1]]))
    ppv=all.ppv@y.values[[1]][wi]
    tss=all.tss[wi]
    acc=all.acc@y.values[[1]][wi]
    kappa=all.kappa[wi]

    # Return evaluation metrics
    weg=c(auc=auc,rmse=rmse,ppv=ppv,tss=tss,acc=acc,kappa=kappa,threshold=as.numeric(thr))

    return(weg)
  }
}
