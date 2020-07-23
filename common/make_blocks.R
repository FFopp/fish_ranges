### =========================================================================
### cv blocks
### =========================================================================
#' Block-wise split data into training and testing
#'
#' Creates a stratum vector based on a data.frame with n columns. If the data.frame
#' has one column strata are created based on clusters separated by quantiles. If
#' the data.frame has two or more columns, strata ere created based on k-medoid
#' clusters (function 'pam' from package cluster). Instead of a data.frame also the
#' argument 'npoints' can be provided, then groups are created by random sampling.
#' An opitimization algorithm (function 'gridSearch' from package NMOF) optimizes
#' for equal stratum sizes.
#'
#' @param nstrata number of approximately equal-sized classes to separate groups
#' in block-cross validation
#' @param df data.frame with n columns containing critera for cluster building.
#' Not necessary if argument npoints is supplied
#' @param nclusters number of clusters based on which strata should be built.
#' Minimum the same number as starta, maxuimum nrow(df)/10
#' @param npoints optional argument if 'df' is not supplied. For how many points
#' should random sampling be made?
#' @return Vector of length nrow(df) or npoints, with integers representing different strata
#' @author Philipp
#' @export
#' @examples
#' ### Test out block generation function
#'
#' # No layers supplied
#' strt.1=make_blocks(npoints=1000)
#' table(strt.1)
#'
#' # Stratified by 1d layer a
#' strt.2=make_blocks(df=env[,1,drop=F],nstrat=5,nclusters=5)
#' table(strt.2)
#'
#' # Stratified by 1d layer b
#' strt.3=make_blocks(df=env[,1,drop=F],nstrat=5,nclusters=15)
#' table(strt.3)
#'
#' # Stratified by 2d layer a
#' strt.4=make_blocks(df=env[,c(1,3)],nstrat=3,nclusters=3)
#' table(strt.4)
#'
#' # Stratified by 2d layer b
#' strt.5=make_blocks(df=env[,c(1,3)],nstrat=5,nclusters=50)
#' table(strt.5)
#'
#' # Stratified by 3d layer
#' strt.6=make_blocks(df=env[,1:3],nstrat=5,nclusters=50)
#' table(strt.6)
#'
#' par(mfrow=c(3,2))
#' plot(env[,c(1,3)],col=strt.1)
#' plot(env[,c(1,3)],col=strt.2)
#' plot(env[,c(1,3)],col=strt.3)
#' plot(env[,c(1,3)],col=strt.4)
#' plot(env[,c(1,3)],col=strt.5)
#' plot(env[,c(1,3)],col=strt.6)
#'
make_blocks<-function(nstrat=4,df=data.frame(),nclusters=nstrat*5,npoints=NA,pres=numeric()){

  ### ------------------------
  ### check input data
  ### ------------------------

  if(nrow(df)==0 & is.na(npoints)){
    stop("Please supply number of points if no data.frame is supplied")
  }

  ### ------------------------
  ### generate clusters
  ### ------------------------

  if(nrow(df)==0){

    ### do ordinary sampling if no strata are supplied
    out.strat=sample(rep(1:nstrat,ceiling(npoints/nstrat)),size=npoints)

  } else {

    # check for reasonable number of boxes
    if(nrow(df)<10*nclusters){
      stop("Too many boxes required!")
    }

    if(ncol(df)==1){
      ### do quantile-based clustering if df contains only one column
      rngi=quantile(df[,1],probs=0:(nclusters)/(nclusters))
      rngi[1]=rngi[1]-1
      rngi[length(rngi)]=rngi[length(rngi)]+1

      clist=as.numeric(cut(df[,1],breaks=rngi,right=TRUE))


    } else {

      # Scale input data
      scd=apply(df,2,scale)


      if(length(pres)==0){

        # do kmedoid clustering for 2 or more columns in df
        kmed=pam(scd,k=nclusters,metric="euclidean")

        # get clusters
        clist=kmed$clustering



      } else {

        # do kmedoid clustering for 2 or more columns in df
        kmed=pam(scd[which(pres==1),],k=nclusters,metric="manhattan")

        knnab=knn(train=scd[which(pres==1),],test=scd[which(pres==0),],cl=kmed$clustering)


        clist=kmed$clustering

        cliful=rep(NA,nrow(scd))
        cliful[which(pres==1)]=kmed$clustering
        cliful[which(pres==0)]=as.numeric(knnab)


      }
    }

    # sort obtained clusters
    tbl=sort(table(clist),decreasing = T)

    ### ------------------------
    ### regularly assign clusters to strata
    ### ------------------------

    if(nclusters != nstrat){

      # prepare strata layers
      grps=rep(list(numeric()),nstrat)

      # for the clusters with many observations
      # distribute them regularly among strata but keep
      # last six clusters for estimating most
      # regular distribution

      if(length(tbl)>6){

        for(i in 1:(length(tbl)-6)){

          # determine to which stratum the cluster should be
          # added
          fl=(floor((i-1)/nstrat))


          if(fl%%2==0){
            j=round(1+nstrat*((i-1)/nstrat-(floor((i-1)/nstrat))))
          } else {
            j=(nstrat+1)-round(1+nstrat*((i-1)/nstrat-(floor((i-1)/nstrat))))
          }
          # add cluster
          grps[[j]]=append(grps[[j]],tbl[i])

        }
      }

      # prepare for optimal distribution of last 6 clusters
      vlis=factor(1:nstrat,levels=1:nstrat)
      prs=rep(list(vlis),min(length(tbl),6))
      sstab=tbl[max(1,(length(tbl)-5)):length(tbl)]

      # Run brute-forcing gridSearch obtimization
      srch=gridSearch(levels=prs,fun=optme,nms=as.vector(sstab),grps=grps,tot=sum(tbl))

      # pull out results
      wi=as.numeric(as.character(srch$minlevels))

      # combine results with predistributed clusters
      for(i in 1:length(grps)){
        grps[[i]]=append(grps[[i]],sstab[wi==i])
      }

      # define vector with output strata
      out.strat=rep(NA,nrow(df))
      for(i in 1:length(grps)){
        if(length(pres)==0){
          out.strat[which(as.character(clist)%in%names(grps[[i]]))]=i
        } else{
          out.strat[which(as.character(cliful)%in%names(grps[[i]]))]=i
        }
      }

    } else {
      # if as many strata as clusters are required, simply return clusters
      if(length(pres)==0){
        out.strat=clist
      } else{
        out.strat=cliful
      }

    }

  }

  # return result
  return(out.strat)
}
