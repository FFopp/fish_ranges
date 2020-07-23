### =========================================================================
### upsample_thin function
### =========================================================================
#' Sample a supsample from a large set of points with a minimum distance
#' constraint
#'
#' Thin a spatial points object by number of points or minimum distance
#'
#' @param spdf SpatialPoints or SpatialPointsDataFrame object
#' @param lim_dist Minimum tolerated distance between the points
#' @param n_tot Minimum tolerated distance between the points
#' @details iteratively samples points and rejects them if min thist
#' to all existing points i the sample are not respected
#' @return SpatialPoints or SpatialPointsDataFrame with at least
#' lim_dist between all points
#' @author Philipp
#' @export
upsample_thin<-function(spdf,lim_dist,n_tot){

  proje=grepl("longlat",spdf@proj4string)

  smps=1:nrow(spdf@coords)
  strfy=spdf[sample(smps,1),]
  mxit=300000
  it=1
  while(nrow(strfy@coords)<n_tot){

    smpi=sample(smps,1)
    cand=spdf[smpi,]
    dsts=spDistsN1(strfy,cand,longlat = proje)
    if(all(dsts>lim_dist)){
      strfy=rbind(strfy,cand)
    }
    # Remove tested point from pool
    smps=smps[-which(smps==smpi)]
    it=it+1
    if(it==mxit){
      stop(paste0("upsample_thin failed, no success after",mxit,"trials!"))
    }
  }

  return(strfy)
}



